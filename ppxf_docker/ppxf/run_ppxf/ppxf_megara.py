#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 11:54:02 2024

@author: mario
"""

import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from os import path
# from urllib import request

from ppxf.ppxf import ppxf, robust_sigma
import ppxf.ppxf_util as util
import ppxf.sps_util as lib
from vorbin.voronoi_2d_binning import voronoi_2d_binning
from plotbin.display_bins import display_bins

import argparse

parser = argparse.ArgumentParser(description='Script for analyzing MEGARA cube')
parser.add_argument('filename', type=str, help='Path to the spectra file to be analyzed')
parser.add_argument('--wave-range', nargs=2, type=float, metavar=('start', 'end'),
                    help='Range of wavelengths to focus on')
parser.add_argument('--target-sn', type=int, default=175, metavar='SN',
                    help='Target signal-to-noise ratio for Voronoi binning (default: 175)')
parser.add_argument('--sps-name', type=str, default='emiles', metavar='SPS_NAME',
                    help='Name of the SPS model to be used (default: emiles)')
parser.add_argument('--mask-file', type=str, help='Path to the mask file')
parser.add_argument('--output-dir', type=str, default='./', metavar='OUTPUT_DIR',
                    help='Directory where the output FITS files will be saved (default: current directory)')
parser.add_argument('--suffix', type=str, default='', help='Optional label to append to the output file names')
parser.add_argument('--redshift', type=float, default=0.0, metavar='REDSHIFT',
                    help='Redshift of the galaxy (default: 0)')
parser.add_argument('--velocity-dispersion', type=float, default=100.0, metavar='VEL_DISP',
                    help='Initial guess for the velocity dispersion of the galaxy (default: 100 km/s)')
parser.add_argument('--sn-range', nargs=2, type=float, metavar=('start', 'end'),
                    help='Range of wavelengths to calculate the signal-to-noise ratio (default: same as --wave-range)')
                    
args = parser.parse_args()

class read_megara_cube(object):
    def __init__(self, spectra_filename, wave_range):
        """
        Read MEGARA cube, log rebin it and compute coordinates of each spaxel.
        
        """
        filename = spectra_filename
        print("=" * 80)
        print(f"Analyzing file:", filename)
        print("=" * 80)
        hdu = fits.open(filename)
        head = hdu[0].header
        cube = hdu[0].data   # cube.shape = (4300, nx, ny)
        
        vph_parameters = {
            'LR-U': (0.672, 4051),
            'LR-B': (0.792, 4800),
            'LR-V': (0.937, 5695),
            'LR-R': (1.106, 6747),
            'LR-I': (1.308, 7991),
            'LR-Z': (1.455, 8900),
            'MR-U': (0.326, 4104),
            'MR-UB': (0.358, 4431),
            'MR-B': (0.395, 4814),
            'MR-G': (0.433, 5213),
            'MR-V': (0.476, 5667),
            'MR-VR': (0.522, 6170),
            'MR-R': (0.558, 6563),
            'MR-RI': (0.608, 7115),
            'MR-I': (0.666, 7767),
            'MR-Z': (0.796, 9262),
            'HR-R': (0.355, 6646),
            'HR-I': (0.462, 8634)
        }
        
        try:
            red = head['VPH']
            print("=" * 80)
            print("Observation VPH:", red)
            print("=" * 80)
            FWHM_VPH, lambda_c = vph_parameters.get(red, (0.937, 5695))  # Default values for missing 'VPH' keys
        except KeyError:
            print("=" * 80)
            print("No VPH information in the header. Using LR-V as default FWHM and central wavelength.")
            print("=" * 80)
            FWHM_VPH, lambda_c = 0.937, 5695
                        
        # Transform cube into 2-dim array of spectra
        npix = cube.shape[0]
        spectra_Jy = cube.reshape(npix, -1) # create array of spectra [npix, nx*ny]
        wave = head['CRVAL3'] + head['CDELT3']*np.arange(npix)
        pixsize = abs(head["CDELT1"])*3600    # 0.2"

        # Only use a restricted wavelength range
        w = (wave > wave_range[0]) & (wave < wave_range[1])
        spectra_Jy = spectra_Jy[w, :]
        wave = wave[w]
        wave_2D = np.tile(wave, (spectra_Jy.shape[1], 1)).T
        spectra_cgs = spectra_Jy*1e20 / (3.33564095e4*(wave_2D**2))   ## Units 10**-20 erg cm-2 s-1 Angstrom-1
        spectra_cgs = np.maximum(spectra_cgs, 0)
        
        # Create coordinates centred on the brightest spectrum
        flux = np.nanmean(spectra_cgs, 0)
        jm = np.argmax(flux)
        row, col = map(np.ravel, np.indices(cube.shape[-2:]))
        x = (col - col[jm])*pixsize
        y = (row - row[jm])*pixsize
        c = 299792.458  # speed of light in km/s
        velscale = c*np.diff(np.log(wave[-2:]))  # Smallest velocity step
        lam_range_temp = [np.min(wave), np.max(wave)]
        spectra_cgs, ln_lam_gal, velscale = util.log_rebin(lam_range_temp, spectra_cgs, velscale=velscale)
        # spectra, ln_lam_gal, velscale = util.log_rebin(lam_range_temp, spectra, velscale=velscale)

        self.original_data = cube
        self.header = head
        self.spectra = spectra_cgs
        self.x = x
        self.y = y
        self.col = col + 1   # start counting from 1
        self.row = row + 1
        self.flux = flux
        self.ln_lam_gal = ln_lam_gal
        self.lam_gal = wave
        self.fwhm_gal = FWHM_VPH
        self.lambda_c = lambda_c
        self.cube_shape = cube.shape

## Function to iteratively clip the outliers

def clip_outliers(galaxy, bestfit, goodpixels):
    """
    Repeat the fit after clipping bins deviants more than 3*sigma
    in relative error until the bad bins don't change any more.
    """
    while True:
        scale = galaxy[goodpixels] @ bestfit[goodpixels]/np.sum(bestfit[goodpixels]**2)
        resid = scale*bestfit[goodpixels] - galaxy[goodpixels]
        err = robust_sigma(resid, zero=1)
        ok_old = goodpixels
        goodpixels = np.flatnonzero(np.abs(bestfit - galaxy) < 3*err)
        if np.array_equal(goodpixels, ok_old):
            break
            
    return goodpixels


## Function to fit the stellar kinematics
#The following function fits the spectrum with `pPXF` while masking the gas emission lines, then iteratively clips the outliers and finally refit the spectrum with `pPXF` on the cleaned spectrum.

def fit_and_clean(templates, galaxy, velscale, start, goodpixels0, lam, lam_temp, plot_title):
    
    print('##############################################################')
    goodpixels = goodpixels0.copy()
    pp = ppxf(templates, galaxy, np.ones_like(galaxy), velscale, start,
              moments=4, degree=4, mdegree=4, lam=lam, lam_temp=lam_temp,
              goodpixels=goodpixels, plot_title=plot_title, reddening=2)
    
    #plt.figure(figsize=(8, 3.7))
    #plt.subplot(121)
    #pp.plot_mine()

    goodpixels = clip_outliers(galaxy, pp.bestfit, goodpixels)

    # Add clipped pixels to the original masked emission lines regions and repeat the fit
    goodpixels = np.intersect1d(goodpixels, goodpixels0)
    pp = ppxf(templates, galaxy, np.ones_like(galaxy), velscale, start,
              moments=4, degree=4, mdegree=4, lam=lam, lam_temp=lam_temp,
              goodpixels=goodpixels, plot_title=plot_title, reddening=2)
    
    #plt.subplot(122)
    pp.plot_mine()


    optimal_template = templates @ pp.weights
    
    return pp, optimal_template

def create_mask(lam, exclude_ranges, apply_redshift, z_galaxy):
    """
    Create a boolean mask to exclude specific wavelength ranges, optionally
    considering the redshift of the galaxy.

    Parameters:
    -----------
    lam : array_like
        Array of observed wavelengths.

    exclude_ranges : list of tuples
        List of tuples specifying the rest-frame wavelength ranges to exclude.
        Each tuple should contain the start and end wavelengths of the range.

    apply_redshift : list of bool
        List of boolean values indicating whether to apply the redshift to
        each exclusion range.

    z_galaxy : float
        Redshift of the galaxy.

    Returns:
    --------
    mask : array_like
        Boolean mask with the same shape as `lam`, where `True` indicates
        wavelengths to be excluded.
    """
    mask = np.ones_like(lam, dtype=bool)
    for i, (start, end) in enumerate(exclude_ranges):
        observed_start = start * (1 + z_galaxy) if apply_redshift[i] else start
        observed_end = end * (1 + z_galaxy) if apply_redshift[i] else end
        mask &= (lam < observed_start) | (lam > observed_end)
    return mask

def read_mask_info(mask_file):
    """
    Reads mask information from a text file.

    Parameters:
    -----------
    mask_file : str
        Name of the file containing the mask information.

    Returns:
    --------
    exclude_ranges : list of tuples
        List of tuples specifying the wavelength ranges to exclude.

    apply_redshift : list of bool
        List of boolean values indicating whether to apply redshift or not to each range.
    """
    exclude_ranges = []
    apply_redshift = []
    with open(mask_file, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                start, end, apply_z = parts[0], parts[1], parts[2]
                exclude_ranges.append((float(start), float(end)))
                apply_redshift.append(apply_z.lower() == 'true')
    return exclude_ranges, apply_redshift


def save_results_as_fits(result_data, original_file, new_filename):
    """
    Saves the result data as a new FITS file with the header extracted from the original file.

    Parameters:
    result_data (array_like): The data to be saved in the FITS file.
    original_file (str): The path to the original FITS file from which to extract the header.
    new_filename (str): The desired filename for the new FITS file.

    Returns:
    None
    """
    # Read the header from the original file
    with fits.open(original_file) as hdul:
        header = hdul[0].header

    # Modify the keyword with the new value
    header["CRVAL3"] = 1
    header["CDELT3"] = 1
    
    # Create a new FITS extension with the result data and the original header
    new_hdu = fits.PrimaryHDU(data=result_data, header=header)

    # Save the FITS file
    new_hdu.writeto(new_filename, overwrite=True)


def calculate_signal_to_noise(spectrum, wavelength, mask_file=None, wavelength_range=None, redshift=None):
    """
    Calculate the signal-to-noise ratio (SNR) of a spectrum within a specified wavelength range.

    Args:
        spectrum (numpy.ndarray): The spectrum data.
        wavelength (numpy.ndarray): The wavelength values corresponding to the spectrum.
        mask_file (str, optional): The path to the mask file containing regions to be excluded.
        wavelength_range (tuple, optional): Tuple specifying the start and end wavelengths of the range.
        redshift (float, optional): Redshift of the galaxy.

    Returns:
        tuple: A tuple containing the signal (np.median(column_spectrum[valid_indices])), and noise (np.sqrt(column_signal)) for each column of the spectrum.
    """
    # Calculate the factor to convert from S/N per pixel to S/N per Angstrom
    factor_aa = 1 / np.sqrt(np.mean(np.diff(wavelength)))

    # Make a copy of the spectrum to avoid modifying the original array
    spectrum_copy = spectrum.copy()

    # Initialize arrays to store signal and noise for each column
    num_columns = spectrum_copy.shape[1]
    signal = np.zeros(num_columns)
    noise = np.zeros(num_columns)

    # Apply mask from mask file if provided
    if mask_file:
        exclude_ranges, apply_redshift = read_mask_info(mask_file)
        mask = create_mask(wavelength, exclude_ranges, apply_redshift, redshift)
        spectrum_copy[mask] = np.nan

    # Apply wavelength range if provided
    if wavelength_range:
        start, end = wavelength_range
        start_index = np.argmax(wavelength >= start)
        end_index = np.argmin(wavelength <= end)
        spectrum_copy[start_index:end_index+1, :] = np.nan

    # Calculate signal and noise for each column
    for i in range(num_columns):
        column_spectrum = spectrum_copy[:, i]
        valid_indices = np.isfinite(column_spectrum)
        column_signal = np.median(column_spectrum[valid_indices])
        column_noise = np.sqrt(column_signal)
        signal[i] = column_signal
        noise[i] = column_noise

    # Calculate the signal-to-noise ratio for each column
    snr = signal / noise

    return signal, noise

def get_age_metal_weights(self, weights):
    """
    Get the age, metallicity, and their corresponding weights.

    Args:
        weights (numpy.ndarray): Array of weights returned by pPXF.

    Returns:
        tuple: A tuple containing three arrays:
            - age_values: Array of age values.
            - metal_values: Array of metallicity values.
            - weight_values: Array of corresponding weights.
    """
    # Initialize lists to store values
    age_values = []
    metal_values = []
    weight_values = []

    # Iterate over each grid point
    for i in range(self.age_grid.shape[0]):
        for j in range(self.age_grid.shape[1]):
            # Get age, metallicity, and weight at grid point (i, j)
            age = self.age_grid[i, j]
            metallicity = self.metal_grid[i, j]
            weight = weights[i, j]

            # Append values to lists
            age_values.append(age)
            metal_values.append(metallicity)
            weight_values.append(weight)

    # Convert lists to numpy arrays
    age_values = np.array(age_values)
    metal_values = np.array(metal_values)
    weight_values = np.array(weight_values)

    return age_values, metal_values, weight_values

def correct_sigma(fwhm_gal, fwhm_template, sigma_ppxf, central_wavelength):
    """
    Correct the sigma values when instrumental dispersion is smaller than
    the template dispersion.

    Args:
        fwhm_gal (float): FWHM of the galaxy in angstrom.
        fwhm_template (float): FWHM of the template in angstrom.
        sigma_ppxf (float): pPXF measured dispersion.
        central_wavelength (float): Central wavelength in angstrom.

    Returns:
        float: Corrected sigma value.
    """
    # Convert FWHM to sigma in km/s
    c = 299792.458  # Speed of light in km/s
    
    sigma_gal = (fwhm_gal / (2 * np.sqrt(2 * np.log(2)))) * (c / central_wavelength)
    sigma_template = (fwhm_template / (2 * np.sqrt(2 * np.log(2)))) * (c / central_wavelength)
    
    # Calculate the difference in sigma
#    sigma_diff_sq = np.sqrt(sigma_gal**2 - sigma_template**2)
    
    # Correct sigma
#    corrected_sigma = np.sqrt(sigma_ppxf**2 - sigma_diff_sq)
    corrected_sigma = np.sqrt(sigma_ppxf**2 - (sigma_gal**2 - sigma_template**2))
    return corrected_sigma

def generate_voronoi_cubes(data_cube, voronoi_bins):
    # Get the dimensions of the data cube
    nz, nx, ny = data_cube.shape
    
    # Create an empty dictionary to store the summed spectra for each Voronoi cell
    voronoi_cells = {}
    # Create an empty dictionary to count the number of spaxels in each Voronoi cell
    voronoi_counts = {}
    
    # Iterate over each spaxel in the Voronoi data
    for x, y, cell in voronoi_bins:
        # If the cell is not in the dictionary, add it with an empty spectrum and count
        if cell not in voronoi_cells:
            voronoi_cells[cell] = np.zeros(nz)
            voronoi_counts[cell] = 0
        
        # Add the spectrum of the current spaxel to the corresponding Voronoi cell
        voronoi_cells[cell] += data_cube[:, x, y]
        # Increment the count of spaxels in the Voronoi cell
        voronoi_counts[cell] += 1
    
    # Create an empty data cube for the Voronoi cells
    voronoi_cube = np.zeros((nz, nx, ny))
    
    # Assign the averaged spectra to the Voronoi cube
    for x, y, cell in voronoi_bins:
        voronoi_cube[:, x, y] = voronoi_cells[cell] / voronoi_counts[cell]
    
    return voronoi_cube

# Read header to get wave range if not provided
if args.wave_range is None:
    with fits.open(args.filename) as hdul:
        header = hdul[0].header
        if 'WAVLIMF1' in header and 'WAVLIMF2' in header:
            args.wave_range = [header['WAVLIMF1'], header['WAVLIMF2']]
        else:
            # Use full wavelength range of observation
            # Extract wavelength range from header
            crval3 = header['CRVAL3']  # starting wavelength
            cdelt3 = header['CDELT3']  # wavelength step
            npix = header['NAXIS3']    # number of pixels in the wavelength direction
            args.wave_range = [crval3, crval3 + (npix - 1) * cdelt3]

print("=" * 80)
print(f"Using wavelength range:", args.wave_range)
print("=" * 80)

## Read the data cube and Voronoi bin the data
#I only extract the cube over the wavelength region where there are emission lines and where the spectrum is less contaminated by sky residuals.

lam_range_temp = args.wave_range
if args.sn_range is None:
    args.sn_range = lam_range_temp
print("=" * 80)
print("Using wavelength range for signal-to-noise calculation:", args.sn_range)
print("=" * 80)

spectra_filename = args.filename  # Choose the spectra to be fitted
s = read_megara_cube(spectra_filename, lam_range_temp)

# In this example I request an excessively large `target_sn=350` to speed up the calculation. This generates only 9 Voronoi bins. But in a real situation the spatially binned data cube will contain over a hundred Voronoi bins e.g. with `target_sn=60`.

#signal = np.median(s.spectra, 0)
#noise = np.sqrt(signal)
# target_sn = 350

signal, noise = calculate_signal_to_noise(s.spectra,np.exp(s.ln_lam_gal), mask_file=args.mask_file, wavelength_range=args.sn_range, redshift=args.redshift)
target_sn = args.target_sn


# Saving files nomenclature

base_filename = os.path.basename(args.filename).replace('.fits', '')
suffix = f"_{args.suffix}" if args.suffix else ''

# Perform Voronoi binning with the method of [Cappellari & Copin (2003)](https://ui.adsabs.harvard.edu/abs/2003MNRAS.342..345C)

plt.figure(figsize=(7,10))
bin_num, x_gen, y_gen, xbin, ybin, sn, nPixels, scale = voronoi_2d_binning(s.x, s.y, signal, noise, target_sn, plot=1, quiet=1)
plt.savefig(path.join(args.output_dir, args.filename.split("/")[-1].replace(".fits", f"_voronoi_binning_sn_{target_sn}{suffix}.pdf")),dpi=600)
# Saving the information of which Voronoi region each pixel belongs to

# Create index matrices
indices_x, indices_y = np.meshgrid(np.arange(s.cube_shape[1]), np.arange(s.cube_shape[2]), indexing='ij')
# Flatten the index matrices
x_flat = indices_x.flatten()
y_flat = indices_y.flatten()
voronoi_bins = np.column_stack((x_flat, y_flat, bin_num))
np.savetxt(path.join(args.output_dir, args.filename.split("/")[-1].replace(".fits", f"_voronoi_binning_info_sn_{target_sn}{suffix}.txt")), voronoi_bins, fmt='%.i %.i %.i', comments='')
# Generate the Voronoi cubes
voronoi_cube_data = generate_voronoi_cubes(s.original_data, voronoi_bins)
# Save the Voronoi cube to a new FITS file with the same header as the original cube
hdu_voronoi = fits.PrimaryHDU(data=voronoi_cube_data, header=s.header)
voronoi_results_FITS = os.path.join(args.output_dir, f"{base_filename}_voronoi_binned_sn_{target_sn}{suffix}.fits")
hdu_voronoi.writeto(voronoi_results_FITS, overwrite=True)
# Save SNR maps to FITS files
snr = (signal/noise).reshape((s.cube_shape[1], s.cube_shape[2]))
hdu_snr_before = fits.PrimaryHDU(data=snr, header=header)
snr_results_FITS = os.path.join(args.output_dir, f"{base_filename}_snr_map{suffix}.fits")
hdu_snr_before.writeto(snr_results_FITS, overwrite=True)

## Setup stellar templates
#The important formula below **defines** the relation between velocity, wavelength and redshift in ``pPXF`` (eq. 8 of [Cappellari 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.466..798C))

#$$V\equiv c\Delta\ln\lambda = c\ln(1+z)$$

c_kms = 299792.458  # speed of light in km/s
velscale = c_kms*np.diff(s.ln_lam_gal[:2])   # eq.(8) of Cappellari (2017)
velscale = velscale[0]   # must be a scalar

# pPXF can be used with any set of SPS population templates. However, I am currently providing (with permission) ready-to-use template files for three SPS. One can just uncomment one of the three models below. The included files are only a subset of the SPS that can be produced with the models, and one should use the relevant software to produce different sets of SPS templates if needed.

# 1. If you use the [fsps v3.2](https://github.com/cconroy20/fsps) SPS model templates, please also cite [Conroy et al. (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJ...699..486C) and [Conroy et al. (2010)](https://ui.adsabs.harvard.edu/abs/2010ApJ...712..833C) in your paper.

# 2. If you use the [GALAXEV v2000](http://www.bruzual.org/bc03/) SPS model templates, please also cite [Bruzual & Charlot (2003)](https://ui.adsabs.harvard.edu/abs/2003MNRAS.344.1000B) in your paper.

# 3. If you use the [E-MILES](http://miles.iac.es/) SPS model templates, please also cite [Vazdekis et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.463.3409V) in your paper. 
# <font color='red'>WARNING: the E-MILES models do not include very young SPS and should not be used for highly star forming galaxies.</font>  


# sps_name = 'fsps'
# sps_name = 'galaxev'
# sps_name = 'emiles'
sps_name = args.sps_name

# Read SPS models file from my GitHub if not already in the pPXF package dir. I am not distributing the templates with pPXF anymore.
# The SPS model files are also available [this GitHub page](https://github.com/micappe/ppxf_data).

ppxf_dir = path.dirname(path.realpath(lib.__file__))
basename = f"spectra_{sps_name}_9.0.npz"
filename = path.join(ppxf_dir, 'sps_models', basename)

# Uncomment these lines to give the option of searching for the SSP model files on GitHub.

if not path.isfile(filename):
    url = "https://raw.githubusercontent.com/micappe/ppxf_data/main/" + basename
    request.urlretrieve(url, filename)
    
FWHM_gal = s.fwhm_gal   # set this to None to skip convolution
sps = lib.sps_lib(filename, velscale, FWHM_gal, norm_range=[5300, 5950])
stars_templates, ln_lam_temp = sps.templates, sps.ln_lam_temp

# The stellar templates are reshaped into a 2-dim array with each spectrum as a column, however we save the original array dimensions, which are needed to specify the regularization dimensions

reg_dim = stars_templates.shape[1:]
stars_templates = stars_templates.reshape(stars_templates.shape[0], -1)

# See the pPXF documentation for the keyword REGUL, for an explanation of the following two lines.

stars_templates /= np.median(stars_templates) # Normalizes stellar templates by a scalar
regul_err = 0.01 # Desired regularization error

z = args.redshift  # redshift estimate from NED
print("=" * 80)
print("Galaxy redshift:", z)
print("=" * 80)
vel0 = c_kms*np.log(1 + z)  # Initial estimate of the galaxy velocity in km/s. eq. (8) of Cappellari (2017)
sigma = args.velocity_dispersion
print("Velocity dispersion initial guess:", sigma, "km/s")
print("=" * 80)
start = [vel0, sigma]  # (km/s), starting guess for [V,sigma]

lam_range_temp = np.exp(ln_lam_temp[[0, -1]])
goodpixels0 = util.determine_goodpixels(s.ln_lam_gal, lam_range_temp, z, width=1000)



## Fit templates and stellar kinematics in Voronoi binned data


nbins = sn.size
# velbin, velbin_error, sigbin, sigbin_error, lg_age_bin, metalbin, nspax = np.zeros((7, nbins))
velbin, velbin_error, sigbin, sigbin_error, h3bin, h3bin_error, h4bin, h4bin_error, lg_age_bin, metalbin, nspax, attbin = np.zeros((12, nbins)) # Uncomment for fitting h3 and h4 and comment the previous line
optimal_templates = np.empty((stars_templates.shape[0], nbins))
# Initialize a cube to store the residuals with the same spatial dimensions as the original data cube
residuals_cube = np.zeros((len(s.ln_lam_gal), s.cube_shape[1], s.cube_shape[2]))
# Initialize a list to store residuals for each Voronoi bin
residuals_list = []

lam_gal = np.exp(s.ln_lam_gal)

# Read mask information from the text file if mask_file argument is provided
if args.mask_file:
    exclude_ranges, apply_redshift = read_mask_info(args.mask_file)
    mask = create_mask(lam_gal, exclude_ranges, apply_redshift, z)
    goodpixels_masked = goodpixels0[mask]
else:
    goodpixels_masked = goodpixels0

# Directory to save the figures
# plot_dir = path.join(args.output_dir, "Plots_spectral_fitting")
plot_dir = os.path.join(args.output_dir, f"Plots_spectral_fitting_{base_filename}_sn_{target_sn}{suffix}")

# Check if the directory exists, and if not, create it
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

# Save the information about stellar pops SSP used for the fitting of each Voronoi bin
voronoi_output_file = path.join(args.output_dir, args.filename.split("/")[-1].replace(".fits", f"_stellar_pops_info_sn_{target_sn}{suffix}.txt"))
voronoi_output_file = open(voronoi_output_file, "w")  # Open file in write mode
voronoi_output_file.write("Age Metallicity Weight\n")

for j in range(nbins):
    w = bin_num == j
    galaxy = np.sum(s.spectra[:, w], 1)
    
    plot_title = str('Bin ' + str(j+1))
    
    pp, bestfit_template = fit_and_clean(stars_templates, galaxy, velscale, start, goodpixels_masked, lam_gal, sps.lam_temp, plot_title)
    # velbin[j], sigbin[j] = pp.sol
    velbin[j], sigbin[j], h3bin[j], h4bin[j] = pp.sol  # Uncomment for fitting h3 and h4 and comment the previous line
    velbin_error[j], sigbin_error[j], h3bin_error[j], h4bin_error[j] = pp.error*np.sqrt(pp.chi2)  # Uncomment for fitting h3 and h4 and comment the previous line
    optimal_templates[:, j] = bestfit_template
    attbin[j] = pp.reddening
    # Calculate the residuals for the current bin
    residuals = galaxy - pp.bestfit
    residuals_list.append(residuals)
    
    # Save figures of the fittings
    
    # Code to generate the figure
    plt.savefig(os.path.join(plot_dir, ('Voronoi_bin_' + str(j+1) + '.pdf')),dpi=600)
    plt.close()
    
    light_weights = pp.weights.reshape(reg_dim)
    lg_age_bin[j], metalbin[j] = sps.mean_age_metal(light_weights)
    
    # Save the SSP templates used in the fitting together with their respectives weights (weight_values_nonzero/np.sum(weight_values))
    
    age_values, metal_values, weight_values = get_age_metal_weights(sps, light_weights)
    
    # Get indices of values with nonzero weight
    nonzero_indices = weight_values > 0

    # Filter age, metallicity, and weight values
    age_values_nonzero = age_values[nonzero_indices]
    metal_values_nonzero = metal_values[nonzero_indices]
    weight_values_nonzero = weight_values[nonzero_indices]/np.sum(weight_values) # Normalized weights
    
    # Combine age, metallicity, and weight arrays into a single 2D array
    sp_data_info = np.column_stack((age_values_nonzero, metal_values_nonzero, weight_values_nonzero))
    
    # Write header
    voronoi_output_file.write("\n" + "Voronoi bin #" + str(j) + "\n")

    # Write data
    for age, metal, weight in sp_data_info:
        voronoi_output_file.write(f"{age:.8f}\t{metal:.8f}\t{weight:.8f}\n")
    
    print("-"*50)
    print("Age values [Gyr]",age_values_nonzero)
    print("Metallicity values [M/H]",metal_values_nonzero)
    print("Weights [Normalized]",weight_values_nonzero)
    print("-"*50)
    
    print(f'Voronoi bin: {j + 1} / {nbins}')
#    plt.title(f"Voronoi bin {j + 1} / {nbins}")
    
    
# Correction of sigma in case the instrumental dispersion is smaller than the dispersion of the SSP models
if FWHM_gal < sps.fwhm_tem[0]:
    corrected_sigbin = correct_sigma(FWHM_gal, sps.fwhm_tem[0], sigbin, s.lambda_c)
    corrected_sigbin_error = sigbin*sigbin_error/corrected_sigbin
else:
    corrected_sigbin = sigbin
    corrected_sigbin_error = sigbin_error


# Assign the residuals to the corresponding spaxels in the residuals cube
for x, y, vor_bin_num in voronoi_bins:
    residuals_cube[:, x, y] = residuals_list[vor_bin_num]


# Save results as FITS files.

# Specify the output file names

# base_filename = os.path.basename(args.filename).replace('.fits', '')
# suffix = f"_{args.suffix}" if args.suffix else ''
kinematics_results_FITS = os.path.join(args.output_dir, f"{base_filename}_kinematics_sn_{target_sn}{suffix}.fits")
stellar_pops_results_FITS = os.path.join(args.output_dir, f"{base_filename}_stellar_pops_sn_{target_sn}{suffix}.fits")
attenuation_results_FITS = os.path.join(args.output_dir, f"{base_filename}_attenuation_sn_{target_sn}{suffix}.fits")
residuals_results_FITS = os.path.join(args.output_dir, f"{base_filename}_residuals_sn_{target_sn}{suffix}.fits")

#kinematics_results_FITS = path.join(args.output_dir, args.filename.split("/")[-1].replace(".fits", "_kinematics.fits"))
#stellar_pops_results_FITS = path.join(args.output_dir, args.filename.split("/")[-1].replace(".fits", "_stellar_pops.fits"))
#attenuation_results_FITS = os.path.join(args.output_dir, args.filename.split("/")[-1].replace(".fits", "_attenuation.fits"))

kinematics_fitting_results = np.concatenate([velbin[bin_num].reshape(-1,s.cube_shape[1],s.cube_shape[2]),
                                             velbin_error[bin_num].reshape(-1,s.cube_shape[1],s.cube_shape[2]),
                                             corrected_sigbin[bin_num].reshape(-1,s.cube_shape[1],s.cube_shape[2]),
                                             corrected_sigbin_error[bin_num].reshape(-1,s.cube_shape[1],s.cube_shape[2]),
                                             h3bin[bin_num].reshape(-1,s.cube_shape[1],s.cube_shape[2]),
                                             h3bin_error[bin_num].reshape(-1,s.cube_shape[1],s.cube_shape[2]),
                                             h4bin[bin_num].reshape(-1,s.cube_shape[1],s.cube_shape[2]),
                                             h4bin_error[bin_num].reshape(-1,s.cube_shape[1],s.cube_shape[2])], axis=0)

stellar_pops_fitting_results = np.concatenate([lg_age_bin[bin_num].reshape(-1,s.cube_shape[1],s.cube_shape[2]),
                                               metalbin[bin_num].reshape(-1,s.cube_shape[1],s.cube_shape[2])], axis=0)

attenuation_fitting_results = attbin[bin_num].reshape(-1,s.cube_shape[1],s.cube_shape[2])

# Call the function to save the results as a FITS file
save_results_as_fits(kinematics_fitting_results, spectra_filename, kinematics_results_FITS)
save_results_as_fits(stellar_pops_fitting_results, spectra_filename, stellar_pops_results_FITS)
save_results_as_fits(attenuation_fitting_results, spectra_filename, attenuation_results_FITS)
save_results_as_fits(residuals_cube, spectra_filename, residuals_results_FITS)
