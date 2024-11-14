# pPXF Cube Analysis Script

This script is designed to analyze MEGARA cube data. It offers various options for customization, allowing users to focus on specific wavelength ranges, adjust signal-to-noise ratios, select spectral population synthesis (SPS) models, and more.


## Usage

The bash_script.sh file contains the commands necessary to run the ppxf_megara.py script with the desired arguments. It should be copied to the shared directory of the container (/home/ppxf/shared_directory/input/) before execution.

An example command to run the script is as follows:

```
python3 /home/ppxf/run_ppxf/ppxf_megara.py /home/ppxf/shared_directory/input/NGC7025_LR-V_final_cube.fits --mask-file /home/ppxf/shared_directory/input/mask.txt --target-sn 200 --redshift 0.016571 --velocity-dispersion 200 --sn-range 5600 5800 --output-dir /home/ppxf/shared_directory/output
```

This command executes the ppxf_megara.py script with the specified input file (NGC7025_LR-V_final_cube.fits), mask file (mask.txt), target signal-to-noise ratio (200), redshift (0.016571), initial guess for velocity dispersion (200), signal-to-noise ratio range (5600 to 5800), and output directory (/home/ppxf/shared_directory/output).

## Input Arguments

filename: Path to the spectra file to be analyzed.

--wave-range start end: Range of wavelengths to focus on (default: full spectral range covered by the observation).

--target-sn SN: Target signal-to-noise ratio for Voronoi binning (default: 175).

--sps-name SPS_NAME: Name of the SPS model to be used (default: emiles).

--mask-file MASK_FILE: Path to the mask file (default: none).

--output-dir OUTPUT_DIR: Directory where the output FITS files will be saved (default: current directory).

--redshift REDSHIFT: Redshift of the galaxy (default: 0).

--velocity-dispersion VEL_DISP`: Initial guess for the velocity dispersion of the galaxy (default: 100 km/s).

--sn-range start end: Range of wavelengths to calculate the signal-to-noise ratio (default: same as --wave-range).

--min-sn MIN_SN: Minimum signal-to-noise ratio for spaxels to be included in the analysis (default: 0).

These arguments allow users to customize the analysis of MEGARA cube data by specifying parameters such as the file paths, wavelength ranges, signal-to-noise ratios, SPS models, output directories, and redshift values.

Feel free to adjust and expand upon this template to better suit your needs!
