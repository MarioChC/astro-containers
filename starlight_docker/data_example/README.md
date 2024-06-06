# Starlight Data Preparation

This folder contains data files necessary for running Starlight analysis.

## NGC7025_LR-V_final_cube.fits

This file (`NGC6027_LR-V_final_cube.fits`) is an observation of the galaxy NGC7025 in a 3D datacube format stored in a FITS file. It represents the spatial and spectral information of the galaxy.

## NGC7025_LR-V_final_cube_kinematics.fits

This file (NGC7025_LR-V_final_cube_kinematics.fits) is an output provided by the pPXF (Penalized Pixel-Fitting) program. It contains the kinematic analysis of the galaxy NGC7025, including information such as velocity, velocity dispersion, skewness (h3) and kurtosis (h4).

## kinematic_information_file_NGC7025_LR-V.txt

This file (located in the "starlight_format" directory) includes the kinematic information from the file NGC7025_LR-V_final_cube_kinematics.fits in a text format and relates it to the nomenclature of the data in the format compatible with Starlight.

## starlight_format

This subfolder contains individual spectra extracted from the data cube (`NGC7025_LR-V_final_cube.fits`). These spectra are stored in text files (`*.txt`) and are formatted for use with Starlight software. Each text file contains the spectrum corresponding to a specific position (x, y) in the data cube, and its filename indicates the coordinates (xpos, ypos) from which the spectrum was extracted.
