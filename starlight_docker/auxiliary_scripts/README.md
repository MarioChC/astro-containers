# Auxiliary Scripts

This folder contains two auxiliary scripts that provide additional functionality for the main project.

## Script 1: convert_3Dfits_to_starlight.py

This script converts 3D FITS files to a format compatible with Starlight.

### Usage

```

python convert_3Dfits_to_starlight.py [options] [-h] [-s 3D DATACUBE] [-o OUTPUT DIRECTORY] [-i OUTPUT FILE]

Arguments:
-h, --help: Show help message and exit.
-s 3D DATACUBE, --spectrum 3D DATACUBE: Path to the 3D datacube.
-k KINEMATICS FITS, --kinematics KINEMATICS FITS: Path to the FITS file with kinematic information.
-o OUTPUT DIRECTORY, --output-directory OUTPUT DIRECTORY: Directory to save the output files.
-i OUTPUT FILE, --output-file OUTPUT FILE: Label added to saved output files.

Example:

python convert_3Dfits_to_starlight.py -s datacube.fits -k kinematics.fits -o output_dir -i output_file
```
### Output files

```
1. Spectral Files: Individual files for each spaxel in the format spectrum_xpos_<x>_ypos_<y>_fiber_number_<n>_<output_file>.txt, where:

<x> and <y> are the spaxel coordinates.
<n> is the fiber number.
<output_file> is the label provided by the -i option.

2. Kinematic Information File: A single file named información_propiedades_cinemáticas_<output_file>.txt that contains the kinematic properties for each spectrum. The columns are:

Spectral File: The name of the corresponding spectral file.
Velocity: The velocity value.
Velocity Error: The error in the velocity value.
Velocity Dispersion: The velocity dispersion.
Velocity Dispersion Error: The error in the velocity dispersion.
h3: The h3 value.
h3 Error: The error in the h3 value.
h4: The h4 value.
h4 Error: The error in the h4 value.
Each line in this file corresponds to one spectrum and its associated kinematic information.
```

## Script 2: starlight_grid_file_assign_spectra.py

This script generates a Starlight input file based on provided parameters and kinematic information.

### Usage

```
starlight_grid_file_assign_spectra.py [-h] output_file [spectrum_names [spectrum_names ...]] [kinematic_file]

Arguments:
output_file: Path to the output file (grid.in).
spectrum_names: Names of spectrum files to include in the analysis.
kinematic_file: Path to the file that contains the information about the kinematics.
-h, --help: Show help message and exit.

Example:

python starlight_grid_file_assign_spectra.py grid.in spectrum1.txt spectrum2.txt spectrum3.txt kinematic_info.txt

python starlight_grid_file_assign_spectra.py grid.in spectra_folder/* kinematic_info.txt
