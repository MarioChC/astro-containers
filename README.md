This project packages the starlight project in a container,
this is to aid with running this project in modern
systems. This works with fortran77 which is more than 12 years
old, still the software is useful and commonly used in the astronomy
community.

## Sources

For any reference and docs go to http://www.starlight.ufsc.br/

## How-to

## Build from the repository
```
# Get the source code from GitHub.
git clone https://github.com/MarioChC/astro-containers.git
cd astro-containers

# Go in the starlight folder.
cd starlight

# Create directories.
mkdir shared_directory
mkdir shared_directory/config_files_starlight
mkdir shared_directory/config_files_starlight/mask
mkdir shared_directory/config_files_starlight/spectrum

# Build the container.
docker build -t astro-containers/starlight .

# Check the image id. 
docker image ls

# Run Starlight in detached mode (-d) from this same folder and leave it running.
docker run -d -v <PATH>/astro-containers/starlight/shared_directory/:/home/starlight/shared_directory/ --name starlight_container <image_id> sleep infinity

# Check the container id.
docker ps

# Copy the ".config" file to the STARLIGHTv04 directory.
docker cp ../config_files_starlight/StCv04.C11.arp220.config <container_id_or_name>:/home/starlight/STARLIGHTv04/

# Copy all the configuration files to the "shared_directory"
cp <PATH>/config_files_starlight/grid_example.in shared_directory/config_files_starlight/
cp <PATH>/config_files_starlight/mask/* shared_directory/config_files_starlight/mask/
cp <PATH>/config_files_starlight/spectrum/* shared_directory/config_files_starlight/spectrum/

# Run the analysis with starlight and the config files
docker exec <container_id> /home/starlight/STARLIGHTv04/bash_script.sh

# The output files will be stored in your computer in the "astro-containers/starlight/shared_directory/output directory"

# Adjust your data files and execute it as you need.
```
