# Voronoi Rebinning Docker Container

This project provides a preconfigured Docker container for running Voronoi rebinning on astronomical data cubes. Voronoi rebinning is a technique used in astronomy to improve the signal-to-noise ratio in their data. By encapsulating this process in a Docker container, we facilitate its execution without worrying about environment dependencies.

## Sources

For any reference and docs go to https://www-astro.physics.ox.ac.uk/~cappellari/software/

## How-to build from the repository
```
# Clone the repository from GitHub:
git clone https://github.com/MarioChC/astro-containers.git
cd astro-containers/voronoi_docker

# Enter the "voronoi" folder:
cd voronoi

# Build the Docker container:
docker build -t voronoi_image .

# Check the image ID:
docker image ls

# Create the shared directories:
mkdir shared_directory

# Create the output directory in the shared volume:
mkdir shared_directory/output

# Copy the directory with yout input files:
cp -r ../input shared_directory/

# Run voronoi container in detached mode (-d) from the same folder, mounting a shared volume between the local machine and the container, and leave it running:

docker run -d --name voronoi_container -v <LOCAL_PATH>/astro-containers/voronoi_docker/voronoi/shared_directory/:/home/voronoi/shared_directory/ <IMAGE_ID> sleep infinity

# Check the container ID:
docker ps

# Run the analysis with Starlight and the configuration files:
docker exec voronoi_container /home/voronoi/shared_directory/input/bash_script.sh

# The output files will be stored in your computer in the "<LOCAL_PATH>/astro-containers/voronoi_docker/voronoi/shared_directory/output" directory.

# Adjust your data files and execute it as you need.
```
