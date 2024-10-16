# Steckmap Docker Container

This project provides a preconfigured Docker container for running the Steckmap astronomy code. Steckmap is a powerful tool used in the astronomical community for stellar population synthesis and spectral fitting. By encapsulating Steckmap in a Docker container, we facilitate its execution without worrying about environment dependencies.

## Sources

For any reference and documentation, please visit the official Steckmap repository or documentation site https://github.com/pocvirk/STECKMAP.

## How-to Build from the Repository
```
# Clone the repository from GitHub:
git clone https://github.com/MarioChC/astro-containers.git
cd astro-containers/steckmap_docker

# Enter the "steckmap" folder:
cd steckmap

# Build the Docker container:
docker build -t steckmap_image .

# Check the image ID:
docker image ls

# Create the shared directories:
mkdir shared_directory

# Create the output directory in the shared volume:
mkdir shared_directory/output

# Copy the directory with your input files:
cp -r ../input shared_directory/

# Run Steckmap container in detached mode (-d) from the same folder, mounting a shared volume between the local machine and the container, and leave it running:
docker run -d --name steckmap_container -v <LOCAL_PATH>/astro-containers/steckmap_docker/steckmap/shared_directory/:/home/steckmap/shared_directory/ <IMAGE_ID> sleep infinity

# Check the container ID:
docker ps

# Run the analysis with Steckmap and the configuration files:
docker exec steckmap_container /home/steckmap/shared_directory/input/run_steckmap.i

# The output files will be stored in your computer in the "<LOCAL_PATH>/astro-containers/steckmap_docker/steckmap/shared_directory/output/" directory.

# Adjust your data files and execute it as you need.

