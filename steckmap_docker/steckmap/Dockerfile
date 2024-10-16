# Use a base image of Ubuntu Focal
FROM ubuntu:focal
MAINTAINER Mario Chamorro <mchamorro@ucm.com>

# Set the timezone non-interactively
ENV DEBIAN_FRONTEND=noninteractive
RUN ln -fs /usr/share/zoneinfo/Europe/Madrid /etc/localtime && \
    apt-get update && apt-get install -y tzdata

# Update the system and install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    emacs \
    yorick \
    unzip \
    xvfb \
    git \
    && rm -rf /var/lib/apt/lists/*

# Clone the full STECKMAP repository from GitHub
RUN git clone https://github.com/pocvirk/STECKMAP.git /root/Yorick/STECKMAP

# Set the working directory
WORKDIR /root/Yorick

# Define the default command to use xvfb with auto-servernum
CMD ["sleep", "infinity"]
