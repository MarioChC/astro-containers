FROM ubuntu:focal
MAINTAINER Mario Chamorro <mchamorro@ucm.com>

# Install Python and pip
RUN apt-get update && \
    apt-get install -y python3 python3-pip

WORKDIR /home/voronoi/

RUN echo "deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe" >> /etc/apt/sources.list
RUN echo "deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe" >> /etc/apt/sources.list
RUN echo "deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe" >> /etc/apt/sources.list
RUN echo "deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe" >> /etc/apt/sources.list

# Install ppxf and other dependencies using pip
RUN pip install --no-cache-dir ppxf vorbin

RUN apt update -y

COPY ./* /home/voronoi/run_voronoi/
RUN patch /usr/local/lib/python3.8/dist-packages/ppxf/ppxf.py run_voronoi/modified_scripts/patch_ppxf.diff
RUN patch /usr/local/lib/python3.8/dist-packages/ppxf/sps_util.py run_voronoi/modified_scripts/patch_sps_util.diff
RUN patch /usr/local/lib/python3.8/dist-packages/vorbin/voronoi_2d_binning.py run_voronoi/modified_scripts/patch_voronoi_2d_binning.diff
