FROM ubuntu:focal
MAINTAINER Mario Chamorro <mchamorro@ucm.com>

# Install Python and pip
RUN apt-get update && \
    apt-get install -y python3 python3-pip

WORKDIR /home/ppxf/

RUN echo "deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe" >> /etc/apt/sources.list
RUN echo "deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe" >> /etc/apt/sources.list
RUN echo "deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe" >> /etc/apt/sources.list
RUN echo "deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe" >> /etc/apt/sources.list


# Install ppxf and other dependencies using pip
RUN pip install --no-cache-dir ppxf vorbin plotbin

RUN apt update -y

COPY ./* /home/ppxf/run_ppxf/
COPY run_ppxf/ssp_models/* /usr/local/lib/python3.8/dist-packages/ppxf/sps_models/
COPY run_ppxf/modified_scripts/ppxf.py /usr/local/lib/python3.8/dist-packages/ppxf/
COPY run_ppxf/modified_scripts/sps_util.py /usr/local/lib/python3.8/dist-packages/ppxf/
COPY run_ppxf/modified_scripts/voronoi_2d_binning.py /usr/local/lib/python3.8/dist-packages/vorbin/