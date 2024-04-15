FROM ubuntu:focal
MAINTAINER Cristina Catalan <cataln.cris@gmail.com>

WORKDIR /home/starlight/STARLIGHTv04

COPY STARLIGHTv04/* /home/starlight/STARLIGHTv04/
COPY shared_directory/* /home/starlight/shared_directory/

RUN ls
RUN pwd

RUN echo "deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe" >> /etc/apt/sources.list
RUN echo "deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe" >> /etc/apt/sources.list
RUN echo "deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe" >> /etc/apt/sources.list
RUN echo "deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe" >> /etc/apt/sources.list

RUN apt update -y
RUN apt install -y gfortran

#COPY grid_example.in /home/starlight/grid_example.in
#COPY out/* /home/starlight/out/
#COPY mask/* /home/starlight/mask/
#COPY spectrum/* /home/starlight/spectrum/
#COPY StCv04.C11.arp220.config /home/starlight/STARLIGHTv04/
#COPY README.md /home/starlight/README.md
#COPY Dockerfile /home/starlight/Dockerfile
#COPY bash_script.sh /home/starlight/STARLIGHTv04/
