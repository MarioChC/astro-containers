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
