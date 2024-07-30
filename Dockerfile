FROM jupyter/scipy-notebook:x86_64-ubuntu-22.04
LABEL maintainer="giovanni.guerrieri@cern.ch"
WORKDIR /
ENV DEBIAN_FRONTEND=noninteractive
USER root

#install the prerequisites (option always yes activated)
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y python3 python3-dev git curl python3-pip \                                           
    && apt-get --yes install \
    dpkg-dev  \
    cmake \
    g++ \
    gcc \
    binutils\
    libx11-dev\
    libxpm-dev\
    libxft-dev\
    libxext-dev\
    libssl-dev\
    nano\
    vim

RUN conda update \
jupyterlab \
notebook

#    ========================== 
#    Installing ROOT 
#    ========================== 
RUN conda install root -c conda-forge

RUN pip3 install --no-cache-dir --upgrade  \
    uproot3 \
    uproot \
    tables \
    matplotlib \
    pydot \
    awkward==1.2.0 \
    vector \
    scikit-learn \
    lmfit \
    jupyter \ 
    ipykernel \ 
    papermill \
    coffea


WORKDIR /home/jovyan/work/

