FROM rootproject/root:6.32.02-ubuntu22.04
LABEL maintainer="giovanni.guerrieri@cern.ch"
WORKDIR /
ENV DEBIAN_FRONTEND=noninteractive
USER root

#install the prerequisites (option always yes activated)
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y python3 python3-dev git curl python3-pip \                                           
    && apt-get --yes install \
    nano\
    vim


RUN pip3 install --no-cache-dir --upgrade  \
    uproot3 \
    uproot \
    tables \
    matplotlib \
    pandas \
    pydot \
    awkward \
    vector \
    scikit-learn \
    lmfit \
    jupyter \ 
    ipykernel \ 
    papermill \
    coffea

WORKDIR /home/jovyan/work/

