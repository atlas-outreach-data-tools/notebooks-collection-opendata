FROM jupyter/scipy-notebook:python-3.11
LABEL author="ATLAS Open Data Team @ CERN 2024"
LABEL maintainer="Giovanni Guerrieri - giovanni.guerrieri@cern.ch"
WORKDIR /opt/app

# Copy the environment.yml file into the container
COPY binder/environment.yml .

RUN conda env update -n base -f environment.yml && \
    conda clean -afy && \
    rm -rf /opt/conda/pkgs/*
SHELL ["bash", "-c"]

# Expose port 8888 for Jupyter Lab
EXPOSE 8888

WORKDIR /home
