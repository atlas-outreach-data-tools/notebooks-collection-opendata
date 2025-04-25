FROM quay.io/jupyter/scipy-notebook:python-3.11.8
LABEL author="ATLAS Open Data Team @ CERN 2024"
LABEL maintainer="Giovanni Guerrieri - giovanni.guerrieri@cern.ch"
WORKDIR /opt/app

# Copy the environment.yml file into the container
COPY binder/environment.yml .

# Install mamba
RUN conda install -n base -c conda-forge mamba && \
    conda clean -afy

# Use mamba to update environment
RUN mamba env update -n base -f environment.yml && \
    conda clean -afy && \
    rm -rf /opt/conda/pkgs/*

SHELL ["bash", "-c"]

# Expose port 8888 for Jupyter Lab
EXPOSE 8888

# Set the home directory as the working directory and default user
WORKDIR $HOME