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

USER root
# Create a new non-root user
ARG NB_USER=atlasuser
ARG NB_UID=1001
ARG NB_GID=100
RUN adduser --disabled-password --gecos "" --uid $NB_UID --gid $NB_GID $NB_USER && \
    mkdir -p /home/$NB_USER && \
    chown -R $NB_USER:$NB_GID /home/$NB_USER /opt/app

# Switch to the created user
USER $NB_USER

# Set the home directory as the working directory
WORKDIR /home/$NB_USER
