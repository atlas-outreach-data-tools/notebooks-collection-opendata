# Keep the Python version in line with binder/environment.yml.
FROM quay.io/jupyter/scipy-notebook:python-3.13.18
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

# Fail the build if the compiled packages don't load with the installed NumPy,
# rather than letting the notebooks fail in CI.
RUN python -c "import numpy, scipy.spatial, matplotlib.pyplot, pandas, sklearn, numba, lmfit, awkward, uproot, vector, ROOT; print('numpy', numpy.__version__)"

# Expose port 8888 for Jupyter Lab
EXPOSE 8888

USER root

RUN chgrp -R 0 /home /opt/app && \
    chmod -R g=u /home /opt/app && \
    chmod -R g+w /home /opt/app

SHELL ["bash", "-c"]

# Set the home directory as the working directory and default user
WORKDIR $HOME

# # Set the default user and group IDs, from https://quay.io/repository/jupyter/scipy-notebook/manifest/sha256:44e4d086d70039ad2ea67d26ed22a2b1fd5ebd24b2ab47364baa77f0ec9549f4
ARG NB_GID=100
ARG NB_UID=1000
ARG NB_USER=jovyan