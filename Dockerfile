FROM continuumio/miniconda3
LABEL author="ATLAS Open Data Team @ CERN 2024"
LABEL maintainer="Giovanni Guerrieri - giovanni.guerrieri@cern.ch"
WORKDIR /opt/app

# Copy the environment.yml file into the container
COPY binder/environment.yml .

RUN conda env update -n base -f environment.yml && \
    conda clean -afy && \
    rm -rf /opt/conda/pkgs/*
SHELL ["conda", "run", "-n", "analysis", "/bin/bash", "-c"]

# # Expose port 8888 for Jupyter Lab
# EXPOSE 8888

# CMD ["jupyter", "lab", "--ip=0.0.0.0", "--allow-root", "--no-browser"]