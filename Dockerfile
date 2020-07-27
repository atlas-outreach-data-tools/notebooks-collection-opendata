# Dockerfile

# Specify the base image that we're building the image on top of
FROM python:3.7-slim

# Copy this ttZ-2lOS-preFit directory to a new directory within the docker image
COPY . /notebooks-collection-opendata/


# Change directory to save input MC files
WORKDIR /notebooks-collection-opendata/13-TeV-examples/uproot_python/Input/4lep/MC/ 

# Update pip, 
# Run some commands to install packages, 
# Have Jupyter notebooks launch without additional command line options
# Get input MC files in docker container
RUN pip install --no-cache-dir -q --upgrade pip && \
    pip install --no-cache-dir -q jupyter \
				  uproot \
				  pandas \
				  numpy \
				  matplotlib \
				  lmfit \
				  scikit-learn \
    				  wget && \
    # Jupyter config commands don't seem to work on GitHub workflows :(
    #jupyter notebook --generate-config && \
    #sed -i -e "/allow_root/ a c.NotebookApp.allow_root = True" ~/.jupyter/jupyter_notebook_config.py && \
    #sed -i -e "/custom_display_url/ a c.NotebookApp.custom_display_url = \'http://localhost:8888\'" ~/.jupyter/jupyter_notebook_config.py && \
    #sed -i -e "/c.NotebookApp.ip/ a c.NotebookApp.ip = '0.0.0.0'" ~/.jupyter/jupyter_notebook_config.py && \
    #sed -i -e "/open_browser/ a c.NotebookApp.open_browser = False" ~/.jupyter/jupyter_notebook_config.py && \
    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/MC/mc_361106.Zee.4lep.root && \
    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/MC/mc_361107.Zmumu.4lep.root && \
    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/MC/mc_341947.ZH125_ZZ4lep.4lep.root && \
    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/MC/mc_341964.WH125_ZZ4lep.4lep.root && \
    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/MC/mc_344235.VBFH125_ZZ4lep.4lep.root && \
    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/MC/mc_345060.ggH125_ZZ4lep.4lep.root && \
    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/MC/mc_363490.llll.4lep.root && \
    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/MC/mc_410000.ttbar_lep.4lep.root && \
    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/MC/mc_307434.RS_G_ZZ_llll_c10_m0500.4lep.root && \
    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/MC/mc_410155.ttW.4lep.root && \
    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/MC/mc_410218.ttee.4lep.root && \
    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/MC/mc_410219.ttmumu.4lep.root

# Change directory to save input 4lep/Data files
WORKDIR /notebooks-collection-opendata/13-TeV-examples/uproot_python/Input/4lep/Data/ 

# Get input 4lep/Data files in docker container
RUN python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/Data/data_A.4lep.root && \
    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/Data/data_B.4lep.root && \
    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/Data/data_C.4lep.root && \
    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/Data/data_D.4lep.root

# Change directory to save input GamGam/Data files
WORKDIR /notebooks-collection-opendata/13-TeV-examples/uproot_python/Input/GamGam/Data/

# Get input GamGam/Data files in docker container
RUN python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/GamGam/Data/data_A.GamGam.root

# Get extra input GamGam/Data files in docker container
#RUN python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/GamGam/Data/data_B.GamGam.root && \
#    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/GamGam/Data/data_C.GamGam.root && \
#    python -m wget https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/GamGam/Data/data_D.GamGam.root

# This sets the default working directory when a container is launched from the image
WORKDIR /notebooks-collection-opendata/13-TeV-examples/uproot_python/

# Add Tini. Tini operates as a process subreaper for jupyter. This prevents kernel crashes.
ENV TINI_VERSION v0.6.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
ENTRYPOINT ["/usr/bin/tini", "--"]
