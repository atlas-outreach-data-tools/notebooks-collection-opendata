[![Docker image](https://github.com/atlas-outreach-data-tools/notebooks-collection-opendata/actions/workflows/docker_push.yml/badge.svg)](https://github.com/atlas-outreach-data-tools/notebooks-collection-opendata/actions/workflows/docker_push.yml) [![Run all Jupyter Notebooks](https://github.com/atlas-outreach-data-tools/notebooks-collection-opendata/actions/workflows/notebooks-health-check.yml/badge.svg)](https://github.com/atlas-outreach-data-tools/notebooks-collection-opendata/actions/workflows/notebooks-health-check.yml)
# ATLAS Open Data Jupyter Notebooks
A set of multiple notebooks using 8 TeV and 13 TeV ATLAS Open Data datasets.

## Run on Binder
To execute in MyBinder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/atlas-outreach-data-tools/notebooks-collection-opendata/HEAD)

For more complex needs, such as running multiple notebooks review the [Binder Documentation](https://mybinder.readthedocs.io/en/latest/#) for guidance on setting up a more personalized Binder environment.

Note: before starting running the code in the jupyter notebooks, click on the up right button "not trusted" in order to get "trusted" displayed. This should lead the JavaScript to be executed, that is useful to visualise interactive histograms. If that doesn't work, simply go to the top of the notebook, find the cell that contains the line of code `%jsroot` and comment out that.

## Run on Docker containers

Docker provides a robust platform for developing, sharing, and running applications within containers.

**Steps to run ATLAS open data using Docker:**

1. **Install Docker**: Download and install Docker from the [official website](https://docs.docker.com/get-docker/). 
2. **Open Docker**
3.  **Run the Docker Image**: Open your terminal and run the Docker image, available on our [GitHub registry](https://github.com/atlas-outreach-data-tools/notebooks-collection-opendata/pkgs/container/notebooks-collection-opendata).
```
docker run -it -p 8888:8888 -v my_volume:/home/jovyan/work ghcr.io/atlas-outreach-data-tools/notebooks-collection-opendata:latest
```
4. **Launch the Notebook**: After running the Docker image, use the link provided in the terminal to access the Jupyter interface. The terminal should look like this:
```
To access the notebook interface, open this file in a browser:
        file:///home/jovyan/.local/share/jupyter/runtime/nbserver-15-open.html
    Or copy and paste one of these URLs:
        http://4c61742ed77c:8888/?token=34b7f124f6783e047e796fea8061c3fca708a062a902c2f9
     or http://127.0.0.1:8888/?token=34b7f124f6783e047e796fea8061c3fca708a062a902c2f9
```

For more, please go to: [opendata.atlas.cern/docs/category/choose-your-enviroment](https://opendata.atlas.cern/docs/category/choose-your-enviroment)

@2025
