[![Docker image](https://github.com/atlas-outreach-data-tools/notebooks-collection-opendata/actions/workflows/docker_push.yml/badge.svg)](https://github.com/atlas-outreach-data-tools/notebooks-collection-opendata/actions/workflows/docker_push.yml) [![Run all Jupyter Notebooks](https://github.com/atlas-outreach-data-tools/notebooks-collection-opendata/actions/workflows/notebooks-health-check.yml/badge.svg)](https://github.com/atlas-outreach-data-tools/notebooks-collection-opendata/actions/workflows/notebooks-health-check.yml)
# ATLAS Open Data Jupyter Notebooks
A set of multiple notebooks using 8 TeV and 13 TeV ATLAS Open Data datasets.

## Running the notebooks
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/atlas-outreach-data-tools/notebooks-collection-opendata/HEAD)  [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/atlas-outreach-data-tools/notebooks-collection-opendata.git)  [![SWAN](https://swan.web.cern.ch/sites/swan.web.cern.ch/files/pictures/open_in_swan.svg)](https://cern.ch/swanserver/cgi-bin/go?projurl=https://github.com/atlas-outreach-data-tools/notebooks-collection-opendata.git)

### Run on Binder
For more complex needs, such as running multiple notebooks review the [Binder Documentation](https://mybinder.readthedocs.io/en/latest/#) for guidance on setting up a more personalized Binder environment.

Note: before starting running the code in the jupyter notebooks, click on the up right button "not trusted" in order to get "trusted" displayed. This should lead the JavaScript to be executed, that is useful to visualise interactive histograms. If that doesn't work, simply go to the top of the notebook, find the cell that contains the line of code `%jsroot` and comment out that.

### Run on Colab
Much of the functionality is the same as in Binder, but with some subtle differences that are explained in the [Colab Documentation](https://colab.research.google.com).

### Run on SWAN
To use SWAN you need a CERN account. Similar to Binder and Colab, but you have access to more memory and other resources, and faster access to the datasets. If you are new to SWAN, you can use the default environment setup to get running. For more information check out the [SWAN Documentation](https://swan.docs.cern.ch).

### Run in Docker containers

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

## Contribution guide
We appreciate the willingness to participate in the ATLAS Open Data project. Below you can find some steps to help you contribute successfully (inspired by [The Turing Way](https://book.the-turing-way.org) contribution guidelines):

### 1. Comment on an existing issue or open a new issue referencing your addition
Before creating a pull request (PR), discuss the change you want to make. That's the easiest way we can check that the work is not overlapping with current work and that it aligns with the goals of the Open Data project. 

In [this blog](https://www.igvita.com/2011/12/19/dont-push-your-pull-requests/) you can find some reasons why putting this work in upfront is so useful to everyone involved.

When opening a new issue be sure to fill all the information necessary. The issue template should help you doing so.

### 2. Fork the repository
Make a fork of the repository in the upper right corner of the repository page in GitHub. In this new copy, changes won't affect anyone else's work, and you can make as many changes as you want. 

Make sure to keep the fork up to date (or develop from a branch based on the main repository's main branch) to avoid conflicts when merging. If you already have conflicts, check [GitHub's guide to resolving a merge conflict](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/addressing-merge-conflicts/resolving-a-merge-conflict-on-github).

### 3. Make the changes you've discussed
Please keep the changes targeted to what was discussed. While making your changes, commit often and write explanatory commit messages.

Please do not re-write history! That is, please do not use the rebase command to edit previous commit messages, combine multiple commits into one, or delete or revert commits that are no longer necessary.

### 4. Submit a pull request
Once you are done with your changes, open a PR. We encourage you to open a PR as early in your contributing process as possible so that everyone can see what is being worked on.

When submitting your PR you will see a template. It will ask:
- To describe the problem you are trying to fix an reference related issues. 
- List the changes you are proposing
- Tell us what we should focus our review on.

Please give us all the details, this makes it easier for us to review the contribution!

@2024
