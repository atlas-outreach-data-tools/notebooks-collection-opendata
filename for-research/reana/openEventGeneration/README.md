# ATLAS Open Data REANA Pipeline

This repository provides a reproducible data analysis pipeline for ATLAS Open Data, using Snakemake as the workflow engine and REANA for scalable execution. The pipeline covers physics event generation and detector simulation using Delphes, containerized via Jupyter and Atlas Outreach Data Tools images.

## Overview

The pipeline automates:
1. **Event generation and extraction** by executing parameterized Jupyter notebooks using Papermill.
2. **Fast detector simulation** with Delphes, producing analysis-ready outputs.

## Files and Structure

- `snakemake/Snakefile`: The Snakemake workflow definition.
- `snakemake/inputs.yaml`: Configuration file specifying analysis parameters (input/output files, dataset IDs, etc).
- `OpenEvgenTutorial.ipynb`: Notebook for event generation/extraction.
- `reana-snakemake.yaml`: REANA workflow specification linking inputs, Snakemake, and outputs.
- `results/`: Directory for output files (plots, processed datasets).

## Workflow Steps

### Rule: extract_EVGEN
- **Purpose:** Runs a Jupyter notebook with specified physics release and dataset ID, to extract EVGEN events and create output files.
- **Container:** `ghcr.io/atlas-outreach-data-tools/notebooks-collection-opendata:latest`
- **Command:** Uses Papermill for parameterized notebook execution.  
  Produces `$evgen_file` and stores results in `results/`.

### Rule: delphes_SIM
- **Purpose:** Runs Delphes simulation on extracted EVGEN events, using a specific detector card.
- **Container:** `ghcr.io/atlas-outreach-data-tools/notebooks-collection-opendata-delphes:latest`
- **Command:** Runs `DelphesHepMC2` with input event file and outputs simulation results.

## Running Locally (Without REANA)

1. Create and activate a local Python virtual environment:
   ```bash
   mkdir snakemake-local-run
   cd snakemake-local-run
   python -m venv atlas-env
   source atlas-env/bin/activate
   pip install snakemake papermill ipykernel
   cp -a ../OpenEvgenTutorial.ipynb .
   ```
2. Run workflow:
   ```bash
   snakemake -s ../snakemake/Snakefile \
     --configfile ../snakemake/inputs.yaml \
     --config notebook=OpenEvgenTutorial.ipynb -p --cores 1
   open results/plot.png
   ```

## Running on REANA

1. Upload all required files: Snakefile, configuration, notebooks.
2. Launch the workflow using REANA client:
   ```bash
   reana-client run -w <workflow-name>
   ```
   REANA executes the Snakemake workflow, handling container orchestration and resource allocation.

3. Retrieve results from the `results/` directory.

## Customization

- Modify `snakemake/inputs.yaml` to change datasets, output locations, or analysis parameters.
- Adapt notebooks in place for specific studies or plots.

## Output

- The pipeline produces analysis plots and simulation files in the `results/` directory, as specified in `inputs.yaml` and the workflow.