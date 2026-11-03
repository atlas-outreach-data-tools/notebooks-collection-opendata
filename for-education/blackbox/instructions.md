# ATLAS Open Data Black Box
A black box dataset is a dataset containing many simulated events, which are initially unknown. The labels of the events are stored separately, making the dataset appear as it would in a real particle physics analysis. This tutorial shows you how you can create your own black box containing the processes and variables you like.


## Set up the environment
First, you have to clone the ATLAS Open Data Git repository. Open your terminal and run:
```
git clone https://github.com/atlas-outreach-data-tools/notebooks-collection-opendata atlas-open-data
cd atlas-open-data/for-education/blackbox
```
You are now in the directory containing all the files required to create your own black box.

To execute the code, several Python packages are required. To ensure that the correct Python version and all required packages are available, you can use `pyenv` to create an environment that meets all requirements. First, install `pyenv` by following [these](https://github.com/pyenv/pyenv#installation/#) instructions. 

The required packages are listed in `environment.txt`. To install them, you can use the `create_environment.sh` script. Make the script executable and run it:
```
chmod +x create_environment.sh
./create_environment.sh environment.txt
```
An environment named `blackbox` is created. To activate it, run: 
```
pyenv activate blackbox
```
To deactivate it, run `pyenv deactivate blackbox`.


## Specify your parameters in the config
In the config file, the content of the black box is specified. You can use the template `tmp_config.json`, while a concrete example is provided in `Zprime_config.json`.
First, you can define the signal and background dataset IDs (DSIDs). An overview of the available $13~\text{TeV}$ 2025 data can be found [here](https://opendata.atlas.cern/docs/data/for_education/13TeV25_metadata/#). 
Next, the desired variables are listed. [Here](https://opendata.atlas.cern/docs/data/for_education/13TeV25_details/#), you can find all available variable names, as well as the possible skims, which act as a preselection for the events.
To define the number of events that enter the black box, the integrated luminosity [ $\text{fb}^{-1}$ ] is defined. It is a measure of the total amount of collisions collected by an experiment over a given period of time. Most datasets are available at an integrated luminosity of $36~\text{fb}^{-1}$. To make the signal process more prominent in the black box, a signal scale factor can be defined. Alternatively, the absolute numbers of signal and background events can be defined in the `samples` block via the keys `signal_nevents` and `background_nevents`.
Finally, the name of the output directory is defined.


## Run the code to create the black box
Now you can run the code with your config:
```
python create_blackbox.py your_config.json
```
A checkpoint directory is created. In its events subdirectory, the events for each DSID are stored in a Parquet file once processing of the DSID is complete. Once all DSIDs have been processed, all events are combined and shuffled. The event variables are saved in blackbox_data.parquet, while the labels are stored in blackbox_labels.parquet.


## Explore the Output
In the interactive Jupyter notebook `explore_blackbox.ipynb` you can explore the structure of the black box. An example funtion to load the black box is also included.

Furthermore, you can train a neural network on the example black box $Z^{\prime} \rightarrow \mu\mu$ to classify signal and background events.


## Example: $Z^{\prime} \rightarrow \mu\mu$
An example black box is created for the signal process $Z^{\prime} \rightarrow \mu\mu$. The $Z^{\prime}$ is a hypothetical boson predicted by theories beyond the Standard Model. For our example, it has the same properties as the Standard Model $Z$ boson but has a much higher mass, in this example $3000~\text{GeV}$. The example config is `Zprime_config.json`. The included background processes are Drell-Yan, tt̄, single-top and diboson processes. Basic muon and jet properties, as well as the missing transverse energy, are included as variables. The luminosity is defined to be $11~\text{fb}^{-1}$, as some datasets have not enough events available to allow a higher value. To significantly enhance the signal, the scale factor is chosen to be 1000. No skim is applied and the output directory is `Zprime_blackbox_checkpoint`.
You can run the code to create the $Z'$ black box:
```
python create_blackbox.py Zprime_config.json
```
You can now find the black box in the output directory.
