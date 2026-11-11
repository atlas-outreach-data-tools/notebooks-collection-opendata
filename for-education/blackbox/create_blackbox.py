import sys
import os
import argparse

import atlasopenmagic as atom

import uproot # for reading .root files
import awkward as ak # for handling complex and nested data structures efficiently
import numpy as np # # for numerical calculations such as histogramming
import matplotlib.pyplot as plt # for plotting
import json


def get_inclusive_yield(lumi, metadata):
    return (
        lumi * 1000
        * metadata["cross_section_pb"]
        * metadata["genFiltEff"]
        * metadata["kFactor"]
    )


def create_objects(data, variables):

    prefixes = {
        "lep": "lep_",
        "jet": "jet_",
        "tau": "tau_",
        "photon": "photon_",
        "largeRJet": "largeRJet_",
        "ScaleFactor": "ScaleFactor_",
        "trig": "trig_",
        "met": "met_"}

    objects = {}

    for variable in variables:
        # special case: met
        if variable == "met":
            if "met" not in objects:
                objects["met"] = {}
            objects["met"]["met"] = data[variable]
            continue

        found = False

        for category, prefix in prefixes.items():
            if variable.startswith(prefix):
                name = variable[len(prefix):]
                if category not in objects:
                    objects[category] = {}
                objects[category][name] = data[variable]
                found = True
                break

        if not found:
            if "info" not in objects:
                objects["info"] = {}
            objects["info"][variable] = data[variable]

    # zip everything
    for category in objects:
        objects[category] = ak.zip(objects[category], depth_limit=1)

    return objects


def create_blackbox(signal_dsid_list, background_dsid_list, variables, signal_nevents, background_nevents, 
                    lumi=36., signal_scale_factor=1., skim="noskim", checkpoint_dir="blackbox_checkpoint",):
    
    print("Variables:", variables)

    # create checkpoint directories
    events_dir = os.path.join(checkpoint_dir, "events")
    completed_file = os.path.join(checkpoint_dir, "completed_dsids.json")

    os.makedirs(events_dir, exist_ok=True)


    # load checkpoint information
    if os.path.exists(completed_file):

        print(f"Load checkpoint: {checkpoint_dir}")

        with open(completed_file, "r") as f:
            completed_dsids = set(json.load(f))

        print(f"{len(completed_dsids)} DSIDs already processed.")
        print(f"Completed DSIDs: {completed_dsids}")

    else:
        print("No checkpoint found. Start from scratch.")
        completed_dsids = set()


    # calculate inclusive yields
    yield_dict = {}

    signal_yield_sum = 0
    background_yield_sum = 0
    nevents_per_sample_dict = {}

    for dsid in signal_dsid_list:

        metadata = atom.get_metadata(dsid)
        inc_yield = get_inclusive_yield(lumi, metadata)

        yield_dict[dsid] = inc_yield * signal_scale_factor
        signal_yield_sum += inc_yield * signal_scale_factor


    for dsid in background_dsid_list:

        metadata = atom.get_metadata(dsid)
        inc_yield = get_inclusive_yield(lumi, metadata)

        yield_dict[dsid] = inc_yield
        background_yield_sum += inc_yield


    # calculate number of events per DSID
    if signal_nevents is not None and background_nevents is not None:
        
        for dsid in signal_dsid_list:
            nevents_per_sample = (signal_nevents / signal_yield_sum * yield_dict[dsid])
            nevents_per_sample_dict[dsid] = nevents_per_sample


        for dsid in background_dsid_list:
            nevents_per_sample = (background_nevents / background_yield_sum * yield_dict[dsid])
            nevents_per_sample_dict[dsid] = nevents_per_sample

    else:
        background_nevents = background_yield_sum
        nevents_per_sample_dict = yield_dict

    # calculate available number of events for each sample to check if is is sufficient
    for dsid in signal_dsid_list + background_dsid_list:
        numevents = 0
        metadata = atom.get_metadata(dsid)
        file_list = atom.get_urls(dsid, skim=skim, protocol='root', cache=False)

        for afile in file_list:
            tree = uproot.open(afile + ":analysis")
            numevents += tree.num_entries

        if numevents < nevents_per_sample_dict[dsid]:
            raise ValueError(f"WARNING: Only found {numevents} for DSID {dsid}, but {nevents_per_sample_dict[dsid]} are required. "
                             "Select a lower luminosity or number of events or choose a looser skim.")


    # save completed DSIDs
    def save_completed_dsids():
        with open(completed_file, "w") as f:
            json.dump(sorted(completed_dsids), f, indent=2)

    # calculate chunk size for efficient processing
    step_size = int(background_nevents) // len(background_dsid_list)
    digits = len(str(abs(step_size)))
    step_size = (step_size // 10**(digits - 2)) * 10**(digits - 2)

    print(f"stepsize: {step_size}")

    step_size = min(step_size, 1000000)
    step_size = max(step_size, 100)


    # process samples
    dsid_labels = ([(dsid, 1) for dsid in signal_dsid_list] + [(dsid, 0) for dsid in background_dsid_list])

    for dsid, label in dsid_labels:

        dsid_str = str(dsid)

        print("\n" + "=" * 60)
        print(f"Current DSID: {dsid}")
        print(f"Label: {label}")
        print(f"Already processed: {dsid_str in completed_dsids}")

        # skip already completed DSID
        if dsid_str in completed_dsids:
            print(f"DSID {dsid} already processed -> SKIP")
            continue

        # number of events to collect
        target_events = int(nevents_per_sample_dict[dsid])

        remaining = target_events
        collected_events = 0

        print(f"Target events: {target_events}")

        # temporary storage for this DSID
        chunks = {}
        label_chunks = []

        # get ROOT files
        file_list = atom.get_urls(dsid, skim, protocol="root", cache=False)
        print("Number of files:", len(file_list))

        # loop over ROOT files
        for file_number, afile in enumerate(file_list, start=1):

            print(f"\nProcessing file " f"{file_number}/{len(file_list)}")

            if collected_events >= target_events:
                break

            # read ROOT file chunk-by-chunk
            for data in uproot.iterate(afile + ":analysis", variables, library="ak", step_size=step_size,):

                if len(data) == 0:
                    continue

                # only take events that are still needed
                if len(data) > remaining:
                    data = data[:remaining]

                n_events = len(data)

                # create objects from variable names
                new_objects = create_objects(data, variables)
                label_chunks.append(ak.Array([label] * n_events))

                # save in chunks
                for name, obj in new_objects.items():
                    if name not in chunks:
                        chunks[name] = []
                    chunks[name].append(obj)

                # update counters
                collected_events += n_events
                remaining = target_events - collected_events

                print(f"Collected events: " f"{collected_events}/{target_events}")

                if collected_events >= target_events:
                    break

        # check whether enough events were found
        if collected_events < target_events:
            raise ValueError(f"WARNING: Only found " f"{collected_events}/{target_events} events " f"for DSID {dsid}.")

        # combine chunks for this DSID
        dsid_arrays = {}

        for name, chunk_list in chunks.items():
            if len(chunk_list) > 0:
                dsid_arrays[name] = ak.concatenate(chunk_list, axis=0)

        # combine labels
        if len(label_chunks) == 0:
            print(f"WARNING: No events collected for DSID {dsid} -> SKIP")
            continue
        dsid_arrays["label"] = ak.concatenate(label_chunks, axis=0)

        # create one Awkward record for this DSID
        dsid_data = ak.zip(dsid_arrays, depth_limit=1)

        # save DSID checkpoint
        checkpoint_file = os.path.join(events_dir, f"dsid_{dsid_str}.parquet")

        print(f"Saving {len(dsid_data)} events to " f"{checkpoint_file}")

        ak.to_parquet(dsid_data, checkpoint_file)

        # mark DSID as completed
        completed_dsids.add(dsid_str)

        save_completed_dsids()


    # load all checkpoint files
    print("\nLoading all checkpoint files...")

    parquet_files = [os.path.join(events_dir, f) for f in os.listdir(events_dir) if f.endswith(".parquet")]

    parquet_files.sort()

    if len(parquet_files) == 0:
        raise RuntimeError("No events found in checkpoint.")

    # load and concatenate
    data = ak.from_parquet(parquet_files)

    print(f"Loaded {len(data)} events " f"from {len(parquet_files)} parquet files.")


    # shuffle
    rng = np.random.default_rng(25)
    indices = rng.permutation(len(data))
    data = data[indices]

    # labels
    labels = data["label"]

    # remove label from the variable dataset
    variable_fields = [field for field in data.fields if field != "label"]
    variables_data = data[variable_fields]

    print(f"Number of events: {len(variables_data)}")
    print(f"Number of labels: {len(labels)}")

    # save final files

    variables_file = os.path.join(checkpoint_dir, "blackbox_data.parquet")
    labels_file = os.path.join(checkpoint_dir, "blackbox_labels.parquet")
    print(f"Saving variables to {variables_file}")
    ak.to_parquet(variables_data, variables_file)
    print(f"Saving labels to {labels_file}")
    ak.to_parquet(labels, labels_file)

    print("Blackbox successfully created!")

    return



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create blackbox dataset")
    parser.add_argument("config", help="Path to the configuration JSON file")
    args = parser.parse_args()

    sys.path += [ f'{os.environ["HOME"]}/.local/lib/python{sys.version_info.major}.{sys.version_info.minor}/site-packages' ]
    atom.set_release('2025e-13tev-beta')

    CONFIG_FILE = args.config

    with open(CONFIG_FILE, "r") as f:
        config = json.load(f)

    signal_dsid_list = config["samples"]["signal"]["dsids"]
    background_dsid_list = config["samples"]["background"]["dsids"]

    variables = config["variables"]

    lumi = config["options"].get("lumi")
    signal_scale_factor = config["options"].get("signal_scale_factor")
    signal_nevents = config["options"].get("signal_nevents")
    background_nevents = config["options"].get("background_nevents")

    skim = config["options"]["skim"]
    checkpoint_dir = config["options"]["checkpoint_dir"]

    # check if DSID inputs are valid
    common_dsids = set(signal_dsid_list) & set(background_dsid_list)
    if common_dsids:
        raise ValueError(f"WARNING: The following DSIDs are present in both signal and background: {sorted(common_dsids)}")

    # check if lumi, scalefactor and number of events are valid inputs
    is_lumi = lumi is not None or signal_scale_factor is not None
    is_nevents = signal_nevents is not None or background_nevents is not None

    if is_lumi and is_nevents:
        raise ValueError("WARNING: Please specify either lumi and signal_scale_factor "
        "or signal_nevents and background_nevents, not both.")

    if is_lumi:
        if lumi is None:
            raise ValueError("WARNING: Please specify a luminosity.")
        if signal_scale_factor is None:
            signal_scale_factor = 1.
        create_blackbox(signal_dsid_list, background_dsid_list, variables, signal_nevents, background_nevents, 
                        lumi, signal_scale_factor, skim=skim, checkpoint_dir=checkpoint_dir)

    if is_nevents:
        if signal_nevents is None or background_nevents is None:
            raise ValueError("WARNING: Please specify both signal_nevents and background_nevents.")
        create_blackbox(signal_dsid_list, background_dsid_list, variables, signal_nevents, background_nevents, skim=skim, checkpoint_dir=checkpoint_dir)

    sys.exit(0)