#!/usr/bin/env python3

# Set up the options properly
import argparse

# For exit codes
import sys

# For getting file names
import glob

# For the eventual output file
import json

# For type hints
import typing

# For reading the input file
import numpy as np
import uproot
import awkward as ak
import coffea.processor as processor
from coffea.nanoevents.schemas.base import BaseSchema, zip_forms
from coffea.nanoevents.methods import base, vector
from coffea.nanoevents import NanoEventsFactory

# For dealing with warnings from uproot that we don't need to worry about
import warnings

# Largely following the ATLAS PHYSLITE tutorial here
# https://gitlab.cern.ch/atlas-analysis-sw-tutorial/plotting/-/blob/master/LXPLUS-plotting_the_z_peak.ipynb?ref_type=heads
# Opening the ntuple and reading ALL of the branches will take unnecessary time and memory
# Let's define a schema to only use the branches we care about

class PHYSLITE_NtupleSchema(BaseSchema):
    def __init__(self, base_form):
        super().__init__(base_form)
        self._form['contents'] = self._build_collections(self._form['contents'])

    # We don't need ALL of the branches from the ntuple
    # We only care about the electron variables
    def _build_collections(self, branch_forms):

        output = {}

        # Event information
        # For data, the relevant event identifiers are run number and event number
        # For MC, the relevant event identifiers are mc channel number (DSID) and MC event number
        output['EventID'] = zip_forms( {'event_number':branch_forms['EventInfoAuxDyn.eventNumber'],
                                        'run_number':branch_forms['EventInfoAuxDyn.runNumber'],
                                        'mc_event_number':branch_forms['EventInfoAuxDyn.mcEventNumber'],
                                        'channel_number':branch_forms['EventInfoAuxDyn.mcChannelNumber']},
                                       'EventID')
        # Jets, small R and large R - just need the 4-vectors
        for jetCollection in ['AnalysisJets','AnalysisLargeRJets']:
            output[jetCollection] = zip_forms( {'pt':branch_forms[jetCollection+'AuxDyn.pt'],
                                                'eta':branch_forms[jetCollection+'AuxDyn.eta'],
                                                'phi':branch_forms[jetCollection+'AuxDyn.phi'],
                                                'mass':branch_forms[jetCollection+'AuxDyn.m']},
                                                'AnalysisJets', 'PtEtaPhiMLorentzVector' )
        # Vertices - just need the 3-vector
        output['PrimaryVertices'] = zip_forms( {'x':branch_forms['PrimaryVerticesAuxDyn.x'],
                                                'y':branch_forms['PrimaryVerticesAuxDyn.y'],
                                                'z':branch_forms['PrimaryVerticesAuxDyn.z']},
                                                  'PrimaryVertices')
        # MET - just need etx and ety, along with the source to be sure we have the right MET
        output['AnalysisMET'] = zip_forms( {'source':branch_forms['MET_Core_AnalysisMETAuxDyn.source'],
                                            'etx':branch_forms['MET_Core_AnalysisMETAuxDyn.mpx'],
                                            'ety':branch_forms['MET_Core_AnalysisMETAuxDyn.mpy']},
                                           'AnalysisMET')
        # Tracks - we use the 5 standard parameters for all four track collections
        for trkCollection in ['InDetTrackParticles','MuonSpectrometerTrackParticles','ExtrapolatedMuonTrackParticles','GSFTrackParticles']:
            track_branch_dict = {'d0':branch_forms[trkCollection+'AuxDyn.d0'],
                                                'z0':branch_forms[trkCollection+'AuxDyn.z0'],
                                                'phi':branch_forms[trkCollection+'AuxDyn.phi'],
                                                'theta':branch_forms[trkCollection+'AuxDyn.theta'],
                                                'qOverP':branch_forms[trkCollection+'AuxDyn.qOverP']}
            if trkCollection+'AuxDyn.numberDoF' in branch_forms:
                track_branch_dict['ndof'] = branch_forms[trkCollection+'AuxDyn.numberDoF']
            if trkCollection+'AuxDyn.chiSquared' in branch_forms:
                track_branch_dict['chi2'] = branch_forms[trkCollection+'AuxDyn.chiSquared']
            output[trkCollection] = zip_forms( track_branch_dict,
                                               trkCollection)
        # Clusters - just the energy, eta, and phi
        output['egammaClusters'] = zip_forms( {'energy':branch_forms['egammaClustersAuxDyn.calE'],
                                               'eta':branch_forms['egammaClustersAuxDyn.calEta'],
                                               'phi':branch_forms['egammaClustersAuxDyn.calPhi']},
                                              'egammaClusters')

        # Muons - eta, phi, type, quality, and some links we need to keep
        output['AnalysisMuons'] = zip_forms( {'Eta':branch_forms['AnalysisMuonsAuxDyn.eta'],
                                              'Phi':branch_forms['AnalysisMuonsAuxDyn.phi'],
                                              'Type':branch_forms['AnalysisMuonsAuxDyn.muonType'],
                                              'Quality':branch_forms['AnalysisMuonsAuxDyn.quality'],
                                              'IDTP_Link':branch_forms['AnalysisMuonsAuxDyn.inDetTrackParticleLink/AnalysisMuonsAuxDyn.inDetTrackParticleLink.m_persIndex'],
                                              'MSTP_Link':branch_forms['AnalysisMuonsAuxDyn.muonSpectrometerTrackParticleLink/AnalysisMuonsAuxDyn.muonSpectrometerTrackParticleLink.m_persIndex'],
                                              'EMTP_Link':branch_forms['AnalysisMuonsAuxDyn.extrapolatedMuonSpectrometerTrackParticleLink/AnalysisMuonsAuxDyn.extrapolatedMuonSpectrometerTrackParticleLink.m_persIndex']},
                                             'AnalysisMuons')
        # Photons - four vector plus cluster link, like the example in the tutorial
        output['AnalysisPhotons'] = zip_forms( {'pt':branch_forms['AnalysisPhotonsAuxDyn.pt'],
                                                'eta':branch_forms['AnalysisPhotonsAuxDyn.eta'],
                                                'phi':branch_forms['AnalysisPhotonsAuxDyn.phi'],
                                                'mass':branch_forms['AnalysisPhotonsAuxDyn.m'],
                                                'Cluster_Links':branch_forms['AnalysisPhotonsAuxDyn.caloClusterLinks']},
                                               'AnalysisPhotons', 'PtEtaPhiMLorentzVector')
        # Electrons - four vector plus cluster and track links
        output['AnalysisElectrons'] = zip_forms( {'pt':branch_forms['AnalysisElectronsAuxDyn.pt'],
                                                'eta':branch_forms['AnalysisElectronsAuxDyn.eta'],
                                                'phi':branch_forms['AnalysisElectronsAuxDyn.phi'],
                                                'mass':branch_forms['AnalysisElectronsAuxDyn.m'],
                                                'Cluster_Links':branch_forms['AnalysisElectronsAuxDyn.caloClusterLinks'],
                                                'Track_Links':branch_forms['AnalysisElectronsAuxDyn.trackParticleLinks']},
                                               'AnalysisElectrons', 'PtEtaPhiMLorentzVector')
        # All done putting together our branches! Now just return the output that we've built
        return output
        
    @property
    def behavior(self):
        behavior = {}
        behavior.update(base.behavior)
        behavior.update(vector.behavior)
        return behavior


# Similar setup, but for flat ntuples provided using the outreach and education ntuple making framework
class Flat_NtupleSchema(BaseSchema):
    def __init__(self, base_form):
        super().__init__(base_form)
        self._form['contents'] = self._build_collections(self._form['contents'])

    # We don't need ALL of the branches from the ntuple
    # We only care about the electron variables
    def _build_collections(self, branch_forms):

        output = {}

        # Event information - None stored in the ntuple by default

        # Jets, small R and large R - just need the 4-vectors
        for jetCollection in ['jet','largeRJet']:
            output[jetCollection] = zip_forms( {'pt':branch_forms[jetCollection+'_pt'],
                                                'eta':branch_forms[jetCollection+'_eta'],
                                                'phi':branch_forms[jetCollection+'_phi'],
                                                'energy':branch_forms[jetCollection+'_e']},
                                                'AnalysisJets', 'PtEtaPhiELorentzVector' )

        # Vertices - None stored in the ntuple by default

        # MET - just need etx and ety
        output['AnalysisMET'] = zip_forms( {'etx':branch_forms['met_mpx'],
                                            'ety':branch_forms['met_mpy']},
                                           'AnalysisMET')

        # Tracks - None stored in the ntuple by default

        # Clusters - None stored in the ntuple by default

        # Photons - four vector plus cluster link, like the example in the tutorial
        output['AnalysisPhotons'] = zip_forms( {'pt':branch_forms['photon_pt'],
                                                'eta':branch_forms['photon_eta'],
                                                'phi':branch_forms['photon_phi'],
                                                'energy':branch_forms['photon_e']},
                                               'AnalysisPhotons', 'PtEtaPhiELorentzVector')

        # Electrons and Muons - four vector plus cluster and track-based info
        output['Leptons'] = zip_forms( {'pt':branch_forms['lep_pt'],
                                        'eta':branch_forms['lep_eta'],
                                        'phi':branch_forms['lep_phi'],
                                        'energy':branch_forms['lep_e'],
                                        'd0':branch_forms['lep_d0'],
                                        'z0':branch_forms['lep_z0'],
                                        'lep_type':branch_forms['lep_type'],
                                        'charge':branch_forms['lep_charge'],
                                        'tight':branch_forms['lep_isTight']},
                                       'AnalysisElectrons', 'PtEtaPhiELorentzVector')
        # All done putting together our branches! Now just return the output that we've built
        return output

    @property
    def behavior(self):
        behavior = {}
        behavior.update(base.behavior)
        behavior.update(vector.behavior)
        return behavior



# Main function that does the work and returns a dictionar
def json_format(files: list[str],
                events: int=10,
                skip: int=0,
                eventListFile: typing.Optional[str]=None) -> dict:
    '''Create a dictionary of events conforming to the Phoenix event display standard
    based on input files from a list.
    files - python list of all files to go through
    events - number of events to output
    skip - number of input events to skip
    eventListFile - file containing a list of run number / event pairs for output
    Phoenix format defined at https://github.com/HSF/phoenix/blob/main/guides/developers/event_data_format.md
    '''

    # Create the output dictionary
    output_dict = {}

    # Get the event list
    eventList = []
    # If an input file was given, make a list of run, event tuples
    if eventListFile is not None:
        with open(eventListFile,'r') as eventListInput:
            for the_line in eventListInput:
                # Allow for comments and commented out lines
                cut_line = the_line.split('#')[0].strip()
                if len(cut_line.strip())==0:
                    continue

                # See if we used a few obvious dividers
                if len(cut_line.split())==2:
                    eventList += [ ( int(cut_line.split()[0].strip()),int(cut_line.split()[1].strip()) ) ]
                elif len(cut_line.split(','))==2:
                    eventList += [ ( int(cut_line.split(',')[0].strip()),int(cut_line.split(',')[1].strip()) ) ]
                else:
                    print(f'Could not parse {cut_line} in {eventListFile} as run, event pair')
        if ntuple:
            print('FYI: When running on ntuples, run numbers are not available, and event numbers correspond to TTree entry numbers')
            # Modify the event list accordingly
            eventList = [ (1,x[1]) for x in eventList ]
    print(f'Will attempt to output {(len(eventList) if len(eventList)>0 else events)} events')

    # Check if this is MC - will do this in not the most elegant way, but can avoid a user setting anything
    # For flat ntuples this doesn't matter, in fact (set below)
    isMC = -1

    # Helpers for muon types and qualities
    muon_types = ['Combined', 'Standalone', 'SegmentTagged', 'CaloTagged', 'SiAssociatedForward']
    muon_quality = ['Tight', 'Medium', 'Loose', 'VeryLoose'];
    # Helpers for jet radii
    jet_radii = {'jet':0.4,'largeRJet':1.0,'AnalysisJets':0.4,'AnalysisLargeRJets':1.0}

    # List of track collections I'll want to save
    trkCollectionList = ['InDetTrackParticles','MuonSpectrometerTrackParticles','ExtrapolatedMuonTrackParticles','GSFTrackParticles']

    # Now the actual work! Start looping over the input files
    for input_file_name in files:
        # Open the input file using uproot and grab the tree we want
        input_file = uproot.open(input_file_name)
        ntuple = True if 'analysis' in input_file else False
        if ntuple and eventListFile is not None:
            print('FYI: When running on ntuples, run numbers are not available, and event numbers correspond to TTree entry numbers')
            # Modify the event list accordingly
            eventList = [ (1,x[1]) for x in eventList ]
        if ntuple:
            isMC = 1
        tree_name = 'analysis' if ntuple else 'CollectionTree'
        collection_tree = input_file[tree_name]

        # Check if we're running on MC simulation or data. Not super elegant, but very functional.
        if isMC<0:
            isMC = 'TruthMuons' in collection_tree.keys()

        # Not all our track collections have all info. Let's check which ones have the extra info we might want.
        trkCollectionsWithNdofChi2 = []
        if not ntuple:
            for trkCollection in trkCollectionList:
                if trkCollection+'AuxDyn.numberDoF' in collection_tree.keys() and trkCollection+'AuxDyn.chiSquared' in collection_tree.keys():
                    trkCollectionsWithNdofChi2 += [trkCollection]

        # This will raise a bunch of warnings about branches that we don't need to worry about
        # In order to avoid worrying the users, we will turn those off
        # Only necessary if we are running on PHYSLITE
        if not ntuple:
            warnings.filterwarnings('ignore', category=UserWarning, module='coffea')
        # Let's use the schema to open our example single ntuple file and make an events object (which is an awkward array)
        events_data = NanoEventsFactory.from_root( input_file,
                                                   treepath=tree_name,
                                                   schemaclass=Flat_NtupleSchema if ntuple else PHYSLITE_NtupleSchema ).events()
        # Now turn warnings back on
        warnings.resetwarnings()

        for processed_events in range(len(events_data)):
            # Output in case we're far enough along
            if (processed_events+1)%100==0:
                print(f'Processed {processed_events} events so far')
            # First priority: skip events if requested
            # Skip events will be a little funny if someone sets an event list - that's fine, that's a weird thing to do
            if skip>0:
                skip -= 1
                continue

            if ntuple:
                event_number = processed_events
                run_number = 1
            else:
                event_number = events_data['EventID'][processed_events].mc_event_number if isMC>0 else events_data['EventID'][processed_events].event_number
                # In MC simulation, the run number indicates the conditions used for simulation
                # The channel number is the dataset ID, which is more intuitive
                run_number = events_data['EventID'][processed_events].channel_number if isMC>0 else events_data['EventID'][processed_events].run_number

            # If we are using a list, check if the event is in our list
            if len(eventList)>0 and (run_number,event_number) not in eventList:
                # Extra catch for MC, because this has tripped up so many people
                conditions_run_number = 1 if ntuple else events_data['EventID'][processed_events].run_number
                if isMC>0 and (conditions_run_number,event_number) in eventList:
                    print(f'List seems to be using the conditions run number instead of MC dataset ID (process ID) - accepting event')
                else:
                    continue

            # Make a key for output - copy what's done by ATLAS
            my_key = f'{event_number}/{run_number}'

            # Start filling the event information
            this_event = { 'event number': event_number,
                           'run number': run_number }

            # Add Jets and Large-radius Jets
            this_event['Jets'] = {}
            for jetCollection in ['jet','largeRJet'] if ntuple else ['AnalysisJets','AnalysisLargeRJets']:
                this_event['Jets'][jetCollection] = []
                # Loop over (small-R) jets in the event
                for ajet,pt in enumerate(events_data[jetCollection][processed_events].pt):
                    # Add the necessary information
                    jet_data = {'eta': events_data[jetCollection][processed_events].eta[ajet],
                                'phi': events_data[jetCollection][processed_events].phi[ajet],
                                'energy': events_data[jetCollection][processed_events].energy[ajet],
                                'coneR': jet_radii[jetCollection]}
                    this_event['Jets'][jetCollection] += [ jet_data ]

            # Add vertices
            if not ntuple:
                this_event['Vertices'] = {}
                this_event['Vertices']['PrimaryVertices'] = []
                # Loop over vertices in the event
                for avert,x in enumerate(events_data['PrimaryVertices'][processed_events].x):
                    vert_data = {'x':events_data['PrimaryVertices'][processed_events].x[avert],
                                 'y':events_data['PrimaryVertices'][processed_events].y[avert],
                                 'z':events_data['PrimaryVertices'][processed_events].z[avert]}
                    this_event['Vertices']['PrimaryVertices'] += [ vert_data ]
                    # We are not linking to tracks here; not all tracks are in PHYSLITE

            # Add MET
            this_event['MissingEnergy'] = {}
            this_event['MissingEnergy']['MET'] = []
            if ntuple:
                this_event['MissingEnergy']['MET'] += [ {'etx': events_data['AnalysisMET'][processed_events].etx,
                                                         'ety': events_data['AnalysisMET'][processed_events].ety} ]
            else:
                for amet,source in enumerate(events_data['AnalysisMET'][processed_events].source):
                    # Check for the magic number that is final MET
                    if source != 69664:
                        continue
                    this_event['MissingEnergy']['MET'] += [ {'etx': events_data['AnalysisMET'][processed_events].etx[amet],
                                                             'ety': events_data['AnalysisMET'][processed_events].ety[amet]} ]

            # Add tracks that we'll want down the line
            this_event['Tracks'] = {}
            if not ntuple:
                # InDetTrackParticles are your vanilla tracks that most objects use
                # MuonSpectrometerTrackParticles are stand-alone muon spectrometer tracks
                # ExtrapolatedMuonTrackParticles are muon tracks extrapolated to the primary vertex
                # GSFTrackParticles are refitted (global sequential fitter) tracks for electrons in particular
                for trkCollection in trkCollectionList:
                    this_event['Tracks'][trkCollection] = []
                    for atrk,phi in enumerate(events_data[trkCollection][processed_events].phi):
                        track_data = {'dparams': [events_data[trkCollection][processed_events].d0[atrk],
                                                  events_data[trkCollection][processed_events].z0[atrk],
                                                  events_data[trkCollection][processed_events].phi[atrk],
                                                  events_data[trkCollection][processed_events].theta[atrk],
                                                  events_data[trkCollection][processed_events].qOverP[atrk] ] }
                        # Info that is only in some of the track collections
                        if trkCollection in trkCollectionsWithNdofChi2:
                            track_data['ndof'] = events_data[trkCollection][processed_events].ndof[atrk]
                            track_data['chi2'] = events_data[trkCollection][processed_events].chi2[atrk]
                        this_event['Tracks'][trkCollection] += [ track_data ]

            # Add e/gamma clusters in preparation for photons and electrons
            this_event['CaloClusters'] = {}
            this_event['CaloClusters']['egammaClusters'] = []
            if not ntuple:
                for aclu,phi in enumerate(events_data['egammaClusters'][processed_events].phi):
                    cluster_data = {'energy': events_data['egammaClusters'][processed_events].energy[aclu],
                                    'eta': events_data['egammaClusters'][processed_events].eta[aclu],
                                    'phi': events_data['egammaClusters'][processed_events].phi[aclu]}
                    this_event['CaloClusters']['egammaClusters'] += [ cluster_data ]

            # Now we can add photons, nice and easy...
            this_event['Photons'] = {}
            this_event['Photons']['AnalysisPhotons'] = []
            for apho,phi in enumerate(events_data['AnalysisPhotons'][processed_events].phi):
                photon_data = {'energy': events_data['AnalysisPhotons'][processed_events].energy[apho],
                               'eta': events_data['AnalysisPhotons'][processed_events].eta[apho],
                               'phi': events_data['AnalysisPhotons'][processed_events].phi[apho]}
                # A photon is allowed to have multiple clusters linked; we want the first index, which 'defines' it
                if ntuple:
                    # If we are running on the flat ntuple, fake the cluster data based on the photon itself
                    cluster_data = {'energy': events_data['AnalysisPhotons'][processed_events].energy[apho],
                                    'eta': events_data['AnalysisPhotons'][processed_events].eta[apho],
                                    'phi': events_data['AnalysisPhotons'][processed_events].phi[apho]}
                    clu_link_index = len(this_event['CaloClusters']['egammaClusters'])
                    this_event['CaloClusters']['egammaClusters'] += [ cluster_data ]
                else:
                    clu_link_index = events_data['AnalysisPhotons'][processed_events].Cluster_Links[apho][0]['m_persIndex']
                if clu_link_index < 1000000:
                    photon_data['LinkedClusters'] = [ f'egammaClusters:{clu_link_index}' ]
                this_event['Photons']['AnalysisPhotons'] += [photon_data]

            # Now we can add muons and electrons
            this_event['Muons'] = {}
            this_event['Muons']['AnalysisMuons'] = []
            this_event['Electrons'] = {}
            this_event['Electrons']['AnalysisElectrons'] = []

            if ntuple:
                this_event['Tracks']['LeptonTracks'] = []
                for alep,leptype in enumerate(events_data['Leptons'][processed_events].lep_type):
                    # Add the track info
                    track_link_index = len(this_event['Tracks']['LeptonTracks'])
                    track_data = {'dparams': [events_data['Leptons'][processed_events].d0[alep],
                                              events_data['Leptons'][processed_events].z0[alep],
                                              events_data['Leptons'][processed_events].phi[alep],
                                              events_data['Leptons'][processed_events].theta[alep], # Conversion thanks to 4-vector
                                              events_data['Leptons'][processed_events].charge[alep] / events_data['Leptons'][processed_events].p[alep] ] }
                    this_event['Tracks']['LeptonTracks'] += [ track_data ]
                    # A bit annoying here that electrons and muons use different formats
                    if leptype==11:
                        electron_data = {'energy': events_data['Leptons'][processed_events].energy[alep],
                                         'eta': events_data['Leptons'][processed_events].eta[alep],
                                         'phi': events_data['Leptons'][processed_events].phi[alep],
                                         'LinkedTracks': [ f'LeptonTracks:{track_link_index}' ]}
                        # If we are running on the flat ntuple, fake the cluster data based on the electron itself
                        cluster_data = {'energy': events_data['Leptons'][processed_events].energy[alep],
                                        'eta': events_data['Leptons'][processed_events].eta[alep],
                                        'phi': events_data['Leptons'][processed_events].phi[alep]}
                        clu_link_index = len(this_event['CaloClusters']['egammaClusters'])
                        this_event['CaloClusters']['egammaClusters'] += [ cluster_data ]
                        electron_data['LinkedClusters'] = [ f'egammaClusters:{clu_link_index}' ]
                        this_event['Electrons']['AnalysisElectrons'] += [ electron_data ]
                    elif leptype==13:
                        muon_data = {'Eta':events_data['Leptons'][processed_events].eta[alep],
                                     'Phi':events_data['Leptons'][processed_events].phi[alep],
                                     'LinkedClusters': None,
                                     'Quality': 'Tight' if events_data['Leptons'][processed_events].tight[alep] else 'Medium',
                                     'LinkedTracks': [ f'LeptonTracks:{track_link_index}' ]}
                        # PassedHighPt, Type not stored in the ntuple
                        this_event['Muons']['AnalysisMuons'] += [ muon_data ]
            else:
                # Loop over all the muons we have
                for amuon,phi in enumerate(events_data['AnalysisMuons'][processed_events].Phi):
                    my_quality = events_data['AnalysisMuons'][processed_events].Quality[amuon]
                    muon_data = {'Eta':events_data['AnalysisMuons'][processed_events].Eta[amuon],
                                 'Phi':events_data['AnalysisMuons'][processed_events].Phi[amuon],
                                 'LinkedClusters': None,
                                 'Type': muon_types[ events_data['AnalysisMuons'][processed_events].Type[amuon] ],
                                 'Quality': muon_quality[ my_quality%8 ],
                                 'PassedHighPt': my_quality>=16,
                                 'LinkedTracks': []}
                    # Add track links. When the link is invalid, the index is very large (4e9)
                    idtp_link_index = events_data['AnalysisMuons'][processed_events].IDTP_Link[amuon]
                    if idtp_link_index < 1000000:
                        muon_data['LinkedTracks'] += [ f'InDetTrackParticles:{idtp_link_index}' ]
                    mstp_link_index = events_data['AnalysisMuons'][processed_events].MSTP_Link[amuon]
                    if mstp_link_index < 1000000:
                        muon_data['LinkedTracks'] += [ f'MuonSpectrometerTrackParticles:{mstp_link_index}' ]
                    emtp_link_index = events_data['AnalysisMuons'][processed_events].EMTP_Link[amuon]
                    if emtp_link_index < 1000000:
                        muon_data['LinkedTracks'] += [ f'ExtrapolatedMuonTrackParticles:{emtp_link_index}' ]
                    # NB if we want to be fancy, we could add a passed high-pT cuts flag as well
                    this_event['Muons']['AnalysisMuons'] += [ muon_data ]
    
                for aelec,phi in enumerate(events_data['AnalysisElectrons'][processed_events].phi):
                    electron_data = {'energy': events_data['AnalysisElectrons'][processed_events].energy[aelec],
                                     'eta': events_data['AnalysisElectrons'][processed_events].eta[aelec],
                                     'phi': events_data['AnalysisElectrons'][processed_events].phi[aelec]}
                    # An electron is allowed to have multiple clusters linked; we want the first index, which 'defines' it
                    clu_link_index = events_data['AnalysisElectrons'][processed_events].Cluster_Links[aelec][0]['m_persIndex']
                    if clu_link_index < 1000000:
                        electron_data['LinkedClusters'] = [ f'egammaClusters:{clu_link_index}' ]
                    track_link_index = events_data['AnalysisElectrons'][processed_events].Track_Links[aelec][0]['m_persIndex']
                    if track_link_index<1000000:
                        electron_data['LinkedTracks'] = [ f'GSFTrackParticles:{track_link_index}' ]
                    this_event['Electrons']['AnalysisElectrons'] += [ electron_data ]

            # Add the event to my output dictionary
            output_dict[my_key] = this_event

            # Did we process enough events?
            if len(output_dict)>=events:
                break
            if len(output_dict)==len(eventList) and len(eventList)>0:
                break

        # Loop over files
        # Did we process enough events?
        if len(output_dict)>=events:
            break
        if len(output_dict)==len(eventList) and len(eventList)>0:
            break

    # Just a little error-checking
    if len(output_dict) < events and len(eventList)==0:
        print(f'Only recorded {len(output_dict)} of {events} requested events')
    if len(eventList)>0 and len(output_dict)!=len(eventList):
        print(f'Of the {len(eventList)} requested run,event number pairs, only found {len(output_dict)} to record')

    # Return the output dictionary that is now full
    return output_dict


# Actual application execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog='write_simple_json.py',
                    description='Writes events in simple json format for the Phoenix event display')
    parser.add_argument('-n', '--nEvents', default=10, type=int, help='Number of events to run over (default 10)')
    parser.add_argument('-s', '--nSkip', default=0, type=int, help='Number of events to skip (default 0)')
    parser.add_argument('-o', '--outputFile', default='my_events.json', help='Name of output json file (default my_events.json)')
    parser.add_argument('-l', '--eventList', default=None, help='File containing run/event number pairs for outputting')
    parser.add_argument('filenames', help='Comma-separated list of files for processing')
    args = parser.parse_args()

    # Get all the files we are meant to go through, looping through the comma-separated list
    all_files = []
    for file_set in args.filenames.split(','):
        # Account for wildcards etc passed by the user
        all_files += glob.glob(file_set)

    print(f'Processing {all_files}')
    print(f'Will write {args.nEvents} events, skipping {args.nSkip}, to {args.outputFile}')
    if args.eventList:
        print(f'Will use event list {args.eventList}')

    # Get our output dictionary
    output_dict = json_format(files=all_files, events=args.nEvents, skip=args.nSkip, eventListFile=args.eventList)

    # Write my json file
    with open(args.outputFile,'w') as json_output_file:
        json.dump(output_dict,json_output_file)

    # Exit cleanly
    sys.exit(0)
