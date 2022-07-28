<CENTER>
    <br><h1>An introductional notebook to RDataFrame in C++</h1></br>
</CENTER>

<p>In this notebook you can find an easy set of commands that show some basic computing techniques commonly used in High Energy Physics (HEP) analyzes.</p>

<p>It also shows how to create an histogram, fill it and draw it. Moreover it is an introduction to [ROOT](https://root.cern.ch/) too. The final output is a plot with the number of leptons.</p>

<p><b>all done with less that 20 lines of code!</b></p>


The Framework used is [RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html#default-branches), a modern, high-level interface to analyse data stored in TTree, CVS and other data formats, through C++ or Python

<h3>Cell #1</h3>
<p>At first we have to include several helpers that will support our analysis:</p>



#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include <Math/Vector4Dfwd.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include "TStyle.h"
#include <string>
#include <ctime>

%jsroot on


<h3>Cell #3</h3>
<p>Next we have to open the data that we want to analyze. It is stored in a _*.root_ file that consists of a tree (mini) having branches and leaves.</p>

In ROOT to get access and open a root file we have to do the following:

TFile *file = TFile::Open("https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/1largeRjet1lep/MC/mc_361106.Zee.1largeRjet1lep.root"); // 13 TeV sample


<b>get the tree:<b>
    
TTree *tree = (TTree*) file->Get("mini");

Then to set the branches e.g.:
    
UInt_t lepton_n = -1;
tree->SetBranchAddress("lep_n", &lepton_n);
    
Finally, loop over all entries to the Histograms, apply a set of selections ..etc



// You set a list of data files in one line
//ROOT::RDataFrame df("mini",{"https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/1largeRjet1lep/MC/mc_361106.Zee.1largeRjet1lep.root","http://opendata.atlas.cern/release/samples/MC/mc_147770.Zee.root"});


// or one can do
//TChain chain("mini");
//chain.Add("https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/1largeRjet1lep/MC/mc_361106.Zee.1largeRjet1lep.root");
//chain.Add("http://opendata.atlas.cern/release/samples/MC/mc_147770.Zee.root");
//ROOT::RDataFrame df(chain);


<b>Let's do a simple thing <b>


ROOT::RDataFrame df("mini","https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/1largeRjet1lep/MC/mc_361106.Zee.1largeRjet1lep.root");


<b>Get the number of Entries in RDataFrame<b>

auto nEntries = *df.Count();
std::cout<<"Number of Entries="<<nEntries<<std::endl;



<b>get the list of column defined in the dataset<b>



//std::vector<std::string> inputCol = df.GetColumnNames();
auto inputCol = df.GetColumnNames();
for(int i=0;i<inputCol.size();++i){
    std::cout<<inputCol.at(i)<<std::endl;
}



%%time
gStyle->SetOptStat(0); gStyle->SetTextFont(42);
auto c = new TCanvas("c", "", 800, 700);


auto h = df.Define("nlepton","lep_n").Filter("lep_n==3").Histo1D("nlepton");
auto h1 = df.Define("nlepton_cut","lep_n").Filter("lep_n <= 2").Histo1D("nlepton_cut");
h->Draw();
h1->Draw("SAME");


c->Draw();










