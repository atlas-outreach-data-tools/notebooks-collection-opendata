#ifndef helperFunctions_h
#define helperFunctions_h
//using namespace Root

#include "ROOT/RDF/RInterface.hxx"
#include <ROOT/RDataSource.hxx>
#include <ROOT/RCsvDS.hxx>
#include <TROOT.h>
#include <TChain.h>
#include <TObjString.h>
#include <TFile.h>
#include <TString.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
//#include <ifstream>
#include <vector>
#include <time.h>
#include <map>
#include <ctime>
#include <dirent.h>
#include "TLorentzVector.h"
#include "TParameter.h"
#include "TSystem.h"
#include <ROOT/RVec.hxx>
#include <bits/stdc++.h> 
using VecF_t = const ROOT::RVec<float>&;
using VecI_t = const ROOT::RVec<int>&;
using VecB_t = const ROOT::VecOps::RVec<bool>;
using VecD_t = const ROOT::RVec<double>&;
using VecUI_t = const ROOT::RVec<UInt_t>&;
std::pair <double,double> getLeptonsFromZ(VecI_t chlep, VecI_t& fllep, VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e);
ROOT::RVec<int> getLeptonsPairsFromZ(VecI_t chlep, VecI_t& fllep, VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e);
float deltaR_ll_4(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, int idx1, int idx2);
float getInvariantMass_ll_4(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, int idx1, int idx2);
int lepchannel(VecI_t& fllep);
#endif
