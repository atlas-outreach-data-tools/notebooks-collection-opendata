//#define helperFunctions_cxx

#include "TSystem.h"
#include "TInterpreter.h"
#include <ROOT/RVec.hxx>

using VecF_t = const ROOT::RVec<float>&;
using VecD_t = const ROOT::RVec<double>&;
using VecI_t = const ROOT::RVec<int>&;
using VecUI_t = const ROOT::RVec<UInt_t>&;
using VecB_t = const ROOT::VecOps::RVec<bool>;

#include "helperFunctions.h"

std::pair <double,double> getLeptonsFromZ(VecI_t chlep, VecI_t& fllep, VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e){
  double diff = 10000000000.0;
  double diff2 = 10000000000.0;
  
  int Zlep1 = -99;
  int Zlep2 = -99;
  int Zlep3 = -999;
  int Zlep4 = -999;
  double Zmass = -1.0;
  double Zmass2 = -1.0;
  bool foundSFOS = false;
  bool foundSFOS2 =false;
  std::pair <double,double> masses;
  for(unsigned int i=0; i<chlep.size(); i++)
    {
      for(unsigned int j=i+1; j<chlep.size(); j++)
        {
          //Opposite-Sign
          if(chlep[i]*chlep[j]<0)
            {
              //Same-Flavor
              if(abs(fllep[i])==abs(fllep[j]))
                {
                  TLorentzVector p1;
                  p1.SetPtEtaPhiE(pt[i], eta[i], phi[i], e[i]);
                  TLorentzVector p2;
                  p2.SetPtEtaPhiE(pt[j], eta[j], phi[j], e[j]);
                  double mass = (p1+p2).M();
                  /*
                    if(mass<5000){
                    continue; // Filtering out the mll<5GeV
                    }
                  */
                  double massdiff = fabs(mass-91187.6);
                  if(massdiff<diff)
                    {
                      diff=massdiff;
                      Zmass=mass;
                      Zlep1 = i;
                      Zlep2 = j;
                      foundSFOS = true;
                    }
                }
            }
        }

    }
  
  if(foundSFOS){
    for(unsigned int i=0; i<chlep.size(); i++)
      {if(i==Zlep1){
        continue;
      }
      if(i==Zlep2){
        continue;
      }
      for(unsigned int j=i+1; j<chlep.size(); j++){
        if(j==Zlep1){
          continue;
        }
        if(j==Zlep2){
          continue;
        }
        if(chlep[i]*chlep[j]<0){
          if(abs(fllep[i])==abs(fllep[j])){
            TLorentzVector p3;
            p3.SetPtEtaPhiE(pt[i], eta[i], phi[i], e[i]);
            TLorentzVector p4;
            p4.SetPtEtaPhiE(pt[j], eta[j], phi[j], e[j]);
            double mass2=(p3+p4).M();
            /*
            if(mass2<5000){
              continue; //To filter out the mll<5GeV
            }
            */
            double massdiff2 = fabs(mass2-91187.6);
            if(massdiff2<diff2){
              diff2=massdiff2;
              Zmass2=mass2;
              Zlep3=i;
              Zlep4=j;
              foundSFOS2=true;
            }
      
          }
        }


      }}
  }
  masses = std::make_pair(Zmass,Zmass2);
    
  return masses;
}

 ROOT::RVec<int> getLeptonsPairsFromZ(VecI_t chlep, VecI_t& fllep, VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e){
  double diff = 10000000000.0;
  double diff2 = 10000000000.0;
  
  int Zlep1 = -99;
  int Zlep2 = -99;
  int Zlep3 = -999;
  int Zlep4 = -999;
  double Zmass = -1.0;
  double Zmass2 = -1.0;
  bool foundSFOS = false;
  bool foundSFOS2 =false;
  std::pair <double,double> masses;
  for(unsigned int i=0; i<chlep.size(); i++)
    {
      for(unsigned int j=i+1; j<chlep.size(); j++)
        {
          //Opposite-Sign
          if(chlep[i]*chlep[j]<0)
            {
              //Same-Flavor
              if(abs(fllep[i])==abs(fllep[j]))
                {
                  TLorentzVector p1;
                  p1.SetPtEtaPhiE(pt[i], eta[i], phi[i], e[i]);
                  TLorentzVector p2;
                  p2.SetPtEtaPhiE(pt[j], eta[j], phi[j], e[j]);
                  double mass = (p1+p2).M();
                  /*
                    if(mass<5000){
                    continue; // Filtering out the mll<5GeV
                    }
                  */
                  double massdiff = fabs(mass-91187.6);
                  if(massdiff<diff)
                    {
                      diff=massdiff;
                      Zmass=mass;
                      Zlep1 = i;
                      Zlep2 = j;
                      foundSFOS = true;
                    }
                }
            }
        }

    }
  
  if(foundSFOS){
    for(unsigned int i=0; i<chlep.size(); i++)
      {
        if(i==Zlep1){
          continue;
        }
        if(i==Zlep2){
          continue;
        }
        for(unsigned int j=i+1; j<chlep.size(); j++){
          if(j==Zlep1){
            continue;
          }
          if(j==Zlep2){
            continue;
          }
          if(chlep[i]*chlep[j]<0){
            if(abs(fllep[i])==abs(fllep[j])){
              TLorentzVector p3;
              p3.SetPtEtaPhiE(pt[i], eta[i], phi[i], e[i]);
              TLorentzVector p4;
              p4.SetPtEtaPhiE(pt[j], eta[j], phi[j], e[j]);
              double mass2=(p3+p4).M();
              /*
                if(mass2<5000){
                continue; //To filter out the mll<5GeV
                }
              */
              double massdiff2 = fabs(mass2-91187.6);
              if(massdiff2<diff2){
                diff2=massdiff2;
                Zmass2=mass2;
                Zlep3=i;
                Zlep4=j;
                foundSFOS2=true;
              }
      
            }
          }


        }
      }
  }
  ROOT::RVec<int> pairlist;
  pairlist.push_back(Zlep1);
  pairlist.push_back(Zlep2);
  pairlist.push_back(Zlep3);
  pairlist.push_back(Zlep4);
  return pairlist;
}

float deltaR_ll_4(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, int idx1, int idx2) {

  const auto size = int(pt.size());
  if(idx1 > size || idx2 > size){
    printf("deltaR_ll::ERROR \t Indices %i and %i are higher than size of vector %i\n",idx1,idx2,size);
    return -1;
  }

  TLorentzVector p1;
  TLorentzVector p2;
  p1.SetPtEtaPhiE(pt[idx1], eta[idx1], phi[idx1], e[idx1]);
  p2.SetPtEtaPhiE(pt[idx2], eta[idx2], phi[idx2], e[idx2]);
  //std::cout<<"deltaR"<<p1.DeltaR(p2)<<std::endl;
  return p1.DeltaR(p2);
}

float getInvariantMass_ll_4(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, int idx1, int idx2){
  TLorentzVector p1;
	p1.SetPtEtaPhiE(pt[idx1], eta[idx1], phi[idx1], e[idx1]);
  TLorentzVector p2;
  p2.SetPtEtaPhiE(pt[idx2], eta[idx2], phi[idx2], e[idx2]);
  double mass = (p1+p2).M();
  return mass;
}

int lepchannel(VecI_t& fllep){
  int totfl = 0;
  const auto size = fllep.size();
  if(size != 4){
  std::cout<<"size of lepton is "<<size<<std::endl;
  }
  for (size_t i=0; i < size; ++i) {
    totfl += fllep[i];
  }
  return totfl;
}







