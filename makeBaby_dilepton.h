#ifndef ScanChain_h
#define ScanChain_h

// C++ includes
#include <string>
#include <vector>

// ROOT includes
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "Math/VectorUtil.h"
#include "Math/Vector4D.h"

#include "/home/users/sanil/CORE/ssSelections.h"
#include "/home/users/sanil/CORE/muonSelections.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class babyMaker {

 public:

  babyMaker() {};
  ~babyMaker() {
    delete BabyFile_;
    delete BabyTree_;
  };

  void ScanChain(TChain* chain, std::string baby_name = "testSample", unsigned int numEvents = 0, float customScale = -1);

  void MakeBabyNtuple(const char *);
  void InitBabyNtuple();
  void FillBabyNtuple();
  void CloseBabyNtuple();

 private:

  TFile *BabyFile_;
  TTree *BabyTree_;

  //baby ntuple variables
  float pf_sumet;
  float pf_tcmet;
  float pf_tcmetPhi;


  float met;
  float metPhi;

  std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > gen_p4;
  std::vector<int> gen_id;
  std::vector<int> gen_idx;
  std::vector<int> gen_id_mother;
  std::vector<vector<int> > gen_lepdaughter_id; 

  //std::vector<LorentzVector> jets_p4;
  std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jets_p4;
  std::vector<float> jets_p4Correction;
  std::vector<int> jets_chargedMulti;
  std::vector<int> jets_neutralMulti;
  std::vector<float> jets_pt;


  int type;

  LorentzVector ll_p4;
  LorentzVector lt_p4;
  LorentzVector total_p4;


  int ll_id;
  int lt_id;
  int ll_charge;
  int lt_charge;
  int ll_index;
  int lt_index;


  bool isSamesignLep;
  bool isStopLep;
  bool isZmetLep;

  bool trackingProblemVeto;
  bool tauVeto;
  bool isoTrackVeto;

  int isRealData;
  int nvtxs;
  int njets;
  int nbtag;

  float scale_1fb;
  std::vector<float> btagDiscriminant;
  
  int numEvents;

  int eventNumber;
  int runNumber;
  int lumiBlock;

bool ele17_ele8;
bool mu17_mu8;  
bool mu17_Tkmu8;
bool mu17_ele8; 
bool mu8_ele17; 

int ll_isFromW;
int lt_isFromW;
  string file;

};

#endif

void babyMaker::MakeBabyNtuple(const char *BabyFilename){

  //
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  BabyFile_ = new TFile(Form("%s", BabyFilename), "RECREATE");
  BabyFile_->cd();
  BabyTree_ = new TTree("tree", "A Baby Ntuple");

  BabyTree_->Branch("pf_sumet", &pf_sumet );
  BabyTree_->Branch("pf_tcmet", &pf_tcmet );
  BabyTree_->Branch("pf_tcmetPhi", &pf_tcmetPhi );


  BabyTree_->Branch("gen_p4", &gen_p4 );
  BabyTree_->Branch("gen_id", &gen_id );
  BabyTree_->Branch("gen_idx", &gen_idx );
  BabyTree_->Branch("gen_id_mother", &gen_id_mother );
  BabyTree_->Branch("gen_lepdaughter_id", &gen_lepdaughter_id );

  BabyTree_->Branch("met", &met );
  BabyTree_->Branch("metPhi", &metPhi );
  BabyTree_->Branch("jets_p4", &jets_p4 );
  BabyTree_->Branch("jets_pt", &jets_pt );
  BabyTree_->Branch("jets_chargedMulti", &jets_chargedMulti );
  BabyTree_->Branch("jets_neutralMulti", &jets_neutralMulti );
  BabyTree_->Branch("jets_p4Correction", &jets_p4Correction );
  BabyTree_->Branch("type", &type);

  BabyTree_->Branch("ll_p4", &ll_p4);
  BabyTree_->Branch("lt_p4", &lt_p4);
  BabyTree_->Branch("total_p4", &total_p4);

  BabyTree_->Branch("ll_id", &ll_id);
  BabyTree_->Branch("lt_id", &lt_id);

  BabyTree_->Branch("ll_charge" , &ll_charge);
  BabyTree_->Branch("lt_charge", &lt_charge);
  BabyTree_->Branch("ll_index", &ll_index);
  BabyTree_->Branch("lt_index", &lt_index);
  
  BabyTree_->Branch("eventNumber", &eventNumber);
  BabyTree_->Branch("runNumber", &runNumber);
  BabyTree_->Branch("lumiBlock", &lumiBlock);
  
  BabyTree_->Branch("scale_1fb", &scale_1fb);
  BabyTree_->Branch("btagDiscriminant", &btagDiscriminant);

  BabyTree_->Branch("isSamesignLep", &isSamesignLep);
  BabyTree_->Branch("isStopLep", &isStopLep);
  BabyTree_->Branch("isZmetLep", &isZmetLep);

  BabyTree_->Branch("trackingProblemVeto", &trackingProblemVeto);
  BabyTree_->Branch("tauVeto", &tauVeto);
  BabyTree_->Branch("isoTrackVeto", &isoTrackVeto);

  BabyTree_->Branch("isRealData", &isRealData );
  BabyTree_->Branch("nvtxs", &nvtxs );
  BabyTree_->Branch("njets", &njets);
  BabyTree_->Branch("nbtag", &nbtag);

  BabyTree_->Branch("numEvents", &numEvents);

  BabyTree_->Branch("file", &file);
  BabyTree_->Branch("ele17_ele8", &ele17_ele8);
  BabyTree_->Branch("mu17_mu8", &mu17_mu8  );
  BabyTree_->Branch("mu17_Tkmu8", &mu17_Tkmu8  );
  BabyTree_->Branch("mu17_ele8", &mu17_ele8 );
  BabyTree_->Branch("mu8_ele17", &mu8_ele17 );
  BabyTree_->Branch("ll_isFromW", &ll_isFromW);
  BabyTree_->Branch("lt_isFromW", &lt_isFromW);

  return;
}

void babyMaker::InitBabyNtuple () {

  pf_sumet = -999.0;
  pf_tcmet = -999.0;
  pf_tcmetPhi = -999.0;

  met = -999.0;
  metPhi = -999.0;

  type = -1;

  ll_id = -1;
  lt_id = -1;

  ll_charge = -999;
  lt_charge = -999;
  ll_index = -1;
  lt_index = -1;

  scale_1fb = 1;

  isRealData = -999;
  nvtxs = -1;
  njets = -1;
  nbtag = -1;

  eventNumber = -1;
  runNumber = -1;
  lumiBlock = 1;

  numEvents = 0;
ll_isFromW = -99;
lt_isFromW = -99;

  return;
}

void babyMaker::FillBabyNtuple(){
  BabyTree_->Fill();
  return;
}

void babyMaker::CloseBabyNtuple(){
  BabyFile_->cd();
  BabyTree_->Write();
  BabyFile_->Close();
  return;
}

bool passHBHEFilter()
{
  if (!cms2.evt_isRealData())                          {return true; }
  if (!cms2.evt_hbheFilter())                          {return false;}
  if (cms2.hcalnoise_isolatedNoiseSumE() >= 50.0)      {return false;}
  if (cms2.hcalnoise_isolatedNoiseSumEt() >= 25.0)     {return false;}
  if (cms2.hcalnoise_numIsolatedNoiseChannels() >= 10) {return false;}
  return true;
}

 
bool goodTau(int index){

  if(cms2.taus_pf_p4()[index].pt() < 20) return false;
  if(fabs(cms2.taus_pf_p4()[index].eta()) > 2.3) return false;
  if(cms2.taus_pf_byDecayModeFinding()[index] < 0.5) return false;
  if(cms2.taus_pf_againstElectronLoose()[index] < 0.5) return false;
  if(cms2.taus_pf_againstMuonTight2()[index] < 0.5) return false;
  if(cms2.taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits()[index] < 0.5) return false;

  return true;
}



 
