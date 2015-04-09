// C++
#include <iostream>
#include <vector>
#include <set>
#include <fstream>

// ROOT
#include "TDirectory.h"
#include "TTreeCache.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"

// CMS2
#include "CMS2.h"
#include "electronSelections.h"
#include "muonSelections.h"
#include "ssSelections.h"
#include "trackSelections.h"
#include "eventSelections.h"
#include "susySelections.h"
#include "susySelections.h"
#include "mcSelections.h"

#include "/home/users/sanil/CMSSW_5_3_2_patch4_V05-03-23/src/CMS2/NtupleMacros/Tools/goodrun.cc"
#include "/home/users/sanil/dilepton/Include.C"

// header
#include "makeBaby_dilepton.h"
#include "/home/users/sanil/CORE/lostlepton/IsoTrackVeto.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std;
using namespace tas;
using namespace ROOT::Math::VectorUtil;
using namespace samesign;


void babyMaker::ScanChain(TChain* chain, std::string baby_name, unsigned int numEvent, float customScale){

  MakeBabyNtuple( Form("/home/users/sanil/makebaby/babies/%s.root", baby_name.c_str()) );

  // File Loop
  unsigned int fileCounter = 0;
  int nDuplicates = 0;
  int nEvents = chain->GetEntries();
  unsigned int nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;



  // set good run list
  set_goodrun_file("/home/users/jgran/analysis/sswh/fakes/json/final_19p49fb_cms2.txt");

  //File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    
    if (fileCounter > numEvent && numEvent !=0) {
        break;
    }
    // Get File Content
    TFile f( currentFile->GetTitle() );
    TTree *tree = (TTree*)f.Get("Events");
    TTreeCache::SetLearnEntries(10);
    tree->SetCacheSize(128*1024*1024);
    cms2.Init(tree);

    // event Loop
    unsigned int nEventsTree = tree->GetEntriesFast();

    for(unsigned int event = 0; event < nEventsTree; ++event) {

      // get event content
      tree->LoadTree(event);
      cms2.GetEntry(event);
      ++nEventsTotal;

       // select good 
      if(evt_isRealData() && !goodrun(evt_run(), evt_lumiBlock())) continue;

      if(evt_isRealData()){
        DorkyEventIdentifier id = {evt_run(), evt_event(), evt_lumiBlock()};
        if (is_duplicate(id)){
	nDuplicates++;
          continue;
        }
      }
  
    // Progress
      CMS2::progress( nEventsTotal, nEventsChain );

      // count number of hypotheses
      int index = 0;



      // Event Cleaning 
      // Filters

      /// Rho
      if ( evt_ww_rho_vor() > 40 )  continue;

      if( evt_isRealData()){
	if(  evt_cscTightHaloId()    )		  continue;
	if( !evt_hbheFilter()        )		  continue;
	if( !filt_hcalLaser()        )		  continue;
	if( !filt_ecalTP()           )		  continue;
	if( !filt_trackingFailure()  )		  continue;
	if( !filt_eeBadSc()          )		  continue;
	if( !passHBHEFilter()	   )		  continue;
      }
	//end Filters


        bool 	trkProbVeto = false, tau_veto = false, foundIsoTrack = false,
		ele17_ele8_ = false, mu17_mu8_   = false, mu17_Tkmu8_ = false, 
		mu17_ele8_  = false, mu8_ele17_  = false;

	int  ll_isFromW_ = 0, lt_isFromW_ = 0;





//Dilepton Selection
	float pt_max = 0;
	int n_good = -1;
for(unsigned int i=0; i<hyp_p4().size(); i++) {

	int lt_id_  = abs(hyp_lt_id().at(i)), lt_idx = hyp_lt_index().at(i), 
	    ll_id_  = abs(hyp_ll_id().at(i)), ll_idx = hyp_ll_index().at(i); 
// opposite sign
	if(hyp_ll_charge().at(i)*hyp_lt_charge().at(i) > 0) continue;
// eta
	if(abs(hyp_ll_p4().at(i).eta()) > 2.4)    continue;
	if(abs(hyp_lt_p4().at(i).eta()) > 2.4)    continue;        
// pt
	if(hyp_ll_p4().at(i).pt()  < 20.0)   continue;        
	if(hyp_lt_p4().at(i).pt()  < 20.0)   continue;        
//electron selections
	if (ll_id_ == 11 && !passElectronSelection_Stop2012_v3(ll_idx) )  continue;
	if (lt_id_ == 11 && !passElectronSelection_Stop2012_v3(lt_idx) )  continue;
//muon selections
	if (ll_id_ == 13 && !muonId(ll_idx, ZMet2012_v1))  continue;
	if (lt_id_ == 13 && !muonId(lt_idx, ZMet2012_v1))  continue;

//find the highest total pt pair
	float pt_sum = hyp_ll_p4().at(i).pt() + hyp_lt_p4().at(i).pt();
	if (pt_sum > pt_max) {
		pt_max = pt_sum;
		n_good = i;
		}//pt_sum
	}//hyp_p4().size

/* if no dilepton pair passed the above selections go to the next event */
	if(n_good == -1) continue;

//JETS
      int n_btag = 0;
      int n_jets = 0;

for(unsigned int c = 0; c < pfjets_p4().size(); c++) {

	float _jetPt = pfjets_p4().at(c).pt() * pfjets_corL1FastL2L3().at(c);

	if(_jetPt < 40) continue;       

	if(abs(pfjets_p4().at(c).eta()) > 3) continue;

	float dr_ll = DeltaR(pfjets_p4().at(c), hyp_ll_p4().at(n_good) );
	float dr_lt = DeltaR(pfjets_p4().at(c), hyp_lt_p4().at(n_good) );

	if(dr_ll < 0.4) continue;      
	if(dr_lt < 0.4) continue;      
	n_jets++;

	//b Tagging
	float bTag_vertex = pfjets_combinedSecondaryVertexBJetTag().at(c);
	if (  bTag_vertex > 0.244) n_btag++;

	//tracking problem veto
	if (( abs(pfjets_p4().at(c).Eta()) > 0.9 && abs(pfjets_p4().at(c).Eta()) < 1.9 && (pfjets_chargedMultiplicity().at(c) - pfjets_neutralMultiplicity().at(c)) > 40) ) {
	trkProbVeto = true;	       
		}

     }//pfjets_p4().size()

//Event Clean-Up//

// Tau Veto //
if(evt_isRealData()){
for(unsigned int i = 0; i < taus_pf_p4().size(); i++){
if(goodTau(i)) tau_veto = true;
      }
}
//if tauVeto == 1 continue


      ///// Iso Track Veto /////

      for (unsigned int ipf = 0; ipf < pfcands_p4().size(); ipf++) {

        if(cms2.pfcands_charge().at(ipf) == 0) continue;
        const bool isLepton = (abs(cms2.pfcands_particleId().at(ipf))==11) || (abs(cms2.pfcands_particleId().at(ipf))==13);
        const float cand_pt = cms2.pfcands_p4().at(ipf).pt();
        if(cand_pt < 5) continue;
        if(!isLepton && (cand_pt < 10)) continue;
      	float dr_lt = DeltaR( cms2.pfcands_p4().at(ipf) , hyp_lt_p4().at(n_good) );
      	float dr_ll = DeltaR( cms2.pfcands_p4().at(ipf) , hyp_ll_p4().at(n_good) );
	if( isLepton && dr_ll < 0.1 )   continue;
	if( isLepton && dr_lt < 0.1 )   continue;
        int itrk = -1;
        float mindz = 999.;
        if (abs(cms2.pfcands_particleId().at(ipf))!=11) {
          itrk = cms2.pfcands_trkidx().at(ipf);
          if( itrk >= (int)cms2.trks_trk_p4().size() || itrk < 0 ) continue;
          mindz=trks_dz_pv(itrk,0).first;
        }   
        if (abs(cms2.pfcands_particleId().at(ipf))==11 && cms2.pfcands_pfelsidx().at(ipf)>=0) {
          itrk = cms2.els_gsftrkidx().at(cms2.pfcands_pfelsidx().at(ipf));
          if( itrk >= (int)cms2.gsftrks_p4().size() || itrk < 0 ) continue;
          mindz=gsftrks_dz_pv(itrk,0).first;
        }
        if(mindz > 0.1) continue;
        float reliso  = TrackIso(ipf) / cand_pt;
        if(isLepton && (reliso < 0.2)){
          foundIsoTrack = true;
          break;
        }
        if(!isLepton && (reliso < 0.1)){
          foundIsoTrack = true;
          break;
        }
      }
//if isotrackveto == 1 continue
      ///// Iso Track Veto /////


//checks the mother of the ll and lt leptons
//returns integer, check mcSelections for reference
if(!evt_isRealData()){
	ll_isFromW_ = leptonIsFromW(hyp_ll_index().at(n_good), hyp_ll_id().at(n_good),false);
	lt_isFromW_ = leptonIsFromW(hyp_lt_index().at(n_good), hyp_lt_id().at(n_good),false);
}




//// Trigger ///
	ele17_ele8_ = 	passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") 	? 1 : 0;
	mu17_mu8_ = 	passUnprescaledHLTTriggerPattern("HLT_Mu17_Mu8_v") 										? 1 : 0; 	
	mu17_Tkmu8_ = 	passUnprescaledHLTTriggerPattern("HLT_Mu17_TkMu8_v") 										? 1 : 0;
	mu17_ele8_ =	passUnprescaledHLTTriggerPattern("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") 						? 1 : 0;
       	mu8_ele17_ = 	passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") 						? 1 : 0;	


// create the ntuple
	InitBabyNtuple();

// set the branch values

// lepton branches
	type  = hyp_type().at(n_good);
	lt_id = hyp_lt_id().at(n_good);
	ll_id = hyp_ll_id().at(n_good);
	ll_p4 = hyp_ll_p4().at(n_good);
	lt_p4 = hyp_lt_p4().at(n_good);
	ll_charge = hyp_ll_charge().at(n_good);
	lt_charge = hyp_lt_charge().at(n_good);


//JETS branches
	njets = n_jets;
	nbtag = n_btag;
	jets_p4 = pfjets_p4();
	jets_p4Correction = pfjets_corL1FastL2L3();
	btagDiscriminant = pfjets_combinedSecondaryVertexBJetTag();
	jets_chargedMulti = pfjets_chargedMultiplicity();
	jets_neutralMulti = pfjets_neutralMultiplicity();

//missing energy branches
	met = evt_pfmet_type1cor();
	metPhi = evt_pfmetPhi_type1cor();
	pf_sumet = evt_pfsumet(); 
	pf_tcmet = evt_pf_tcmet();  
	pf_tcmetPhi = evt_pf_tcmetPhi(); 
      
//clean-up bools 
	trackingProblemVeto = trkProbVeto;
	tauVeto = tau_veto;
	isoTrackVeto = foundIsoTrack;

//event information
	eventNumber = evt_event();
	runNumber = evt_run();
	lumiBlock = evt_lumiBlock();
	scale_1fb = evt_scale1fb();
	isRealData = evt_isRealData();
	nvtxs = evt_nvtxs();

//triggers
	ele17_ele8 = ele17_ele8_;
	mu17_mu8   = mu17_mu8_ ;
	mu17_Tkmu8 = mu17_Tkmu8_ ;
	mu17_ele8  = mu17_ele8_ ;
	mu8_ele17  = mu8_ele17_ ;

//MC branches
if(!evt_isRealData()){
	ll_isFromW = ll_isFromW_;
	lt_isFromW = lt_isFromW_;
	gen_p4		   = genps_p4();
	gen_id		   = genps_id();
	gen_id_mother	   = genps_id_mother();
	gen_lepdaughter_id = genps_lepdaughter_id();
}
//extra branches

	file = Form("%s", currentFile->GetTitle());

	FillBabyNtuple();

    } // end of loop over events in file

    delete tree;
    f.Close();

  } // end of loop over files

  if (nEventsChain != nEventsTotal){
      cout << "WARNING: The number of events added is not equal to the total number of events!" << endl;
  }

  cout << nDuplicates << " duplicate events were skipped." << endl;

  CloseBabyNtuple();


  return;
}

