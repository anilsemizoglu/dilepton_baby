doBaby_mctt(int numEvents = 0, std::string file_name = "file_name")
{
  gSystem->AddIncludePath(Form("-I%s/CORE", gSystem->Getenv("HOME")));
  gSystem->Load(Form("%s/CORE/libCMS2NtupleMacrosCORE.so", gSystem->Getenv("HOME")));
  gROOT->ProcessLine(".L IsoTrackVeto.cc+");
  gROOT->ProcessLine(".L makeBaby_dilepton.C+");

  babyMaker *baby = new babyMaker();
  
  TChain *tt = new TChain("Events"); 
for(unsigned int i=1;i<=68;i++)
tt->Add(Form("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_%d.root",i));
//skipping merged_ntuple_69, corrupted file
for(unsigned int i=70;i<=139;i++)
tt->Add(Form("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_%d.root",i));

  baby->ScanChain(tt, "ttbar_v4", numEvents); 
 
}
