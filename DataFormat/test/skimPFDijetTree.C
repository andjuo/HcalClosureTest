#include "../interface/GammaJetFitData.h"
#include "link_PFDijetTree.h"

void skimPFDijetTree(TString inpFileName,
		     Long64_t maxEntries=-1,
		     TString outFileName="skim_dijet.root") {
  dijet_PFNtuple inpData(inpFileName);

  //inpData.DeactivateBranches();
  //inpData.ActivateBranches_forFitSkim();

  //
  GammaJetEvent_t::Class()->IgnoreTObjectStreamer();
  GammaJetEventAuxInfo_t::Class()->IgnoreTObjectStreamer();

  GammaJetEvent_t *dt= new GammaJetEvent_t();
  GammaJetEventAuxInfo_t aux;

  TFile *fout=new TFile(outFileName,"recreate");
  if (!fout->IsOpen()) {
    std::cout << "Failed to create the file <" << outFileName << ">\n";
    return;
  }

  TTree *tree= new TTree("gjet_data","gjet_data");
  tree->Branch("gjet_data","GammaJetEvent_t",&dt);

  // read in the file
  Long64_t nEntries= inpData.fChain->GetEntries();
  if (maxEntries<0) maxEntries=nEntries;
  Long64_t nBytes=0, passedCount=0;

  for (Long64_t iEntry=0; (iEntry<nEntries) && (iEntry<maxEntries);
       iEntry++) {
    if (inpData.LoadTree(iEntry) < 0) break;
    Long64_t nb = inpData.GetEntry(iEntry);
    nBytes += nb;
    if (iEntry%10000==0) std::cout << " ... reading entry " << iEntry << "\n";
    //std::cout << "ientry=" << iEntry << "\n";

    passedCount++;

    aux.SetEventNo(-1);
    aux.SetRunNo(-1);
    aux.SetGenE(inpData.tpfjet_genE,inpData.ppfjet_genE);
    dt->SetAuxInfo(aux);

    // Note: the dijet event has the tag jet (idx=1) and the probe jet (0)

    double tagjet_ecalE= inpData.getSumEcalE(1,1);
    double tagjet_hcalE_noRecHits= inpData.getSumHcalE_trackDiffEcal(1);
    dt->SetTagEtaPhiEn(inpData.tpfjet_eta, inpData.tpfjet_phi,
		       tagjet_ecalE + tagjet_hcalE_noRecHits,
		       inpData.getHcalEMap(1,1e-4));

    double ecalE= inpData.getSumEcalE(0,1);
    double hcalE_noRecHits= inpData.getSumHcalE_trackDiffEcal(0);
    dt->SetProbeEtaPhiEn(inpData.ppfjet_eta, inpData.ppfjet_phi,
			 ecalE+hcalE_noRecHits,
			 inpData.getHcalEMap(0,1e-4));
    tree->Fill();
  }

  tree->Write();
  fout->Close();
  std::cout << "file <" << fout->GetName() << "> created\n";
  std::cout << "nEntries=" << nEntries << ", passedCount=" << passedCount << "\n";
  return;
}
