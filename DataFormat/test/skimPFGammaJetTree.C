#include "../interface/GammaJetFitData.h"
#include "pf_gammajettree.h"

void skimPFGammaJetTree(TString inpFileName,
			Long64_t maxEntries=-1,
			TString outFileName="skim.root") {
  pf_gammajettree inpData(inpFileName);

  inpData.DeactivateBranches();
  inpData.ActivateBranches_forFitSkim();

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

    aux.SetEventNo(inpData.EventNumber);
    aux.SetRunNo(inpData.RunNumber);
    aux.SetGenE(inpData.tagPho_genEnergy,inpData.ppfjet_genE);
    if (!aux.SetPhotonQuality(inpData.tagPho_idLoose, inpData.tagPho_idTight)){
      std::cout << "photon quality error detected\n";
      return;
    }
    aux.SetJetQuality(inpData.passTightJetID(1));
    dt->SetAuxInfo(aux);

    dt->SetTagEEtaPhi(inpData.tagPho_energy,
		      inpData.tagPho_eta, inpData.tagPho_phi);
    double ecalE= inpData.getSumEcalE(0,1);
    double hcalE_noRecHits= inpData.getSumHcalE_trackDiffEcal(1);
    dt->SetProbeEtaPhiEn(inpData.ppfjet_eta, inpData.ppfjet_phi,
			 ecalE+hcalE_noRecHits, inpData.getHcalEMap(1,1e-4));

    if (0)
    if (fabs(dt->CalcDiffEt())/dt->GetTagEt() > 0.1) {
      std::cout << "event iEntry=" << iEntry << ", dEt/tagEt= "
		<< dt->CalcDiffEt() << "/" << dt->GetTagEt() << " = "
		<< (dt->CalcDiffEt()/dt->GetTagEt())
		<< " .. skipping\n";
      continue;
    }

    passedCount++;

    tree->Fill();
  }

  tree->Write();
  fout->Close();
  std::cout << "file <" << fout->GetName() << "> created\n";
  std::cout << "nEntries=" << nEntries << ", passedCount=" << passedCount << "\n";
  return;
}
