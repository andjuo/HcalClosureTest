// Oct 10, 2014. DijetRespCorrDatum::CandTrack* info is not saved

#include "DijetRespCorrDataExtended.h"
#include "dijet_PFNtuple.h"
#include "helper.hh"

void skimmer(Long64_t maxEntries=10) {

  //dijet_PFNtuple inpData("dijet_tree.root");
  dijet_PFNtuple inpData("dijet_tree_5809400A.root");
  inpData.DeactivateBranches();

  std::vector<TString> tagV;
  tagV.reserve(30);

  tagV.push_back("tpfjet_E");
  tagV.push_back("tpfjet_p");
  tagV.push_back("tpfjet_eta");
  tagV.push_back("tpfjet_phi");
  tagV.push_back("tpfjet_unkown_E");
  tagV.push_back("tpfjet_electron_E");
  tagV.push_back("tpfjet_muon_E");
  tagV.push_back("tpfjet_photon_E");
  tagV.push_back("tpfjet_had_n");
  tagV.push_back("tpfjet_had_E");
  tagV.push_back("tpfjet_had_EcalE");
  tagV.push_back("tpfjet_had_id");
  tagV.push_back("tpfjet_had_candtrackind");
  tagV.push_back("tpfjet_had_ntwrs");
  tagV.push_back("tpfjet_candtrack_px");
  tagV.push_back("tpfjet_candtrack_py");
  tagV.push_back("tpfjet_candtrack_pz");
  tagV.push_back("tpfjet_twr_hade");
  tagV.push_back("tpfjet_twr_frac");
  tagV.push_back("tpfjet_twr_ieta");
  tagV.push_back("tpfjet_twr_iphi");
  tagV.push_back("tpfjet_twr_clusterind");
  tagV.push_back("tpfjet_cluster_dR");
  tagV.push_back("tpfjet_genE");
  tagV.push_back("tpfjet_genpt");
  
  tagV.push_back("tpfjet_cluster_eta");
  tagV.push_back("tpfjet_cluster_phi");
  tagV.push_back("tpfjet_twr_ieta");
  tagV.push_back("tpfjet_twr_iphi");

  inpData.ActivateBranches(tagV);

  for (unsigned int i=0; i<tagV.size(); ++i) {
    tagV[i].ReplaceAll("tpfjet","ppfjet");
  }
  inpData.ActivateBranches(tagV);

  inpData.ActivateBranches(2,"pf_thirdjet_px","pf_thirdjet_py");
  inpData.ActivateBranches(3,"pf_dijet_deta","pf_dijet_dphi","pf_weight");


  Long64_t nEntries= inpData.fChain->GetEntries();
  if (maxEntries<0) maxEntries=nEntries;

  DijetRespCorrDatumExtended_t *datum= new DijetRespCorrDatumExtended_t();
  DijetRespCorrDatum::Class()->IgnoreTObjectStreamer();
  DijetRespCorrDatumExtended_t::Class()->IgnoreTObjectStreamer();

  TFile outFile("skim.root","RECREATE");
  TTree *outTree = new TTree("dijet_data","dijet_data");
  outTree->Branch("dijet_data","DijetRespCorrDatumExtended_t",&datum);

  Long64_t nBytes=0, passedCount=0;
  for (Long64_t iEntry=0; (iEntry<nEntries) && (iEntry<maxEntries);
       iEntry++) {
    if (inpData.LoadTree(iEntry) < 0) break;
    Long64_t nb = inpData.GetEntry(iEntry);
    nBytes += nb;
    if (iEntry%10000==0) std::cout << " ... reading entry " << iEntry << "\n";

    //HERE("got"); std::cout << inpData << "\n";

    datum->Clear();

    datum->SetOrigIEntry(iEntry);
    datum->SetWeight(inpData.pf_weight);
    datum->SetThirdJetPx(inpData.pf_thirdjet_px);
    datum->SetThirdJetPy(inpData.pf_thirdjet_py);
    datum->SetCandTrackN(0);
    datum->SetDijetDEtaDPhi(inpData.pf_dijet_deta,inpData.pf_dijet_dphi);
    datum->SetJetEPEtaPhi(1,
			  inpData.tpfjet_E,inpData.tpfjet_p,
			  inpData.tpfjet_eta,inpData.tpfjet_phi);
    datum->SetJetEPEtaPhi(0,
			  inpData.ppfjet_E,inpData.ppfjet_p,
			  inpData.ppfjet_eta,inpData.ppfjet_phi);

    double tjet_nonHadEn=
      inpData.tpfjet_unkown_E + inpData.tpfjet_electron_E +
      inpData.tpfjet_muon_E + inpData.tpfjet_photon_E;
    double pjet_nonHadEn=
      inpData.ppfjet_unkown_E + inpData.ppfjet_electron_E +
      inpData.ppfjet_muon_E + inpData.ppfjet_photon_E;

    if (0 &&  !datum->PassCuts()) continue;
    //else std::cout << "passed iEntry=" << iEntry << "\n";
    passedCount++;

    datum->SetJetTowers(1,
			tjet_nonHadEn,
			*inpData.tpfjet_twr_ieta,
			*inpData.tpfjet_twr_iphi,
			*inpData.tpfjet_twr_hade,
			*inpData.tpfjet_twr_frac,
			*inpData.tpfjet_had_EcalE,
			*inpData.tpfjet_twr_clusterind,
			*inpData.tpfjet_cluster_dR,
			*inpData.tpfjet_had_id,
			*inpData.tpfjet_had_ntwrs,
			*inpData.tpfjet_had_candtrackind,
			*inpData.tpfjet_candtrack_px,
			*inpData.tpfjet_candtrack_py,
			*inpData.tpfjet_candtrack_pz);
    datum->SetJetTowers(0,
			pjet_nonHadEn,
			*inpData.ppfjet_twr_ieta,
			*inpData.ppfjet_twr_iphi,
			*inpData.ppfjet_twr_hade,
			*inpData.ppfjet_twr_frac,
			*inpData.ppfjet_had_EcalE,
			*inpData.ppfjet_twr_clusterind,
			*inpData.ppfjet_cluster_dR,
			*inpData.ppfjet_had_id,
			*inpData.ppfjet_had_ntwrs,
			*inpData.ppfjet_had_candtrackind,
			*inpData.ppfjet_candtrack_px,
			*inpData.ppfjet_candtrack_py,
			*inpData.ppfjet_candtrack_pz);

    datum->SetGenInfo(inpData.tpfjet_genE,inpData.tpfjet_genpt,
		      inpData.ppfjet_genE,inpData.ppfjet_genpt);

    const int gen_debug=2;
    if (gen_debug==1) {
      std::cout << "iEntry=" << iEntry << ", ppfjet_genE=" << inpData.ppfjet_genE << "\n";
      std::cout << "  ppfjet_twr_ieta " << inpData.ppfjet_twr_ieta << "\n";
      std::cout << "  ppfjet_twr_hade " << inpData.ppfjet_twr_hade << "\n";
      std::cout << "  ppfjet_twr_frac " << inpData.ppfjet_twr_frac << "\n";
      std::cout << "  ppfjet_twr_clusterind " << inpData.ppfjet_twr_clusterind << "\n";
      std::cout << "  ppfjet_cluster_dR " << inpData.ppfjet_cluster_dR << "\n";
    }
    else if (gen_debug==2) {
      int reverse= (inpData.tpfjet_genE < inpData.ppfjet_genE) ? 1:0;
      for (int iLoop=0; iLoop<2; iLoop++) {
	if (((iLoop==0) && !reverse) ||
	    ((iLoop==1) &&  reverse)) {
	  std::cout << "genE=" << inpData.tpfjet_genE << ", pfE=" << inpData.tpfjet_E << ", rhE=" << datum->CalcRecHits_TagJetEn() << "\n";
	}
	else {
	  std::cout << "genE=" << inpData.ppfjet_genE << ", pfE=" << inpData.ppfjet_E << ", rhE=" << datum->CalcRecHits_ProbeJetEn() << "\n";
	}
      }
    }


    if (maxEntries!=nEntries) {
      if (0) {
	inpData.PrintSelectedFields();
	std::cout << (*datum) << "\n";
      }
      datum->PrintBaseClass(1);
    }

    outTree->Fill();
  }

  //outTree->Write();
  outFile.Write();
  delete outTree;
  outFile.Close();
  std::cout << "File <" << outFile.GetName() << "> recreated\n";
  std::cout << "passed " << passedCount << " of " << nEntries << " entries\n";
  return;
}
