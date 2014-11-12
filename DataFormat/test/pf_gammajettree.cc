#define pf_gammajettree_cxx
#include "pf_gammajettree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void pf_gammajettree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L pf_gammajettree.C
//      Root > pf_gammajettree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}

// --------------------------------------------------------------

void pf_gammajettree::PrintSelectedFields(int selection) {
  std::cout << "PrintSelectedFields(selection=" << selection << ")\n";
  std::cout << " ... not ready\n";
}


// --------------------------------------------------------------

void pf_gammajettree::ActivateBranches_forRecHitsEnergyCalc() {
  ActivateBranches(3, "tagPho_energy","tagPho_eta","tagPho_phi");
  ActivateBranches(2, "ppfjet_eta","ppfjet_phi");
  ActivateBranches(4, "ppfjet_unkown_E","ppfjet_electron_E",
			   "ppfjet_muon_E","ppfjet_photon_E");
  ActivateBranches(4,"ppfjet_had_EcalE","ppfjet_had_id",
		   "ppfjet_had_ntwrs", "ppfjet_had_candtrackind");
  ActivateBranches(4, "ppfjet_twr_ieta","ppfjet_twr_hade",
			 "ppfjet_twr_frac","ppfjet_twr_clusterind");
  ActivateBranches(1,"ppfjet_cluster_dR");
  ActivateBranches(3,"ppfjet_candtrack_px","ppfjet_candtrack_py",
			 "ppfjet_candtrack_pz");
}


// --------------------------------------------------------------

void pf_gammajettree::ActivateBranches_forFitSkim() {
  ActivateBranches(3, "EventNumber","RunNumber","EventWeight");
  ActivateBranches(1, "tagPho_pt");
  ActivateBranches(2, "tagPho_idTight","tagPho_idLoose");
  ActivateBranches(2, "nPhotons","nPFJets");
  ActivateBranches(1, "pfjet2_pt");
  ActivateBranches_forRecHitsEnergyCalc();
}

// --------------------------------------------------------------

void pf_gammajettree::ActivateBranches_genBasicSet() {
  ActivateBranches(5,"tagPho_genPt","tagPho_genEnergy","tagPho_genEta",
		   "tagPho_genPhi","tagPho_genDeltaR");
  ActivateBranches(4,"ppfjet_genpt","ppfjet_genp","ppfjet_genE",
		   "ppfjet_gendr");
}

// --------------------------------------------------------------

Double_t pf_gammajettree::getSumEcalE(int tag, int includeOthers) const {
  double sum=0;
  if (tag==1) {
    sum= tagPho_energy;
  }
  else {
    for (unsigned int i=0; i< ppfjet_had_EcalE->size(); i++) {
      if (ppfjet_had_id->at(i) < 2) { // hadron not in HF
	sum+= ppfjet_had_EcalE->at(i);
      }
      else {
	std::cout << "ppfjet_had_id[" <<i << "]=" << ppfjet_had_id->at(i)
		  << ", EcalE=" << ppfjet_had_EcalE->at(i) << "\n";
      }
    }
    if (includeOthers) {
      sum+= ppfjet_unkown_E + ppfjet_electron_E + ppfjet_muon_E +
	ppfjet_photon_E;
    }
  }
  return sum;
}

// --------------------------------------------------------------

Double_t pf_gammajettree::getSumHcalE_trackDiffEcal(int leadingJet) const {
  if (!leadingJet) {
    std::cout << "pf_gammajettree::getSumHcalE_trackOnly: "
	      << "(!leadingJet) is not ready\n";
    return 0;
  }

  Double_t sum=0;
  for (unsigned int i=0; i<ppfjet_had_id->size(); ++i) {
    int trackIdx= ppfjet_had_candtrackind->at(i);
    if ((ppfjet_had_id->at(i) == 0) && // charged hadron
	(trackIdx>-1)    // has a track
	&& (ppfjet_had_ntwrs->at(i) == 0)  // has no recHits
	) {
      sum += sqrt( pow(ppfjet_candtrack_px->at(trackIdx),2) +
		   pow(ppfjet_candtrack_py->at(trackIdx),2) +
		   pow(ppfjet_candtrack_pz->at(trackIdx),2) )
	- ppfjet_had_EcalE->at(i);
    }
  }
  return sum;
}

// --------------------------------------------------------------

std::map<Int_t,Double_t> pf_gammajettree::getHcalEMap
    (int leadingJet, double thrContribution) const {
  std::map<Int_t,Double_t> hcalE;
  if (!leadingJet) {
    std::cout << "pf_gammajettree::getHcalE: (!leadingJet) is not ready\n";
    return hcalE;
  }

  for (unsigned int i=0; i<ppfjet_twr_hade->size(); ++i) {
    if (ppfjet_twr_hade->at(i)<=0) continue;
    int clusterIdx= ppfjet_twr_clusterind->at(i);
    if ((clusterIdx<0) || (ppfjet_cluster_dR->at(clusterIdx)<0.5)) {
      int iEta= ppfjet_twr_ieta->at(i);
      const int MAXIETA=41;
      if ((iEta<-MAXIETA) || (iEta>MAXIETA) || (iEta==0)) {
	std::cout << "getHcalE: bad iEta=" << iEta << "\n";
      }
      double deposit= ppfjet_twr_hade->at(i);
      double fraction= ppfjet_twr_frac->at(i);
      std::cout << "iEta=" << iEta << ", deposit=" << deposit << ", fraction=" << fraction << "\n";
      if (deposit*fraction < thrContribution) {
	std::cout << " .. skip\n";
	continue;
      }
      hcalE[iEta] += deposit*fraction;
    }
  }

  return hcalE;
}

// --------------------------------------------------------------
