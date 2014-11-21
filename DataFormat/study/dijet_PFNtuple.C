#define dijet_PFNtuple_cxx
#include "dijet_PFNtuple.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "helper.hh"

void dijet_PFNtuple::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L dijet_PFNtuple.C
//      Root > dijet_PFNtuple t
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

void dijet_PFNtuple::PrintSelectedFields(int selection) {
  std::cout << "selection=" << selection << "\n";
  std::cout << "dijet (dEta,dPhi)=("
	    << pf_dijet_deta << "," << pf_dijet_dphi << ")\n";
  std::cout << "tag jet  p=" << tpfjet_p << ", E=" << tpfjet_E
	    << ", eta=" << tpfjet_eta << ", phi=" << tpfjet_phi <<"\n";
  std::cout << "probe jet  p=" << ppfjet_p << ", E=" << ppfjet_E
	    << ", eta=" << ppfjet_eta << ", phi=" << ppfjet_phi <<"\n";
  std::cout << "tag unkown_E+electron_E+muon_E+photon_E= "
	    << tpfjet_unkown_E << " + " << tpfjet_electron_E << " + "
	    << tpfjet_muon_E << " + " << tpfjet_photon_E << " = "
	    << (tpfjet_unkown_E + tpfjet_electron_E +
		tpfjet_muon_E + tpfjet_photon_E) << "\n";
  std::cout << "probe unkown_E+electron_E+muon_E+photon_E= "
	    << ppfjet_unkown_E << " + " << ppfjet_electron_E << " + "
	    << ppfjet_muon_E << " + " << ppfjet_photon_E << " = "
	    << (ppfjet_unkown_E + ppfjet_electron_E +
		ppfjet_muon_E + ppfjet_photon_E) << "\n";
#ifdef helper_HH
  printVec("tpfjet_twr_ieta",*tpfjet_twr_ieta);
  printVec("tpfjet_twr_iphi",*tpfjet_twr_iphi);
  printVec("tpfjet_twr_hade",*tpfjet_twr_hade);
  printVec("tpfjet_twr_frac",*tpfjet_twr_frac);
  printVec("tpfjet_had_EcalE",*tpfjet_had_EcalE);

  printVec("ppfjet_twr_ieta",*ppfjet_twr_ieta);
  printVec("ppfjet_twr_iphi",*ppfjet_twr_iphi);
  printVec("ppfjet_twr_hade",*ppfjet_twr_hade);
  printVec("ppfjet_twr_frac",*ppfjet_twr_frac);
  printVec("ppfjet_had_EcalE",*ppfjet_had_EcalE);
#endif
}


// --------------------------------------------------------------
// --------------------------------------------------------------

Double_t dijet_PFNtuple::getSumEcalE(int tagJet, int includeOthers) const {
  double sum=0;

  if (tagJet==1) {
    for (unsigned int i=0; i< tpfjet_had_EcalE->size(); i++) {
      sum+= tpfjet_had_EcalE->at(i);
      if (0 && (tpfjet_had_id->at(i) >= 2)) {
	std::cout << "tpfjet_had_id[" << i << "]=" << tpfjet_had_id->at(i)
		  << ", EcalE=" << tpfjet_had_EcalE->at(i) << "\n";
      }
    }
    if (includeOthers) {
      sum+= tpfjet_unkown_E + tpfjet_electron_E + tpfjet_muon_E +
	tpfjet_photon_E;
    }
  }
  else {
    for (unsigned int i=0; i< ppfjet_had_EcalE->size(); i++) {
      sum+= ppfjet_had_EcalE->at(i);
      if (0 && (ppfjet_had_id->at(i) >= 2)) {
	std::cout << "ppfjet_had_id[" << i << "]=" << ppfjet_had_id->at(i)
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

Double_t dijet_PFNtuple::getSumHcalE_trackDiffEcal(int tagJet) const {
  Double_t sum=0;
  if (tagJet) {
    for (unsigned int i=0; i<tpfjet_had_id->size(); ++i) {
      int trackIdx= tpfjet_had_candtrackind->at(i);
      if ((tpfjet_had_id->at(i) == 0) && // charged hadron
	  (trackIdx>-1)    // has a track
	  && (tpfjet_had_ntwrs->at(i) == 0)  // has no recHits
	  ) {
	sum += sqrt( pow(tpfjet_candtrack_px->at(trackIdx),2) +
		     pow(tpfjet_candtrack_py->at(trackIdx),2) +
		     pow(tpfjet_candtrack_pz->at(trackIdx),2) )
	  - tpfjet_had_EcalE->at(i);
      }
    }
  }
  else {
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
  }
 return sum;
}

// --------------------------------------------------------------

std::map<Int_t,Double_t> dijet_PFNtuple::getHcalEMap
    (int tagJet, double thrContribution) const
{
  std::map<Int_t,Double_t> hcalE;

  if (tagJet) {
    for (unsigned int i=0; i<tpfjet_twr_hade->size(); ++i) {
      if (tpfjet_twr_hade->at(i)<=0) continue;
      int clusterIdx= tpfjet_twr_clusterind->at(i);
      if (0) {
	std::cout << "had_i=" << i << ", clusterIdx=" << clusterIdx;
	if (clusterIdx>=0) {
	  std::cout << ", cluster_dR=" << tpfjet_cluster_dR->at(clusterIdx);
	  if (tpfjet_cluster_dR->at(clusterIdx)>=0.5) {
	    std::cout << ", don't use " << tpfjet_twr_hade->at(i)
		      << " x " << tpfjet_twr_frac->at(i) << "\n";
	  }
	}
	std::cout << "\n";
      }
      if ((clusterIdx<0) || (tpfjet_cluster_dR->at(clusterIdx)<0.5)) {
	int iEta= tpfjet_twr_ieta->at(i);
	const int MAXIETA=41;
	if ((iEta<-MAXIETA) || (iEta>MAXIETA) || (iEta==0)) {
	  std::cout << "getHcalE: bad iEta=" << iEta << "\n";
	}
	double deposit= tpfjet_twr_hade->at(i);
	double fraction= tpfjet_twr_frac->at(i);
	//std::cout << "iEta=" << iEta << ", deposit=" << deposit << ", fraction=" << fraction << "\n";
	if (deposit*fraction < thrContribution) {
	  //std::cout << " .. skip\n";
	  continue;
	}
	hcalE[iEta] += deposit*fraction;
      }
    }
  }
  else {
    for (unsigned int i=0; i<ppfjet_twr_hade->size(); ++i) {
      if (ppfjet_twr_hade->at(i)<=0) continue;
      int clusterIdx= ppfjet_twr_clusterind->at(i);
      if (0) {
	std::cout << "had_i=" << i << ", clusterIdx=" << clusterIdx;
	if (clusterIdx>=0) {
	  std::cout << ", cluster_dR=" << ppfjet_cluster_dR->at(clusterIdx);
	  if (ppfjet_cluster_dR->at(clusterIdx)>=0.5) {
	    std::cout << ", don't use " << ppfjet_twr_hade->at(i)
		      << " x " << ppfjet_twr_frac->at(i) << "\n";
	  }
	}
	std::cout << "\n";
      }
      if ((clusterIdx<0) || (ppfjet_cluster_dR->at(clusterIdx)<0.5)) {
	int iEta= ppfjet_twr_ieta->at(i);
	const int MAXIETA=41;
	if ((iEta<-MAXIETA) || (iEta>MAXIETA) || (iEta==0)) {
	  std::cout << "getHcalE: bad iEta=" << iEta << "\n";
	}
	double deposit= ppfjet_twr_hade->at(i);
	double fraction= ppfjet_twr_frac->at(i);
	//std::cout << "iEta=" << iEta << ", deposit=" << deposit << ", fraction=" << fraction << "\n";
	if (deposit*fraction < thrContribution) {
	  //std::cout << " .. skip\n";
	  continue;
	}
	hcalE[iEta] += deposit*fraction;
      }
    }
  }

  return hcalE;
}

// --------------------------------------------------------------
