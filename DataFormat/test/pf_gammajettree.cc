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
  ActivateBranches(3, "ppfjet_eta","ppfjet_phi","ppfjet_pt");
  ActivateBranches(4, "ppfjet_unkown_E","ppfjet_electron_E",
			   "ppfjet_muon_E","ppfjet_photon_E");
  ActivateBranches(4,"ppfjet_had_EcalE","ppfjet_had_id",
		   "ppfjet_had_ntwrs", "ppfjet_had_candtrackind");
  ActivateBranches(2, "ppfjet_twr_hadind", "ppfjet_had_n");
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
  ActivateBranches(3, "tagPho_idTight","tagPho_idLoose","tagPho_pixelSeed");
  ActivateBranches(2, "nPhotons","nPFJets");
  ActivateBranches(1, "pfjet2_pt");
  ActivateBranches(2, "tagPho_genEnergy", "ppfjet_genE");
  ActivateBranches_forRecHitsEnergyCalc();
  ActivateBranches_jetID(1);
}

// --------------------------------------------------------------

void pf_gammajettree::ActivateBranches_genBasicSet() {
  ActivateBranches(5,"tagPho_genPt","tagPho_genEnergy","tagPho_genEta",
		   "tagPho_genPhi","tagPho_genDeltaR");
  ActivateBranches(4,"ppfjet_genpt","ppfjet_genp","ppfjet_genE",
		   "ppfjet_gendr");
}

// --------------------------------------------------------------

void pf_gammajettree::ActivateBranches_jetID(int leadingJet) {
  TString base=(leadingJet) ? "ppfjet_" : "pfjet2_";
  std::vector<TString> brV;
  brV.reserve(6);
  brV.push_back(base + TString("NeutralHadronFrac"));
  brV.push_back(base + TString("NeutralEMFrac"));
  brV.push_back(base + TString("nConstituents"));
  brV.push_back(base + TString("ChargedHadronFrac"));
  brV.push_back(base + TString("ChargedMultiplicity"));
  brV.push_back(base + TString("ChargedEMFrac"));
  ActivateBranches(brV);
}

// --------------------------------------------------------------

TString explain_failingTightJetIDFlag(int flag, int leadingJet,
				      const pf_gammajettree &e) {
  TString ans=Form("flag=%d,leadingJet=%d ",flag,leadingJet);
  if ((flag&1)!=0) ans.Append(" NeutralHadronFrac>=0.90 ");
  if ((flag&2)!=0) ans.Append(" NeutralEMFrac>=0.90 ");
  if ((flag&4)!=0) ans.Append(" nConstituents<=1 ");
  if ((flag&8)!=0) ans.Append(" (Barrel)ChargedHadronFrac=0" );
  if ((flag&16)!=0) ans.Append(" (Barrel)ChargedMultiplicity=0 ");
  if ((flag&32)!=0) ans.Append(" (Barrel)ChargedEMFrac>=0.99 ");
  if (0 && leadingJet) std::cout << e.ppfjet_nConstituents << "\n";
  return ans;
}

// --------------------------------------------------------------

int pf_gammajettree::passTightJetID(int leadingJet) const {
  int flag=0;
  if (leadingJet) {
    if (ppfjet_NeutralHadronFrac >= 0.90) flag|=1;
    if (ppfjet_NeutralEMFrac >= 0.90) flag|=2;
    if (ppfjet_nConstituents <= 1) flag|=4;
    if (fabs(ppfjet_eta) < 2.4) {
      // additional requirements for the barrel jets
      if (ppfjet_ChargedHadronFrac == Float_t(0)) flag|=8;
      if (ppfjet_ChargedMultiplicity == Float_t(0)) flag|=16;
      if (ppfjet_ChargedEMFrac >= 0.99) flag|=32;
    }
  }
  else {
    if (pfjet2_NeutralHadronFrac >= 0.90) flag|=1;
    if (pfjet2_NeutralEMFrac >= 0.90) flag|=2;
    if (pfjet2_nConstituents <= 1) flag|=4;
    if (fabs(pfjet2_eta) < 2.4) {
      // additional requirements for the barrel jets
      if (pfjet2_ChargedHadronFrac == Float_t(0)) flag|=8;
      if (pfjet2_ChargedMultiplicity == Float_t(0)) flag|=16;
      if (pfjet2_ChargedEMFrac >= 0.99) flag|=32;
    }
  }
  if (0 && (flag!=0)) std::cout << "passTightJetID failed: "
	     << explain_failingTightJetIDFlag(flag,leadingJet,*this) << "\n";
  return (flag==0) ? 1:0;
}

// --------------------------------------------------------------

TString explain_failingPhotonJetFlag(int flag, const pf_gammajettree &e) {
  TString ans = Form("flag=%d: ",flag);
  if ((flag&1)!=0) ans.Append(" tagPho_pt<=10 ");
  if ((flag&2)!=0) ans.Append(" ppfjet_pt<=10 ");
  if ((flag&4)!=0) {
    double pfjet2_pt= e.pfjet2_pt;
    double tagPho_pt= e.tagPho_pt;
    ans.Append(Form(" pfjet2_pt/tagPho_pt=%6.4lf/%6.4lf=%6.4lf ",
		    pfjet2_pt,tagPho_pt,pfjet2_pt/tagPho_pt));
  }
  if ((flag&8)!=0) {
    double tagPho_phi= e.tagPho_phi;
    double ppfjet_phi= e.ppfjet_phi;
    ans.Append(Form("|dPhi|=|%4.2lf-%4.2lf|=%4.2lf<=2.95",
		    tagPho_phi,ppfjet_phi, fabs(tagPho_phi-ppfjet_phi)));
  }
  if ((flag&16)!=0) ans.Append(" !tagPho_idLoose ");
  if ((flag&32)!=0) ans.Append(" !tagPho_idTight ");
  if ((flag&64)!=0) ans.Append(" tagPho_pixelSeed!=0 ");
  return ans;
}

// --------------------------------------------------------------

int pf_gammajettree::passPhotonJetRequirements(int theSet) const {
  if (theSet==0) return 1;
  int flag=0;
  if (tagPho_pt <= 10.) flag|=1;
  if (ppfjet_pt <= 10.) flag|=2;
  if ((tagPho_pt>0.) && (pfjet2_pt/tagPho_pt>=0.2)) flag|=4;
  if (fabs(tagPho_phi - ppfjet_phi) <= 2.95) flag|=8;
  if (theSet==1) {
    if (tagPho_idLoose==0) flag|=16;
  }
  else if (theSet==2) {
    if (tagPho_idTight==0) flag|=32;
  }
  //if (tagPho_pixelSeed != 0) flag|=64;
  if (0 && (flag!=0)) std::cout << "passPhotonJetRequirements failed: "
			<< explain_failingPhotonJetFlag(flag,*this) << "\n";
  return (flag==0) ? 1:0;
}

// --------------------------------------------------------------

Double_t pf_gammajettree::getSumEcalE(int tag, int includeOthers) const {
  double sum=0;
  if (tag==1) {
    sum= tagPho_energy;
  }
  else if (tag==0) {
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
  else if (tag==2) {
    for (unsigned int i=0; i< pfjet2_had_EcalE->size(); i++) {
      sum+= pfjet2_had_EcalE->at(i);
      if (0 && (pfjet2_had_id->at(i) >= 2)) { // hadron in HF
	std::cout << "pfjet2_had_id[" << i << "]=" << pfjet2_had_id->at(i)
		  << ", EcalE=" << pfjet2_had_EcalE->at(i) << "\n";
      }
    }
    if (includeOthers) {
      sum+= pfjet2_unkown_E + pfjet2_electron_E + pfjet2_muon_E +
	pfjet2_photon_E;
    }
  }
  else {
    std::cout << "getSumEcalE: wrong value of tag=" << tag << "\n";
  }
  return sum;
}

// --------------------------------------------------------------

Double_t pf_gammajettree::getSumHcalE_trackDiffEcal(int leadingJet) const {
  Double_t sum=0;
  if (leadingJet) {
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
  else {
    for (unsigned int i=0; i<pfjet2_had_id->size(); ++i) {
      int trackIdx= pfjet2_had_candtrackind->at(i);
      if ((pfjet2_had_id->at(i) == 0) && // charged hadron
	  (trackIdx>-1)    // has a track
	  && (pfjet2_had_ntwrs->at(i) == 0)  // has no recHits
	  ) {
	sum += sqrt( pow(pfjet2_candtrack_px->at(trackIdx),2) +
		     pow(pfjet2_candtrack_py->at(trackIdx),2) +
		     pow(pfjet2_candtrack_pz->at(trackIdx),2) )
	  - pfjet2_had_EcalE->at(i);
      }
    }
  }
  return sum;
}

// --------------------------------------------------------------

std::map<Int_t,Double_t> pf_gammajettree::getHcalEMap
    (int leadingJet, double thrContribution) const {
  std::map<Int_t,Double_t> hcalE;

  if (leadingJet) {
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
  else {
    for (unsigned int i=0; i<pfjet2_twr_hade->size(); ++i) {
      if (pfjet2_twr_hade->at(i)<=0) continue;
      int clusterIdx= pfjet2_twr_clusterind->at(i);
      if ((clusterIdx<0) || (pfjet2_cluster_dR->at(clusterIdx)<0.5)) {
	int iEta= pfjet2_twr_ieta->at(i);
	const int MAXIETA=41;
	if ((iEta<-MAXIETA) || (iEta>MAXIETA) || (iEta==0)) {
	  std::cout << "getHcalE: bad iEta=" << iEta << "\n";
	}
	double deposit= pfjet2_twr_hade->at(i);
	double fraction= pfjet2_twr_frac->at(i);
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
