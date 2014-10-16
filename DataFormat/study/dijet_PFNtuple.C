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
