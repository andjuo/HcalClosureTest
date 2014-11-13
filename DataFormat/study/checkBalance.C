#include "dijet_PFNtuple.h"
#include <TH1D.h>
#include <TCanvas.h>
#include <TLeaf.h>

void checkBalance(int theCase=1, int imax=1, int targetEntry=-1) {
  dijet_PFNtuple inpData;

  inpData.DeactivateBranches();
  //inpData.ActivateBranches(4,"tpfjet_pt","ppfjet_pt","pf_thirdjet_px","pf_thirdjet_py");

  std::vector<TString> brTargetVec;
  std::vector<TString> brAddVec;
  brAddVec.reserve(10);

  switch(theCase) {
  case 1:
    brTargetVec.push_back("tpfjet_E");
    brAddVec.push_back("tpfjet_unkown_E");
    brAddVec.push_back("tpfjet_electron_E");
    brAddVec.push_back("tpfjet_muon_E");
    brAddVec.push_back("tpfjet_photon_E");
    brAddVec.push_back("tpfjet_had_n");
    brAddVec.push_back("tpfjet_had_E");
    break;
  case 2:
  case 3:
    brTargetVec.push_back("tpfjet_E");
    brTargetVec.push_back("tpfjet_phi");
  default:
    std::cout << "not ready for theCase=" << theCase << "\n";
    return;
  }

  inpData.ActivateBranches(brTargetVec);
  inpData.ActivateBranches(brAddVec);

  Long64_t nentries = inpData.fChain->GetEntries();
  if (nentries<=0) return;
  if (imax<0) imax=int(nentries);

  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; (jentry<nentries) && (jentry<imax); jentry++) {
    Long64_t ientry = inpData.LoadTree(jentry);
    if (ientry < 0) break;
    nb = inpData.GetEntry(jentry);   nbytes += nb;
    if (jentry%100000==0) std::cout << "... reading entry " << jentry << "\n";
    // if (Cut(ientry) < 0) continue;

    double targetVal=0;
    TLeaf *targetLeaf= inpData.fChain->GetLeaf(brTargetVec[0]);
    targetLeaf->GetBranch()->GetEntry(jentry);
    targetVal= targetLeaf->GetValue();

    if (0 && (theCase==1)) {
      std::cout << " code check. Compare: " << targetLeaf->GetValue() << " vs " << inpData.tpfjet_E << "\n";
    }

    double addedVal=0.;
    for (unsigned int i=0; i<brAddVec.size(); ++i) {
      TString brName= brAddVec[i];
      TLeaf *leaf= inpData.fChain->GetLeaf(brName);
      leaf->GetBranch()->GetEntry(jentry);
      std::cout << "branch=" << brName << ": ";
      if (brName.Index("_had_")==-1) {
	double val= leaf->GetValue();
	std::cout << "  add: " << val << "\n";
	addedVal+=val;
      }
      else {
	/*
	  if (!leaf->GetValuePointer()) {
	  std::cout << "null ptr\n";
	  return;
	}
	*/
	if (brName.Index("_had_n")!=-1) {
	  std::cout << "value=" << leaf->GetValue() << " ";
	  std::cout << "skipped\n";
	}
	else {
	  std::vector<float> *fVec= NULL;
	  if (brName="tpfjet_had_E") fVec=inpData.tpfjet_had_E;
	  else if (brName="tpfjet_had_px") fVec=inpData.tpfjet_had_px;
	  else if (brName="tpfjet_had_py") fVec=inpData.tpfjet_had_py;
	  if (!fVec) { std::cout << "null *_had_* branch\n"; return; }
	  std::cout << " size=" << fVec->size() << +"\n";
	  for (unsigned int iv=0; iv<fVec->size(); ++iv) {
	    double val= fVec->at(iv);
	    std::cout << " .. val=" << val << "\n";
	    addedVal+= val;
	  }
	}
      }
    }
    std::cout << "targetVal=" << targetVal << ", addedVal=" << addedVal << "\n";

    //double jet3pt= sqrt(pow(inpData.pf_thirdjet_px,2) + pow(inpData.pf_thirdjet_py,2));
    //double jet12ptAve= 0.5*(inpData.tpfjet_pt + inpData.ppfjet_pt);
    //double ratio=(jet12ptAve==0) ? -9.99 : jet3pt/jet12ptAve;
    //if (ratio<0) badRatioCount++;

  }
}
