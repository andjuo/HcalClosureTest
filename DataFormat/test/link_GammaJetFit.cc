#include "link_GammaJetFit.h"
#include "../src/GammaJetFitData.cc"

void link_GammaJetFit() {
  std::cout << "link_GammaJetFit\n";
}


// -----------------------------------------------------------

int LoadGammaJetEvents(const TString fname,
		       const GammaJetCuts_t &cuts,
		       GammaJetFitter_t &fitter,
		       Long64_t maxEntries)
{
  GammaJetEvent_t *dt= new GammaJetEvent_t();
  TFile *inpF= new TFile(fname,"read");
  if (!inpF || !inpF->IsOpen()) {
    std::cout << "LoadGammaJetEvents: failed to open the file <"
	      << fname << ">\n";
    return 0;
  }
  TTree *inpTree= (TTree*)inpF->Get("gjet_data");
  if (!inpTree) {
    std::cout << "null inpTree\n";
    return 0;
  }
  inpTree->SetBranchAddress("gjet_data",&dt);

  Long64_t nEntries= inpTree->GetEntriesFast();
  if (maxEntries<0) maxEntries=nEntries;

  for (Long64_t iEntry=0; (iEntry<nEntries) && (iEntry<maxEntries);
       ++iEntry) {
    inpTree->GetEntry(iEntry);
    if (cuts.passCuts(*dt)) {
      fitter.push_back(*dt);
    }
  }
  inpF->Close();
  HERE("file closed");
  delete inpF;
  HERE("file deleted");
  //delete inpTree;
  //HERE("tree deleted");
  return 1;
}

// -----------------------------------------------------------

int GetEmptyTowers(const GammaJetFitter_t &fitter,
		   std::vector<Int_t> &towers,
		   double minWeight) {
  towers.clear();
  towers.reserve(NUMTOWERS);
  std::vector<Double_t> countV(NUMTOWERS,0);

  const std::vector<GammaJetEvent_t*> *d= & fitter.GetData();
  for (unsigned int i=0; i<d->size(); ++i) {
    const GammaJetEvent_t *e= d->at(i);
    double w= e->GetWeight();
    for (int iProbe=0; iProbe<2; ++iProbe) {
      const std::map<Int_t,Double_t>* hMap= 
	(iProbe==1) ? &e->GetProbeHcalE() : &e->GetTagHcalE();
      if (hMap->size()==0) continue;
      for (std::map<Int_t,Double_t>::const_iterator it= hMap->begin();
	   it!=hMap->end(); it++) {
	countV[ it->first + MAXIETA ] += w;
      }
    }
  }

  for (unsigned int i=0; i<countV.size(); ++i) {
    if (countV[i]<minWeight) {
      towers.push_back(i-MAXIETA);
    }
  }
  return 1;
}

// -----------------------------------------------------------
