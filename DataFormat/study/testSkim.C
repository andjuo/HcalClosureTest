#include "DijetRespCorrDataExtended.h"
#include "helper.hh"
#include <TFile.h>
#include <TTree.h>
#include <assert.h>

void testSkim(Long64_t maxEntries=10) {
  DijetRespCorrDatumExtended_t *d= new DijetRespCorrDatumExtended_t();
  TFile inpFile("skim.root","READ");
  TTree *inpTree = (TTree*)inpFile.Get("dijet_data");
  assert(inpTree);
  inpTree->SetBranchAddress("dijet_data",&d);
  
  Long64_t nEntries= inpTree->GetEntriesFast();
  if (maxEntries<0) maxEntries=nEntries;

  for (Long64_t iEntry=0; (iEntry<nEntries) && (iEntry<maxEntries);
       ++iEntry) {
    inpTree->GetEntry(iEntry);
    std::cout << "iEntry=" << iEntry << "\n";
    if (0) {
      std::cout << (*d) << "\n";
    }
    else d->PrintBaseClass(1);
  }
  inpFile.Close();
}
