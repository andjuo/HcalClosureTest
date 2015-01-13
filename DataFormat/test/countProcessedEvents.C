#include "misc_gammajettree.h"

void countProcessedEvents(TString fname, Long64_t maxEntries=-1)
{
  misc_gammajettree inpMiscTree(fname);

    // read in the file
  Long64_t nEntries= inpMiscTree.fChain->GetEntries();
  if (maxEntries<0) maxEntries=nEntries;
  if (nEntries==0) {
    std::cout << "nEntries=0. Stop\n";
    return;
  }
  Long64_t nBytes=0; //, passedCount=0;
  ULong64_t nProcessed=0;

  for (Long64_t iEntry=0; (iEntry<nEntries) && (iEntry<maxEntries);
       iEntry++) {
    if (inpMiscTree.LoadTree(iEntry) < 0) break;
    std::cout << " ... reading entry " << iEntry << "\n";
    Long64_t nb = inpMiscTree.GetEntry(iEntry);
    nBytes += nb;

    nProcessed+= inpMiscTree.nProcessed;
    std::cout << "nProcessed=" << inpMiscTree.nProcessed
	      << ", total=" << nProcessed << "\n";
  }
  return;
}
