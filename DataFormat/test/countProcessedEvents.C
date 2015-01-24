#include "pf_gammajettree.h"
#include "misc_gammajettree.h"

void countProcessedEvents(TString fname, Long64_t maxEntries=-1,
			  int countSelectedEvents=0)
{
  ULong64_t nProcessed=0;
  ULong64_t nChecked=0, nSelected=0;

  std::cout << "\tNote dual use due to countSelectedEvents=0 or 1\n";

  if (!countSelectedEvents) {
    // count processed events
    misc_gammajettree inpMiscTree(fname);

    // read in the file
    Long64_t nEntries= inpMiscTree.fChain->GetEntries();
    std::cout << "nEntries=" << nEntries << "\n";
    if (maxEntries<0) maxEntries=nEntries;
    if (nEntries==0) {
      std::cout << "nEntries=0. Stop counting processed events\n";
      return;
    }
    Long64_t nBytes=0; //, passedCount=0;

    for (Long64_t iEntry=0; (iEntry<nEntries) && (iEntry<maxEntries);
	 iEntry++) {
      if (inpMiscTree.LoadTree(iEntry) < 0) break;
      std::cout << " ... reading entry " << iEntry << "\n";
      Long64_t nb = inpMiscTree.GetEntry(iEntry);
      nBytes += nb;

      nProcessed+= inpMiscTree.nProcessed;
      std::cout << "... nProcessed=" << inpMiscTree.nProcessed
		<< ", total=" << nProcessed << "\n";
    }

  }

  // Count selected events
  if (countSelectedEvents==1) {
    pf_gammajettree inpData(fname);

    inpData.DeactivateBranches();

    // photon
    inpData.ActivateBranches
      (4,
       "tagPho_pt","tagPho_phi",
       "tagPho_idLoose","tagPho_idTight");
    // jet id
    inpData.ActivateBranches
      (9,
       "ppfjet_pt","ppfjet_eta",
       "ppfjet_phi",
       "ppfjet_NeutralHadronFrac","ppfjet_NeutralEMFrac",
       "ppfjet_nConstituents", "ppfjet_ChargedHadronFrac",
       "ppfjet_ChargedMultiplicity","ppfjet_ChargedEMFrac");
    // alpha
    inpData.ActivateBranches
      (1, "pfjet2_pt");

    // read in the file
    Long64_t nEntries= inpData.fChain->GetEntries();
    std::cout << "nEntries=" << nEntries << "\n";
    if (maxEntries<0) maxEntries=nEntries;
    Long64_t nBytes=0;

    for (Long64_t iEntry=0; (iEntry<nEntries) && (iEntry<maxEntries);
	 iEntry++) {
      if (inpData.LoadTree(iEntry) < 0) {
	std::cout << "failed to get the entry iEntry=" << iEntry << "\n";
	break;
      }
      Long64_t nb = inpData.GetEntry(iEntry);
      nBytes += nb;
      nChecked++;
      if (iEntry%10000==0) {
	std::cout << " ... reading entry " << iEntry
		  << " (selected so far=" << nSelected << ")\n";
      }
      //std::cout << "ientry=" << iEntry << "\n";

      if (!inpData.passCuts(20.,int(_phoLooseID),1)) continue; // recommended use
      //if ( inpData.pfjet2_pt/inpData.tagPho_pt > 0.05) continue; // stricter alpha cut

      nSelected++;
    }
  }

  std::cout << "\nFile name=<" << fname << ">\n";
  if (!countSelectedEvents) {
    std::cout << " total number of processed events=" << nProcessed << " (into ntuple)\n";
  }
  if (countSelectedEvents) {
    std::cout << " total number of   checked events=" << nChecked << "\n";
    std::cout << " total number of  selected events=" << nSelected << " ("
	      << (nSelected*100/double(nChecked)) << "%)\n";
  }

  return;
}
