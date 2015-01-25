#include "pf_gammajettree.h"
#include "misc_gammajettree.h"
#include "helper.h"

void countProcessedEvents(TString fname, Long64_t maxEntries=-1,
			  int countSelectedEvents=0, int listTriggers=1,
			  TString restrictTriggers="|")
{
  ULong64_t nProcessed=0;
  ULong64_t nChecked=0, nSelected=0;

  std::cout << "\tNote dual use due to countSelectedEvents=0 or 1\n";

  std::vector<int> chkPhoTriggers,chkJetTriggers;
  if (!identifyTriggerIndices(restrictTriggers,
			      chkPhoTriggers,chkJetTriggers,0)) return;



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
      if (listTriggers) {
	std::cout << " ignoreHLT=" << inpMiscTree.ignoreHLT << "\n";
	std::cout << " doPFJets=" << inpMiscTree.doPFJets << "\n";
	std::cout << " doGenJets=" << inpMiscTree.doGenJets << "\n";
	std::cout << " workOnAOD=" << inpMiscTree.workOnAOD << "\n";
	if (inpMiscTree.photonTriggerNames) {
	  std::cout << " " << inpMiscTree.photonTriggerNames->size()
		    << " photon triggers:\n";
	  for (unsigned int i=0; i<inpMiscTree.photonTriggerNames->size();i++){
	    std::cout << " " << i << " "
		      << inpMiscTree.photonTriggerNames->at(i);
	    if (std::find(chkPhoTriggers.begin(),chkPhoTriggers.end(), int(i))
		!= chkPhoTriggers.end()) std::cout << " -- marked";
	    std::cout << "\n";
	  }
	}
	else {
	  std::cout << " photonTriggerNames is null\n";
	}
	if (inpMiscTree.jetTriggerNames) {
	  std::cout << " " << inpMiscTree.jetTriggerNames->size()
		    << " jet triggers :\n";
	  for (unsigned int i=0; i<inpMiscTree.jetTriggerNames->size();i++){
	    std::cout << " " << i << " "
		      << inpMiscTree.jetTriggerNames->at(i);
	    if (std::find(chkJetTriggers.begin(),chkJetTriggers.end(), int(i))
		!= chkJetTriggers.end()) std::cout << " -- marked";
	    std::cout << "\n";
	  }
	}
	else {
	  std::cout << " jetTriggerNames is null\n";
	}
      }
    }

  }

  // Count selected events
  if (countSelectedEvents>=1) {
    pf_gammajettree inpData(fname);

    inpData.DeactivateBranches();

    // trigger
    inpData.ActivateBranches (2, "photonTrig_fired","jetTrig_fired");
    // photon
    inpData.ActivateBranches
      (5,
       "tagPho_pt","tagPho_phi","tagPho_eta",
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

      if (!inpData.passCuts(20.,int(_phoTightID),1)) continue; // recommended use
      if ((countSelectedEvents>1) &&
	  ((fabs(inpData.tagPho_eta)>2.4) || (fabs(inpData.ppfjet_eta)>2.4))) {
	continue;
      }

      if ((restrictTriggers.Length()>1) &&
	  !inpData.hltTriggerFired(chkPhoTriggers,chkJetTriggers)) continue;

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
