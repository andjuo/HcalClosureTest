#ifndef link_GammaJetFit_H
#define link_GammaJetFit_H

#define __localRun
#include "../interface/GammaJetFitData.h"

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TString.h>

#include "helper.h"

// -----------------------------------------------------------

int LoadGammaJetEvents(const TString fname,
		       const GammaJetCuts_t &cuts,
		       GammaJetFitter_t &fitter,
		       Long64_t maxEntries=-1);

// -----------------------------------------------------------

int GetEmptyTowers(const GammaJetFitter_t &fitter,
		   std::vector<Int_t> &towers,
		   double minWeight=0.);

// -----------------------------------------------------------

#endif
