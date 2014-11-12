#ifndef GammaJetFitAnalyzer_H
#define GammaJetFitAnalyzer_H

#ifdef __localRun
#  include "../interface/GammaJetFitData.h"
#else
#  include "HcalClosureTest/DataFormat/interface/GammaJetFitData.h"
#endif

#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>

// --------------------------------------------------------------

class GammaJetFitAnalyzer_t {
  const GammaJetFitter_t* fData;
 public:
  GammaJetFitAnalyzer_t(const GammaJetFitter_t *d)
    : fData(d)
    {}

 public:
  TH2D* plot_EtVsEt(const char *hName, const char *hTitle,
		    const TArrayD *hcalCorrCf=NULL,
		    int nBins=100,
		    double EtMin=0., double EtMax=1000) const;

  TH2D* plot_TowerEn(const char *hNameBase, const char *hTitle,
		     unsigned int idxMin, unsigned int idxMax,
		     int plotProbe=1,
		     const TArrayD *hcalCorrCf=NULL) const;

  TH2D* plot_TowerFitProfile(const char *hName, const char *hTitle,
			     int normalized,
			     int nBins=50,
			     double cfMin=0., double cfMax=2.,
			     const std::vector<double> *setCfs=NULL) const;
};

// --------------------------------------------------------------

#endif
