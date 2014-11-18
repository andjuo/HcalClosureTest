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
		    int nBins=100, double EtMin=0., double EtMax=1000) const;

  TH1D* plot_Energy(const char *hName, const char *hTitle,
		    int tag=0, int transverse=0,
		    const TArrayD *hcalCorrCf=NULL,
		    int nBins=100, double EtMin=0., double EtMax=1000) const;

  TH1D* plot_EnergyOverGenEnergy
                   (const char *hName, const char *hTitle,
		    int tag=0, int transverse=-1, // transverse is ignored
		    const TArrayD *hcalCorrCf=NULL,
		    int nBins=100, double RatioMin=0.,double RatioMax=2.) const;

  TH2D* plot_TowerEn(const char *hNameBase, const char *hTitle,
		     unsigned int idxMin, unsigned int idxMax,
		     int plotProbe=1,
		     const TArrayD *hcalCorrCf=NULL) const;

  TH2D* plot_TowerFitProfile(const char *hName, const char *hTitle,
			     int normalized,
			     int nBins=50,
			     double cfMin=0., double cfMax=2.,
			     const std::vector<double> *setCfs=NULL) const;

  TH1D* plot_TowerEntryCount(const char *hName, const char *hTitle,
			     int weighted,
			     double minWeight=0., int fractionFromMax=1,
			     std::vector<Double_t> *towers_ptr=NULL) const;

};

// --------------------------------------------------------------
// --------------------------------------------------------------

void plotTwoHistos(TCanvas *cx, TString title,
		   TH1D* h1, TString label1, int color1, TString drawOpt1,
		   TH1D* h2, TString label2, int color2, TString drawOpt2);

// --------------------------------------------------------------

#endif
