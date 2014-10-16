#ifndef Analyzer_hh
#define Analyzer_hh

#include <TROOT.h>
#include "DijetRespCorrDataExtended.h"
#include "CPlot.hh"

// -------------------------------------------------------

class Cuts_t {
protected:
  TString fName;
  double fJetEtDiff_min, fJetEtDiff_max;
  double fJetRecHitsEtDiff_min, fJetRecHitsEtDiff_max;
  double fDijet_dPhi_min, fDijet_dPhi_max;
  double fDijetEtDiffRatio_min, fDijetEtDiffRatio_max; // sumHad over nonHad
public:
  Cuts_t(const TString &new_name="cuts");
  Cuts_t(const TString &new_name, const Cuts_t &c);

  TString name() const { return fName; }

  void swapEtFields(const TString &search_for,
		    const TString &replace_with);

  void setEtDiff(double new_etDiff_min, double new_etDiff_max,
		 double new_RHEtDiff_min, double new_RHEtDiff_max) {
    fJetEtDiff_min= new_etDiff_min;
    fJetEtDiff_max= new_etDiff_max;
    fJetRecHitsEtDiff_min=new_RHEtDiff_min;
    fJetRecHitsEtDiff_max=new_RHEtDiff_max;
  }

  void setEtDiffRatio(double sumHad_over_nonHad_min,
		      double sumHad_over_nonHad_max) {
    fDijetEtDiffRatio_min = sumHad_over_nonHad_min;
    fDijetEtDiffRatio_max = sumHad_over_nonHad_max;
  }

  void setDijetDiPhi(double new_dPhi_min, double new_dPhi_max) {
    fDijet_dPhi_min= new_dPhi_min;
    fDijet_dPhi_max= new_dPhi_max;
  }

  int passCuts(const DijetRespCorrDatumExtended_t &d) const;
};

// -------------------------------------------------------
/*
class Analyzer_t {
protected:
  TString fName;
  Cuts_t fCuts;
  TH1D* fhTow;
  
public:
  Analyzer_t(TString name, const Cuts_t &cut);
protected:
  int create();

public:
  int Fill(const DijetRespCorrDatumExtended_t &d);
  
};
*/
// -------------------------------------------------------
// -------------------------------------------------------

//typedef enum { _fillTower_count, _fillTower_enFrac } TFillTower_t;

void FillTowerCount(TH1D* histo, const DijetRespCorrDatum &d, int jetId=-1,
		    int useWeight=1);

// -------------------------------------------------------

#endif
