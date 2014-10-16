#include "Analyzer.hh"

// -------------------------------------------------------

template <class T>
inline void SWAP(T &a, T &b) {
  T tmp=a; a=b; b=tmp;
}

// -------------------------------------------------------
// -------------------------------------------------------

Cuts_t::Cuts_t(const TString &new_name) :
  fName(new_name),
  fJetEtDiff_min(-1), fJetEtDiff_max(-1),
  fJetRecHitsEtDiff_min(-1), fJetRecHitsEtDiff_max(-1),
  fDijet_dPhi_min(-1), fDijet_dPhi_max(-1),
  fDijetEtDiffRatio_min(-1e9), fDijetEtDiffRatio_max(1e9)
{}

// -------------------------------------------------------

Cuts_t::Cuts_t(const TString &new_name, const Cuts_t &c) :
  fName(new_name),
  fJetEtDiff_min(c.fJetEtDiff_min),
  fJetEtDiff_max(c.fJetEtDiff_max),
  fJetRecHitsEtDiff_min(c.fJetRecHitsEtDiff_min),
  fJetRecHitsEtDiff_max(c.fJetRecHitsEtDiff_max),
  fDijet_dPhi_min(c.fDijet_dPhi_min),
  fDijet_dPhi_max(c.fDijet_dPhi_max),
  fDijetEtDiffRatio_min(c.fDijetEtDiffRatio_min),
  fDijetEtDiffRatio_max(c.fDijetEtDiffRatio_max)
{}

// -------------------------------------------------------

void Cuts_t::swapEtFields(const TString &search_for,
			  const TString &replace_with) {
  fName.ReplaceAll(search_for,replace_with);
  SWAP(fJetEtDiff_min, fJetRecHitsEtDiff_min);
  SWAP(fJetEtDiff_max, fJetRecHitsEtDiff_max);
}

// -------------------------------------------------------

int Cuts_t::passCuts(const DijetRespCorrDatumExtended_t &d) const {
  int fail=0;;
  if (!fail && (fJetEtDiff_max>0)) {
    double etDiff= fabs(d.CalcEtDiff());
    if (etDiff < fJetEtDiff_min) fail=1;
    if (etDiff > fJetEtDiff_max) fail=1;
  }

  if (!fail && (fJetRecHitsEtDiff_max>0)) {
    double etRHDiff= fabs(d.CalcRecHitsEtDiff());
    if (etRHDiff < fJetRecHitsEtDiff_min) fail=1;
    else if (etRHDiff > fJetRecHitsEtDiff_max) fail=1;
  }

  if (!fail && (fDijet_dPhi_max>0)) {
    double dPhi=fabs(d.GetDPhi());
    //std::cout << "dPhi=" << dPhi << "\n";
    if (fDijet_dPhi_min > dPhi) fail=1;
    else if (fDijet_dPhi_max < dPhi) fail=2;
    //if (fail) std::cout << "fails=" << fail << " (" << fDijet_dPhi_min << " .. " << fDijet_dPhi_max << ")\n";
  }

  if (!fail) {
    double ratio= d.CalcRecHits_EnergyDiff(_hadEtDiff_over_nonHadEtDiff);
    if (ratio< fDijetEtDiffRatio_min) fail=1;
    else if (ratio> fDijetEtDiffRatio_max) fail=2;
  }

  if (fail) return 0;
  return 1;
}

// -------------------------------------------------------
// -------------------------------------------------------
/*
Analyzer_t::Analyzer_t(TString name, const Cuts_t &cuts) :
  fName(name), fCuts(cuts),
  fhTow(NULL)
{
  if (!this->create()) std::cout << "error in Analyzer contructor\n";
}

// -------------------------------------------------------

int Analyzer_t::create() {
  fhTow=new TH1D("hTow_" + fName, "hTow_" + fName,
		 NUMTOWERS,0,NUMTOWERS+1);
  return 1;
}

// -------------------------------------------------------

int Analyzer_t::Fill(const DijetRespCorrDatumExtended_t &d) {
  if (!fCuts.passCuts(d)) return 0;
  return 1;
}
*/
// -------------------------------------------------------
// -------------------------------------------------------

void FillTowerCount(TH1D* histo, const DijetRespCorrDatum &d, int jetId,
		    int useWeight) {
  double w=(useWeight) ? d.GetWeight() : 1;
  for (int iJet=0; iJet<1; ++iJet) {
    if ((iJet==jetId) || (jetId==-1)) {
      typedef std::map<Int_t,Double_t> TMap_t;
      TMap_t heMap;
      if (iJet==0) d.GetTagHcalE(heMap);
      else d.GetProbeHcalE(heMap);
      for (TMap_t::const_iterator it= heMap.begin();
	   it!=heMap.end(); it++) {
	histo->Fill(it->first,w);
      }
    }
  }
  return;
}

// -------------------------------------------------------
// -------------------------------------------------------
