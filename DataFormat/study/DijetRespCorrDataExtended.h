#ifndef DijetRespCorrDataExtended_H
#define DijetRespCorrDataExtended_H

//#include "../interface/DijetRespCorrData.h"
#include <iostream>
#include "link_DijetRespCorrData.h"

// ------------------------------------------------------

typedef enum { _sumHadEt, _sumNonHadEt,
	       _sumHadEt_divAvg, _sumNonHadEt_divAvg,
	       _hadEtDiff_over_nonHadEtDiff
} TEnergyDiff_t;

// ------------------------------------------------------

class DijetRespCorrDatumExtended_t : public DijetRespCorrDatum {
 public:
  Long64_t fOrigIEntry;
  Float_t fDijet_dEta,fDijet_dPhi;
  Float_t fTPFJetE, fTPFJetP;
  Float_t fPPFJetE, fPPFJetP;
  Float_t fTJet_SumHadE, fTJet_SumHadOtherE, fTJet_SumNonHadOtherE;
  Float_t fPJet_SumHadE, fPJet_SumHadOtherE, fPJet_SumNonHadOtherE;
  Float_t fTPFJetGenE, fTPFJetGenPt;
  Float_t fPPFJetGenE, fPPFJetGenPt;
  std::vector<int> *fTPfJet_had_twr_ieta;
  std::vector<int> *fTPfJet_had_twr_iphi;
  std::vector<float> *fTPfJet_had_twr_hade;
  std::vector<float> *fTPfJet_had_twr_frac;
  std::vector<float> *fTPfJet_had_EcalE;
  std::vector<float> *fTPfJet_had_EcalEcorrection;

  std::vector<int> *fPPfJet_had_twr_ieta;
  std::vector<int> *fPPfJet_had_twr_iphi;
  std::vector<float> *fPPfJet_had_twr_hade;
  std::vector<float> *fPPfJet_had_twr_frac;
  std::vector<float> *fPPfJet_had_EcalE;
  std::vector<float> *fPPfJet_had_EcalEcorrection;
 public:
  DijetRespCorrDatumExtended_t();
  DijetRespCorrDatumExtended_t(const DijetRespCorrDatumExtended_t &d);
  ~DijetRespCorrDatumExtended_t();
  void Clear();

 private:
  int __allocate(); // create the pointers

 public:
  Long64_t GetOrigIEntry() const { return fOrigIEntry; }
  Float_t GetDPhi() const { return fDijet_dPhi; }
  Float_t GetDEta() const { return fDijet_dEta; }
  Float_t CalcEnDiff() const { return fTPFJetE-fPPFJetE; }

  Float_t CalcEt_TagJet() const { return fTPFJetE/cosh(this->GetTagEta()); }
  Float_t CalcEt_ProbeJet() const { return fPPFJetE/cosh(this->GetProbeEta());}
  Float_t CalcEtDiff() const
  { return (CalcEt_TagJet() - CalcEt_ProbeJet()); }
  Float_t CalcEtSum() const
  { return (CalcEt_TagJet() + CalcEt_ProbeJet()); }

  Float_t CalcRecHits_TagJetEn() const
  { return (fTJet_SumHadE + fTJet_SumHadOtherE + fTJet_SumNonHadOtherE); }
  Float_t CalcRecHits_ProbeJetEn() const
  { return (fPJet_SumHadE + fPJet_SumHadOtherE + fPJet_SumNonHadOtherE); }
  Float_t CalcRecHitsEn(int tagJet) const
  { return (tagJet) ? CalcRecHits_TagJetEn() : CalcRecHits_ProbeJetEn(); }
  Float_t CalcRecHitsEnDiff() const
  { return (CalcRecHits_TagJetEn() - CalcRecHits_ProbeJetEn()); }

  Float_t CalcRecHits_TagJetEt() const
  { return CalcRecHits_TagJetEn()/cosh(this->GetTagEta()); }
  Float_t CalcRecHits_ProbeJetEt() const
  { return CalcRecHits_ProbeJetEn()/cosh(this->GetProbeEta()); }
  Float_t CalcRecHitsEt(int tagJet) const
  { return (tagJet) ? CalcRecHits_TagJetEt() : CalcRecHits_ProbeJetEt(); }
  Float_t CalcRecHitsEtDiff() const
  { return (CalcRecHits_TagJetEt() - CalcRecHits_ProbeJetEt()); }
  Float_t CalcRecHitsEtSum() const
  { return (CalcRecHits_TagJetEt() + CalcRecHits_ProbeJetEt()); }

  Float_t Calc_EtBalance() const
  { return 2.*CalcEtDiff()/CalcEtSum(); }
  Float_t Calc_RecHitsEtBalance() const
  { return 2.*CalcRecHitsEtDiff()/CalcRecHitsEtSum(); }

  Float_t CalcRecHitsEnDiscrepancy_TagJet() const
  { return (fTPFJetE - CalcRecHits_TagJetEn()); }
  Float_t CalcRecHitsEnDiscrepancy_ProbeJet() const
  { return (fPPFJetE - CalcRecHits_ProbeJetEn()); }
  Float_t CalcRecHitsEnDiscrepancy(int tagJet) const
  { return (tagJet) ?
      CalcRecHitsEnDiscrepancy_TagJet():CalcRecHitsEnDiscrepancy_ProbeJet(); }

  Float_t CalcRecHits_EnergyDiff(TEnergyDiff_t enDiff) const;

  void SetOrigIEntry(Long64_t iEntry) { fOrigIEntry=iEntry; }
  void SetDijetDEtaDPhi(float dEta, float dPhi) {
    fDijet_dEta=dEta; fDijet_dPhi=dPhi;
  }

  void SetJetEPEtaPhi(int tagJet, float E, float P, float eta, float phi);

  void SetJetTowers(int tagJet,
		    double NonHadronic_EcalE,
		    const std::vector<int> &vec_twr_iEta,
		    const std::vector<int> &vec_twr_iPhi,
		    const std::vector<float> &vec_twr_hadE,
		    const std::vector<float> &vec_twr_frac,
		    const std::vector<float> &vec_had_EcalE,
		    // auxiliary quantities
		    const std::vector<int> &vec_twr_ClusterInd,
		    const std::vector<float> &vec_ci_Cluster_dR,
		    const std::vector<int>   &vec_had_Id,
		    const std::vector<int>   &vec_had_nTwrs,
		    const std::vector<int>   &vec_had_candtrackind,
		    const std::vector<float> &vec_cti_candtrack_px,
		    const std::vector<float> &vec_cti_candtrack_py,
		    const std::vector<float> &vec_cti_candtrack_pz
		    );

  void SetGenInfo(Float_t lead_genE, Float_t lead_genPt,
		  Float_t sublead_genE, Float_t sublead_genPt);

  int PassCuts(double minSumJetEt=40, double minJetEt=20,
	       double maxThirdJetEt=15, double maxDeltaEta=0.5,
	       double maxDeltaPhi=999.99, int *failFlag=NULL) const;

  void PrintBaseClass(int printExtra=0) const;

  friend std::ostream& operator<<(std::ostream& out,
				  const DijetRespCorrDatumExtended_t &d);

  ClassDef(DijetRespCorrDatumExtended_t, 1)
};


// ------------------------------------------------------
// ------------------------------------------------------

#endif
