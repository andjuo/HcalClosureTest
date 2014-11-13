#include "DijetRespCorrDataExtended.h"
#include "helper.hh"
#include <TString.h>
//#include <cmath>

#ifdef __CINT__
#pragma link C++ class DijetRespCorrDatumExtended_t;
#endif



// ----------------------------------------------------------------

DijetRespCorrDatumExtended_t::DijetRespCorrDatumExtended_t() :
  DijetRespCorrDatum()
{
  if (!__allocate())
    std::cout << "error in DijetRespCorrDatumExtended_t" << std::endl;
  else this->Clear();
}

// ----------------------------------------------------------------

DijetRespCorrDatumExtended_t::DijetRespCorrDatumExtended_t(const DijetRespCorrDatumExtended_t &d) :
  DijetRespCorrDatum(d)
{
  if (!__allocate())
    std::cout << "error in DijetRespCorrDatumExtended_t" << std::endl;
  std::cout << "DijetRespCorrDatumExtended_t(same) lacks full copying\n";
}

// ----------------------------------------------------------------

DijetRespCorrDatumExtended_t::~DijetRespCorrDatumExtended_t() {
  if (fTPfJet_had_twr_ieta) delete fTPfJet_had_twr_ieta;
  if (fTPfJet_had_twr_iphi) delete fTPfJet_had_twr_iphi;
  if (fTPfJet_had_twr_hade) delete fTPfJet_had_twr_hade;
  if (fTPfJet_had_twr_frac) delete fTPfJet_had_twr_frac;
  if (fTPfJet_had_EcalE) delete fTPfJet_had_EcalE;
  if (fTPfJet_had_EcalEcorrection) delete fTPfJet_had_EcalEcorrection;

  if (fPPfJet_had_twr_ieta) delete fPPfJet_had_twr_ieta;
  if (fPPfJet_had_twr_iphi) delete fPPfJet_had_twr_iphi;
  if (fPPfJet_had_twr_hade) delete fPPfJet_had_twr_hade;
  if (fPPfJet_had_twr_frac) delete fPPfJet_had_twr_frac;
  if (fPPfJet_had_EcalE) delete fPPfJet_had_EcalE;
  if (fPPfJet_had_EcalEcorrection) delete fPPfJet_had_EcalEcorrection;

}

// ----------------------------------------------------------------

int DijetRespCorrDatumExtended_t::__allocate() {
  fTPfJet_had_twr_ieta = new std::vector<int>();
  fTPfJet_had_twr_iphi = new std::vector<int>();
  fTPfJet_had_twr_hade = new std::vector<float>();
  fTPfJet_had_twr_frac = new std::vector<float>();
  fTPfJet_had_EcalE = new std::vector<float>();
  fTPfJet_had_EcalEcorrection= new std::vector<float>();
  fPPfJet_had_twr_ieta = new std::vector<int>();
  fPPfJet_had_twr_iphi = new std::vector<int>();
  fPPfJet_had_twr_hade = new std::vector<float>();
  fPPfJet_had_twr_frac = new std::vector<float>();
  fPPfJet_had_EcalE = new std::vector<float>();
  fPPfJet_had_EcalEcorrection= new std::vector<float>();
  int ok= (fTPfJet_had_twr_ieta && fTPfJet_had_twr_iphi &&
	   fTPfJet_had_twr_hade && fTPfJet_had_twr_frac &&
	   fTPfJet_had_EcalE && fTPfJet_had_EcalEcorrection) ? 1:0;
  if (ok) {
    ok= (fPPfJet_had_twr_ieta && fPPfJet_had_twr_iphi &&
	 fPPfJet_had_twr_hade && fPPfJet_had_twr_frac &&
	 fPPfJet_had_EcalE && fPPfJet_had_EcalEcorrection) ? 1:0;
  }
  return ok;
}

// ----------------------------------------------------------------

void DijetRespCorrDatumExtended_t::Clear() {
  DijetRespCorrDatum::Clear();
  fWeight=0;
  fTagEta=0;
  fTagPhi=0;
  fTagHcalE.clear();
  fTagEcalE=0;
  fProbeEta=0;
  fProbePhi=0;
  fProbeHcalE.clear();
  fProbeEcalE=0;
  fThirdJetPx=0;
  fThirdJetPy=0;
  fCandTrackN=0;
  fCandTrackP.clear();
  fCandTrackEcalE.clear();
  fCandTrackHcalE.clear();

  fOrigIEntry=-999;
  fDijet_dEta=-999;
  fDijet_dPhi=-999;
  fTPFJetE=-999; fTPFJetP=-999;
  fPPFJetE=-999; fPPFJetP=-999;
  fTJet_SumHadE=-999;
  fTJet_SumHadOtherE=-999;
  fTJet_SumNonHadOtherE=-999;
  fPJet_SumHadE=-999;
  fPJet_SumHadOtherE=-999;
  fPJet_SumNonHadOtherE=-999;
  fTPFJetGenE=0;
  fTPFJetGenPt=0;
  fPPFJetGenE=0;
  fPPFJetGenPt=0;
  fTPfJet_had_twr_ieta->clear();
  fTPfJet_had_twr_iphi->clear();
  fTPfJet_had_twr_hade->clear();
  fTPfJet_had_twr_frac->clear();
  fTPfJet_had_EcalE->clear();
  fTPfJet_had_EcalEcorrection->clear();
  fPPfJet_had_twr_ieta->clear();
  fPPfJet_had_twr_iphi->clear();
  fPPfJet_had_twr_hade->clear();
  fPPfJet_had_twr_frac->clear();
  fPPfJet_had_EcalE->clear();
  fPPfJet_had_EcalEcorrection->clear();
  //std::cout << "\nafter Clear\n";
  //this->PrintBaseClass();
  //std::cout << "--------\n";
}

// ----------------------------------------------------------------

Float_t DijetRespCorrDatumExtended_t::CalcRecHits_EnergyDiff
  (TEnergyDiff_t enDiff) const
{
  double val=0;
  double t_hadEt=fTJet_SumHadE/cosh(this->GetTagEta());
  double p_hadEt=fPJet_SumHadE/cosh(this->GetProbeEta());
  double t_nonHadEt=
    (fTJet_SumHadOtherE+fTJet_SumNonHadOtherE)/cosh(this->GetTagEta());
  double p_nonHadEt=
    (fPJet_SumHadOtherE+fPJet_SumNonHadOtherE)/cosh(this->GetProbeEta());
  double t_totEt= t_hadEt + t_nonHadEt;
  double p_totEt= p_hadEt + p_nonHadEt;
  double totEt= t_totEt + p_totEt;

  switch(enDiff) {
  case _sumHadEt:
  case _sumHadEt_divAvg:
    {
      val= (t_hadEt - p_hadEt);
      if (enDiff==_sumHadEt_divAvg) val/= totEt;
    }
    break;
  case _sumNonHadEt:
  case _sumNonHadEt_divAvg:
    {
      val= (t_nonHadEt - p_nonHadEt);
      if (enDiff==_sumNonHadEt_divAvg) val/= totEt;
    }
    break;
  case _hadEtDiff_over_nonHadEtDiff:
    val= (t_hadEt - p_hadEt)/(t_nonHadEt - p_nonHadEt);
    break;
  }
  return Float_t(val);
}

// ----------------------------------------------------------------

void DijetRespCorrDatumExtended_t::SetJetEPEtaPhi(int tagJet,
	       float E, float P, float eta, float phi) {
  if (tagJet==1) {
    fTPFJetE=E; fTPFJetP=P;
    this->SetTagEta(eta);
    this->SetTagPhi(phi);
  }
  else {
    fPPFJetE=E; fPPFJetP=P;
    this->SetProbeEta(eta);
    this->SetProbePhi(phi);
  }
}


// ----------------------------------------------------------------

void DijetRespCorrDatumExtended_t::SetJetTowers
           (int tagJet,
	    double NonHadronic_EcalE,
	    const std::vector<int> &vec_twr_iEta,
	    const std::vector<int> &vec_twr_iPhi,
	    const std::vector<float> &vec_twr_hadE,
	    const std::vector<float> &vec_twr_frac,
	    const std::vector<float> &vec_had_EcalE,
       // auxiliary quantities
	    const std::vector<int>   &vec_twr_ClusterInd,
	    const std::vector<float> &vec_ci_Cluster_dR,
	    const std::vector<int>   &vec_had_Id,
	    const std::vector<int>   &vec_had_nTwrs,
	    const std::vector<int>   &vec_had_candtrackind,
	    const std::vector<float> &vec_cti_candtrack_px,
	    const std::vector<float> &vec_cti_candtrack_py,
	    const std::vector<float> &vec_cti_candtrack_pz
						) {
  std::vector<float> *ecalCorr= 
    (tagJet) ? fTPfJet_had_EcalEcorrection : fPPfJet_had_EcalEcorrection;
  ecalCorr->clear();
  ecalCorr->reserve(vec_had_EcalE.size());
  for (unsigned int i=0; i<vec_had_EcalE.size(); ++i) ecalCorr->push_back(0);

  double sumHE=0;
  for (unsigned int iTwr=0; iTwr<vec_twr_iEta.size(); iTwr++) {
    const int ci=vec_twr_ClusterInd[iTwr];
    if ((vec_twr_hadE[iTwr]>0.) &&
	((ci<0) || (vec_ci_Cluster_dR[ci]<0.5))) {
      double term= vec_twr_hadE[iTwr] * vec_twr_frac[iTwr];
      sumHE += term;
      if (tagJet==1) this->AddTagHcalE( term, vec_twr_iEta[iTwr] );
      else this->AddProbeHcalE( term, vec_twr_iEta[iTwr] );
    }
  }

  double sumOE=0., sumENoRecHits=0.;
  for (unsigned int iHad=0; iHad<vec_had_EcalE.size(); ++iHad) {
    sumOE += vec_had_EcalE[iHad];
    const int candTrackIdx=vec_had_candtrackind[iHad];
    double corrTerm=0;
    if ( (vec_had_nTwrs[iHad]==0) && (candTrackIdx>-1)) {
      const int cti=candTrackIdx;
      double tpx= vec_cti_candtrack_px[cti];
      double tpy= vec_cti_candtrack_py[cti];
      double tpz= vec_cti_candtrack_pz[cti];
      corrTerm= sqrt( tpx*tpx + tpy*tpy + tpz*tpz ) - vec_had_EcalE[iHad];
      sumENoRecHits += corrTerm;
    }
    (*ecalCorr)[iHad] = corrTerm;
  }

  if (tagJet==1) {
    fTJet_SumHadE= sumHE;
    fTJet_SumHadOtherE= sumOE + sumENoRecHits;
    fTJet_SumNonHadOtherE= NonHadronic_EcalE;
    this->SetTagEcalE( fTJet_SumHadOtherE + fTJet_SumNonHadOtherE );
  }
  else {
    fPJet_SumHadE= sumHE;
    fPJet_SumHadOtherE= sumOE + sumENoRecHits;
    fPJet_SumNonHadOtherE= NonHadronic_EcalE;
    this->SetProbeEcalE( fPJet_SumHadOtherE + fPJet_SumNonHadOtherE );
  }

  if (tagJet==1) {
    (*fTPfJet_had_twr_ieta) = vec_twr_iEta;
    (*fTPfJet_had_twr_iphi) = vec_twr_iPhi;
    (*fTPfJet_had_twr_hade) = vec_twr_hadE;
    (*fTPfJet_had_twr_frac) = vec_twr_frac;
    (*fTPfJet_had_EcalE) = vec_had_EcalE;
  }
  else {
    (*fPPfJet_had_twr_ieta) = vec_twr_iEta;
    (*fPPfJet_had_twr_iphi) = vec_twr_iPhi;
    (*fPPfJet_had_twr_hade) = vec_twr_hadE;
    (*fPPfJet_had_twr_frac) = vec_twr_frac;
    (*fPPfJet_had_EcalE) = vec_had_EcalE;
  }
}

// ----------------------------------------------------------------------

void DijetRespCorrDatumExtended_t::SetGenInfo
    (Float_t lead_genE, Float_t lead_genPt,
     Float_t sublead_genE, Float_t sublead_genPt) {
  fTPFJetGenE = lead_genE; fTPFJetGenPt = lead_genPt;
  fPPFJetGenE = sublead_genE; fPPFJetGenPt = sublead_genPt;
}

// ----------------------------------------------------------------------

int DijetRespCorrDatumExtended_t::PassCuts
   (double minSumJetEt,
    double minJetEt,
    double maxThirdJetEt,
    double maxDeltaEta,
    double maxDeltaPhi,
    int *failFlag) const
{
  double tEta= this->GetTagEta();
  double pEta= this->GetProbeEta();
  double tEt= fTPFJetE / cosh(tEta);
  double pEt= fPPFJetE / cosh(pEta);

  int fail=0;
  if (tEt + pEt < minSumJetEt) fail |= 0x1;
  if ((tEt<minJetEt) || (pEt<minJetEt)) fail |= 0x2;
  if (sqrt(pow(this->GetThirdJetPx(),2) + pow(this->GetThirdJetPy(),2))
      > maxThirdJetEt) fail |= 0x4;
  if (fDijet_dEta > maxDeltaEta) fail |= 0x8;
  const double cPI=static_cast<double>(4)*atan(static_cast<double>(1));
  double deltaPhi=fabs(fDijet_dPhi)-cPI;
  if (fabs(deltaPhi) > maxDeltaPhi) fail |= 0x16;

  if (fDijet_dEta<0) std::cout << "Negative fDijet_dEta\n";
  //if (fDijet_dPhi<0) std::cout << "Negative fDijet_dPhi\n";

  if (failFlag) *failFlag = fail;
  return (fail) ? 0:1;
}

// ----------------------------------------------------------------

void DijetRespCorrDatumExtended_t::PrintBaseClass(int printExtra) const {
  std::cout << *(DijetRespCorrDatum*)this;
  if (printExtra) {
    std::cout << " extra info:\n";
    std::cout << "     origIEntry=" << fOrigIEntry << "\n";
    std::cout << "     tag   sumHadE=" << fTJet_SumHadE
	      << ", sumHadOtherE=" << fTJet_SumHadOtherE
	      << ", sumNonHadOtherE=" << fTJet_SumNonHadOtherE << "\n";
    std::cout << "     probe sumHadE=" << fPJet_SumHadE
	      << ", sumHadOtherE=" << fPJet_SumHadOtherE
	      << ", sumNonHadOtherE=" << fPJet_SumNonHadOtherE << "\n";
  }
}


// ----------------------------------------------------------------

std::ostream& operator<<(std::ostream& out,
			 const DijetRespCorrDatumExtended_t &d) {
  if (out != std::cout) {
    out << "cannot print object DijetRespCorrDatumExtended_t\n";
    return out;
  }

  out << "origIEntry=" << d.fOrigIEntry << "\n";
  out << "dijet_dEta=" << d.fDijet_dEta
      << ", dijet_dPhi=" << d.fDijet_dPhi << "\n";
  out << "tag jet E=" << d.fTPFJetE << ", p=" << d.fTPFJetP
      << ", eta=" << d.GetTagEta() << ", phi=" << d.GetTagPhi() << "\n";
  out << "     sumHadE=" << d.fTJet_SumHadE
      << ", sumHadOtherE=" << d.fTJet_SumHadOtherE
      << ", sumNonHadOtherE=" << d.fTJet_SumNonHadOtherE << "\n";
  out << "probe jet E=" << d.fPPFJetE << ", p=" << d.fPPFJetP
      << ", eta=" << d.GetProbeEta() << ", phi=" << d.GetProbePhi() << "\n";
  out << "     sumHadE=" << d.fPJet_SumHadE
      << ", sumHadOtherE=" << d.fPJet_SumHadOtherE
      << ", sumNonHadOtherE=" << d.fPJet_SumNonHadOtherE << "\n";
  out << "     jet1 genE=" << d.fTPFJetGenE << ", genPt="
      << d.fTPFJetGenPt << ", jet2: " << d.fPPFJetGenE << ", genPt="
      << d.fPPFJetGenPt << "\n";
  printVec("tagJet had_twr_ieta ", *d.fTPfJet_had_twr_ieta);
  printVec("tagJet had_twr_iphi ", *d.fTPfJet_had_twr_iphi);
  printVec("tagJet had_twr_hade ", *d.fTPfJet_had_twr_hade);
  printVec("tagJet had_twr_frac ", *d.fTPfJet_had_twr_frac);
  printVec("tagJet had_twr_EcalE ", *d.fTPfJet_had_EcalE);
  printVec("tagJet had_twr_EcalEcorr ", *d.fTPfJet_had_EcalEcorrection);

  printVec("probeJet had_twr_ieta ", *d.fPPfJet_had_twr_ieta);
  printVec("probeJet had_twr_iphi ", *d.fPPfJet_had_twr_iphi);
  printVec("probeJet had_twr_hade ", *d.fPPfJet_had_twr_hade);
  printVec("probeJet had_twr_frac ", *d.fPPfJet_had_twr_frac);
  printVec("probeJet had_twr_EcalE ", *d.fPPfJet_had_EcalE);
  printVec("probeJet had_twr_EcalEcorr ", *d.fPPfJet_had_EcalEcorrection);
  return out;
}


// ----------------------------------------------------------------
// ----------------------------------------------------------------

