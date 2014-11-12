#ifdef __localRun
#  include "../interface/GammaJetFitData.h"
#else
#  include "HcalClosureTest/DataFormat/interface/GammaJetFitData.h"
#endif

#ifdef __localRun
#ifdef __CINT__
#pragma link C++ class GammaJetEventAuxInfo_t;
#pragma link C++ class GammaJetEvent_t;
#pragma link C++ class GammaJetFitter_t;
#endif
#endif

// ----------------------------------------------------------------

TMinuit *myMinuit= NULL;

// ----------------------------------------------------------------
// ----------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const std::map<Int_t,Double_t> &m){
  out << "map.size=" << m.size() << ": ";
  for (std::map<Int_t,Double_t>::const_iterator it=m.begin();
       it!=m.end(); it++) {
    out << "(" << it->first << "," << it->second << ")";
  }
  return out;
}

// ----------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const std::vector<Double_t> &vec) {
  out << "vec[" << vec.size() << "]: ";
  for (unsigned int i=0; i<vec.size(); ++i) {
    out << " " << vec[i];
  }
  return out;
}

// ----------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const TArrayD &arr) {
  out << "arr[" << arr.GetSize() << "]: ";
  for (int i=0; i<arr.GetSize(); ++i) {
    out << " " << arr[i];
  }
  return out;
}

// ----------------------------------------------------------------

TArrayD convert(const std::vector<Double_t> &vec) {
  TArrayD arr(vec.size());
  for (unsigned int i=0; i<vec.size(); ++i) {
    arr[i]=vec[i];
  }
  return arr;
}

// ----------------------------------------------------------------
// ----------------------------------------------------------------

void GammaJetEventAuxInfo_t::Assign(const GammaJetEventAuxInfo_t &e) {
  fEventNo=e.fEventNo;
  fRunNo=e.fRunNo;
}

// ----------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const GammaJetEventAuxInfo_t &e) {
  out << "eventNo=" << e.fEventNo << ", runNo=" << e.fRunNo;
  return out;
}

// ----------------------------------------------------------------
// ----------------------------------------------------------------
// ----------------------------------------------------------------

GammaJetEvent_t::GammaJetEvent_t()
  : TObject(),
    fWeight(1.),
    fTagEta(0.), fTagPhi(0.), fTagEcalE(0.),
    fTagHcalE(),
    fTagE(0.), fTagEt(0.),
    fProbeEta(0.), fProbePhi(0.), fProbeEcalE(0.),
    fProbeHcalE(),
    fProbeE(0.), fProbeEt(0.),
    fAInfo()
{ }

// ----------------------------------------------------------------

GammaJetEvent_t::GammaJetEvent_t(const GammaJetEvent_t &e)
  : TObject(e),
    fWeight(e.fWeight),
    fTagEta(e.fTagEta), fTagPhi(e.fTagPhi), fTagEcalE(e.fTagEcalE),
    fTagHcalE(e.fTagHcalE),
    fTagE(e.fTagE), fTagEt(e.fTagEt),
    fProbeEta(e.fProbeEta), fProbePhi(e.fProbePhi), fProbeEcalE(e.fProbeEcalE),
    fProbeHcalE(e.fProbeHcalE),
    fProbeE(e.fProbeE), fProbeEt(e.fProbeEt),
    fAInfo(e.fAInfo)
{}

// ----------------------------------------------------------------

GammaJetEvent_t::~GammaJetEvent_t() {
  // use defaults
}

// ----------------------------------------------------------------

Double_t GammaJetEvent_t::GetTagHcalE(Int_t ieta) const {
  int found=0;
  double val=0.;
  for (std::map<Int_t,Double_t>::const_iterator it= fTagHcalE.begin();
       !found && (it!=fTagHcalE.end()); ++it) {
    if (it->first == ieta) {
      found=1;
      val=it->second;
    }
  }
  return val;
}

// ----------------------------------------------------------------

void GammaJetEvent_t::SetTagEEtaPhi(double valE, double valEta, double valPhi){
  fTagEta=valEta; fTagPhi=valPhi;
  fTagEcalE= valE;
  fTagHcalE.clear();
  fTagE= valE; fTagEt= valE/cosh(valEta);
}

// ----------------------------------------------------------------

void GammaJetEvent_t::SetTagEtaPhiEn(double valEta, double valPhi,
				     double valEcalE,
				     const std::map<Int_t,Double_t> &valHcalE){
  fTagEta=valEta; fTagPhi=valPhi;
  fTagEcalE= valEcalE;
  fTagHcalE.clear();
  fTagHcalE= valHcalE;
  fTagE= calc_HcalE(valHcalE);
  fTagE += valEcalE;
  fTagEt= fTagE/cosh(valEta);
}

// ----------------------------------------------------------------

Double_t GammaJetEvent_t::GetProbeHcalE(Int_t ieta) const {
  int found=0;
  double val=0.;
  for (std::map<Int_t,Double_t>::const_iterator it= fProbeHcalE.begin();
       !found && (it!=fProbeHcalE.end()); ++it) {
    if (it->first == ieta) {
      found=1;
      val=it->second;
    }
  }
  return val;
}

// ----------------------------------------------------------------

void GammaJetEvent_t::SetProbeEtaPhiEn(double valEta, double valPhi,
				       double valEcalE,
				     const std::map<Int_t,Double_t> &valHcalE){
  fProbeEta=valEta; fProbePhi=valPhi;
  fProbeEcalE= valEcalE;
  fProbeHcalE.clear();
  fProbeHcalE= valHcalE;
  fProbeE= calc_HcalE(valHcalE);
  fProbeE += valEcalE;
  fProbeEt= fProbeE/cosh(valEta);
}

// ----------------------------------------------------------------

void GammaJetEvent_t::MyPrint(int detail, std::ostream& out) const {
  out << fAInfo << "\n";
  out << "tagEt=" << fTagEt << ", probeEt=" << fProbeEt << "\n";
  if (detail) {
    // print tag info
    out << "tag(Eta,Phi) = (" << fTagEta << "," << fTagPhi << ")\n";
    out << "tag(EcalE,HcalE, E) = (" << fTagEcalE << ","
	<< (fTagE - fTagEcalE) << ", " << fTagE << "\n";
    if (detail==2) {
      out << "tagHcalE " << fTagHcalE << ": sum="
	  << calc_HcalE(fTagHcalE) << "\n";
    }
    // print probe info
    out << "probe(Eta,Phi) = (" << fProbeEta << "," << fProbePhi << ")\n";
    out << "probe(EcalE,HcalE, E) = (" << fProbeEcalE << ","
	<< (fProbeE - fProbeEcalE) << ", " << fProbeE << "\n";
    if (detail==2) {
      out << "probeHcalE " << fProbeHcalE << ": sum="
	  << calc_HcalE(fProbeHcalE) << "\n";
    }
  }
}

// ----------------------------------------------------------------
// ----------------------------------------------------------------

GammaJetCuts_t::GammaJetCuts_t()
  : fEtDiffMin(-1.), fEtDiffMax(1e16)
{}

// ----------------------------------------------------------------

GammaJetCuts_t::GammaJetCuts_t(const GammaJetCuts_t &c)
  : fEtDiffMin(c.fEtDiffMin), fEtDiffMax(c.fEtDiffMax)
{}

// ----------------------------------------------------------------

int GammaJetCuts_t::passCuts(const GammaJetEvent_t &e) const {
  int ok=1;
  double tagProbeEtDiff= fabs(e.CalcDiffEt());
  if ((tagProbeEtDiff<fEtDiffMin) || (tagProbeEtDiff>fEtDiffMax)) ok=0;
  return ok;
}

// ----------------------------------------------------------------
// ----------------------------------------------------------------

GammaJetFitter_t::GammaJetFitter_t()
  : fData(),
    fFittingProcedure(0),
    fResolution(0.1),
    fPrintLevel(0),
    fErrorDef(0.5),
    fParStep(0.1), fParMin(0.), fParMax(2.),
    fEcalRes(0.07), fHcalRes(1.15), fHfRes(1.35)
{ }

// ----------------------------------------------------------------

GammaJetFitter_t::~GammaJetFitter_t() {
  for (unsigned int i=0; i<fData.size(); ++i) delete fData[i];
  fData.clear();
}

// ----------------------------------------------------------------

Double_t GammaJetFitter_t::EstimateResolution(int setValue) {
  if (setValue) fResolution=0.;
  double sumW=0, sumWX=0, sumWXX=0.;
  int allZeros=1;
  for (unsigned int i=0; i<fData.size(); ++i) {
    const GammaJetEvent_t *e= fData[i];
    Double_t w= e->GetWeight();
    Double_t x= this->CalcFitValue(e);
    if (allZeros && (x!=Double_t(0))) allZeros=0;
    sumW += w;
    sumWX += w*x;
    sumWXX += w*x*x;
  }
  if (allZeros) {
    std::cout << "EstimateResolution: allZeros detected" << std::endl;
    return 0.;
  }
  Double_t sigmaSqr= sumWXX/sumW - pow(sumWX/sumW,2);
  if (sigmaSqr<=Double_t(0)) {
    std::cout << "EstimateResolution: sigmaSqr<=0" << std::endl;
    return 0.;
  }
  Double_t resolution= sqrt(sigmaSqr);
  if (setValue) fResolution = resolution;
  return resolution;
}


// ----------------------------------------------------------------

Double_t GammaJetFitter_t::GetChi2(const TArrayD &hcalCorrCf) const {
  double chi2=0;

  //std::cout << "GetChi2 : ";
  //std::cout << "hcalCorrCf " << hcalCorrCf << "\n";

  for (unsigned int i=0; i<fData.size(); ++i) {
    const GammaJetEvent_t *e= fData[i];
    double w = e->GetWeight();
    double diff= e->CalcDiffEt(hcalCorrCf);
    //std::cout << "i=" << i << ", diff=" << diff << ", w=" << w << "\n";
    chi2 += pow(diff*w/fResolution,2);
  }
  //std::cout << "chi2=" << chi2 << "\n";
  return chi2;
}

// ----------------------------------------------------------------

Double_t GammaJetFitter_t::GetChi2_jj(const TArrayD &hcalCorrCf) const {
  double chi2=0;
  for (unsigned int i=0; i<fData.size(); ++i) {
    const GammaJetEvent_t *e= fData[i];
    double w = e->GetWeight();
    double diff= e->CalcDiffEt_jj(hcalCorrCf);
    chi2 += pow(diff*w/fResolution,2);
  }
  return chi2;
}

// ----------------------------------------------------------------

Double_t GammaJetFitter_t::GetFitValue(const GammaJetEvent_t *e,
				       const TArrayD &hcalCorrCf) const {
  Double_t val=0;
  switch(fFittingProcedure) {
  case 0: val= e->CalcDiffEt(hcalCorrCf); break;
  case 1: val= e->CalcDiffEt_jj(hcalCorrCf); break;
  default:
    std::cout << "GetFitValue(GammaJetEvent,TArrayD) "
	      << "is not ready for fittingProcedure="
	      << fFittingProcedure << std::endl;
  }
  return val;
}

// ----------------------------------------------------------------

Double_t GammaJetFitter_t::GetFitValue(const TArrayD &hcalCorrCf) const {
  Double_t val=0;
  switch(fFittingProcedure) {
  case 0: val=GetChi2(hcalCorrCf); break;
  case 1: val=GetChi2_jj(hcalCorrCf); break;
  default:
    std::cout << "GetFitValue(TArrayD) is not ready for fittingProcedure="
	      << fFittingProcedure << std::endl;
  }
  //std::cout << "GetFitValue=" << val << "\n";
  return val;
}

// ----------------------------------------------------------------

TH1D* GammaJetFitter_t::doFit(const char *histoName, const char *histoTitle,
			      std::vector<Double_t> &hcalCorrCf,
			      const std::vector<Int_t> *fixTowers) {
  int corrArrSize=NUMTOWERS;

  if (this->EstimateResolution()==Double_t(0)) return NULL;

  // Create the Minuit object
  myMinuit = new TMinuit(corrArrSize);
  myMinuit->SetPrintLevel(fPrintLevel);
  myMinuit->SetErrorDef(fErrorDef);
  myMinuit->SetFCN(gammaJet_FCN);
  myMinuit->SetObjectFit(this);

  // define the parameters
  for (int i=0; i<NUMTOWERS; i++) {
    int ieta= -MAXIETA + i;
    myMinuit->DefineParameter(i, Form("Tower ieta: %d",ieta),
			     hcalCorrCf[i],
			     fParStep, fParMin, fParMax);
    if (fixTowers &&
	(std::find(fixTowers->begin(),fixTowers->end(), ieta) !=
	 fixTowers->end())) {
      std::cout << " -- fix parameter\n";
      myMinuit->FixParameter(i);
    }
  }

  // scan the parameters
  //myMinuit->mcscan();

  // Minimize
  myMinuit->Migrad();

  // return the results
  TH1D* histo=new TH1D(histoName,histoTitle,
		       NUMTOWERS,-MAXIETA-0.5,MAXIETA+0.5);
  for(int i=1; i<=NUMTOWERS; i++) {
    Double_t val, error;
    myMinuit->GetParameter(i-1, val, error);
    histo->SetBinContent(i, val);
    histo->SetBinError(i, error);
  }

  for (unsigned int i=0; i<hcalCorrCf.size(); ++i) {
    Double_t val,error;
    myMinuit->GetParameter(i, val,error);
    hcalCorrCf[i]=val;
    if (int(i)>=NUMTOWERS) {
      std::cout << "spec param #" << i << " " << val << " +- " << error <<"\n";
    }
  }

  return histo;

}

// ----------------------------------------------------------------

void GammaJetFitter_t::PrintFitValueChange(
            const std::vector<Double_t> &finalCf,
	    const std::vector<Double_t> *iniCf) const {
  double iniVal=0;
  if (iniCf) {
    TArrayD cf= convert(*iniCf);
    iniVal= GetFitValue(cf);
  }
  else {
    TArrayD cf(finalCf.size());
    cf.Reset(1.);
    iniVal= GetFitValue(cf);
  }
  TArrayD cf1= convert(finalCf);
  double finVal= GetFitValue(cf1);
  std::cout << "fit value changed from " << iniVal << " to " << finVal << "\n";
}

// ----------------------------------------------------------------


void gammaJet_FCN(Int_t &npar, Double_t* gin, Double_t &f,
		  Double_t *par, Int_t iflag) {
  if (0 && (npar==0)) std::cout << "npar=0\n"; // get rid of compiler msg
  if (0 && iflag) std::cout << "iflag=" << iflag << "\n";
  if (!gin) std::cout << "gin is NULL\n";
  const GammaJetFitter_t* fitter=
    dynamic_cast<const GammaJetFitter_t*>(myMinuit->GetObjectFit());
  TArrayD corrCf;
  int arrSize= fitter->GetParameterCount();
  corrCf.Set(arrSize, par);
  //std::cout << "corrCf=" << corrCf << "\n";
  f = fitter->GetFitValue(corrCf);
  //std::cout << "f=" << f << "\n";
  return;
}

// ----------------------------------------------------------------
// ----------------------------------------------------------------
// ----------------------------------------------------------------
