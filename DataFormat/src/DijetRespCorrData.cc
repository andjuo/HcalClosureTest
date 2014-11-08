#ifndef __localRun
#  include "HcalClosureTest/DataFormat/interface/DijetRespCorrData.h"
#else
#  include "../interface/DijetRespCorrData.h"
#endif

#include "TMinuit.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"

#include <iostream>
#include <sstream>
#include <cassert>
#include <cmath>

ClassImp(DijetRespCorrDatum)
ClassImp(DijetRespCorrData)

TMinuit *xMinuit=NULL;

DijetRespCorrDatum::DijetRespCorrDatum() {
  fTagEta=fProbeEta=fTagEcalE=fProbeEcalE=0.0;
}

DijetRespCorrDatum::~DijetRespCorrDatum() {}

Double_t DijetRespCorrDatum::GetWeight(void) const
{
  return fWeight;
}

Double_t DijetRespCorrDatum::GetTagEta(void) const
{
  return fTagEta;
}

Double_t DijetRespCorrDatum::GetTagPhi(void) const
{
  return fTagPhi;
}

Double_t DijetRespCorrDatum::GetTagHcalE(Int_t ieta)
{
  Double_t v=fTagHcalE[ieta];
  return v;
}

void DijetRespCorrDatum::GetTagHcalE(std::map<Int_t, Double_t>& m) const
{
  m=fTagHcalE;
  return;
}

Double_t DijetRespCorrDatum::GetTagEcalE(void) const
{
  return fTagEcalE;
}

Double_t DijetRespCorrDatum::GetProbeEta(void) const
{
  return fProbeEta;
}

Double_t DijetRespCorrDatum::GetProbePhi(void) const
{
  return fProbePhi;
}

Double_t DijetRespCorrDatum::GetProbeHcalE(Int_t ieta)
{
  Double_t v=fProbeHcalE[ieta];
  return v;
}

void DijetRespCorrDatum::GetProbeHcalE(std::map<Int_t, Double_t>& m) const
{
  m=fProbeHcalE;
  return;
}

Double_t DijetRespCorrDatum::GetProbeEcalE(void) const
{
  return fProbeEcalE;
}

Double_t DijetRespCorrDatum::GetThirdJetPx(void) const
{
  return fThirdJetPx;
}

Double_t DijetRespCorrDatum::GetThirdJetPy(void) const
{
  return fThirdJetPy;
}

Int_t DijetRespCorrDatum::GetCandTrackN(void) const
{
  return fCandTrackN;
}

void DijetRespCorrDatum::SetWeight(Double_t v)
{
  fWeight = v;
  return;
}

void DijetRespCorrDatum::SetTagEta(Double_t v)
{
  fTagEta = v;
  return;
}

void DijetRespCorrDatum::SetTagPhi(Double_t v)
{
  fTagPhi = v;
  return;
}

void DijetRespCorrDatum::SetTagHcalE(Double_t v, Int_t ieta)
{
  assert(ieta<=MAXIETA && ieta>=-MAXIETA && ieta!=0);
  fTagHcalE[ieta] = v;
  return;
}

void DijetRespCorrDatum::AddTagHcalE(Double_t v, Int_t ieta)
{
  assert(ieta<=MAXIETA && ieta>=-MAXIETA && ieta!=0);
  fTagHcalE[ieta] += v;
  return;
}

void DijetRespCorrDatum::SetTagEcalE(Double_t v)
{
  fTagEcalE = v;
  return;
}

void DijetRespCorrDatum::SetProbeEta(Double_t v)
{
  fProbeEta = v;
  return;
}

void DijetRespCorrDatum::SetProbePhi(Double_t v)
{
  fProbePhi = v;
  return;
}


void DijetRespCorrDatum::SetProbeHcalE(Double_t v, Int_t ieta)
{
  assert(ieta<=MAXIETA && ieta>=-MAXIETA && ieta!=0);
  fProbeHcalE[ieta] = v;
  return;
}

void DijetRespCorrDatum::AddProbeHcalE(Double_t v, Int_t ieta)
{
  assert(ieta<=MAXIETA && ieta>=-MAXIETA && ieta!=0);
  fProbeHcalE[ieta] += v;
  return;
}

void DijetRespCorrDatum::SetProbeEcalE(Double_t v)
{
  fProbeEcalE = v;
  return;
}

void DijetRespCorrDatum::SetThirdJetPx(Double_t v)
{
  fThirdJetPx = v;
  return;
}

void DijetRespCorrDatum::SetThirdJetPy(Double_t v)
{
  fThirdJetPy = v;
  return;
}

void DijetRespCorrDatum::SetCandTrackN(Int_t v)
{
  fCandTrackN = v;
  return;
}

void DijetRespCorrDatum::AddCandTrackP(Double_t v)
{
  fCandTrackP.push_back(v);
  return;
}

void DijetRespCorrDatum::AddCandTrackEcalE(Double_t v)
{
  fCandTrackEcalE.push_back(v);
  return;
}

void DijetRespCorrDatum::AddCandTrackHcalE(std::map<Int_t, Double_t> v)
{
  fCandTrackHcalE.push_back(v);
  return;
}

void DijetRespCorrDatum::GetTagEnergies(const TArrayD& respcorr, Double_t& ecal, Double_t& hcal, Double_t& hf) const
{
  ecal=GetTagEcalE();
  hcal=0.0;
  hf=0.0;

  std::map<Int_t,Double_t> tagmap;
  GetTagHcalE(tagmap);
  for(std::map<Int_t, Double_t>::const_iterator mapit=tagmap.begin(); mapit!=tagmap.end(); ++mapit) {
    int ieta=mapit->first;
    double energy=mapit->second;
    int index=ieta+MAXIETA;
    if(std::abs(ieta)>29)
      hf += respcorr[index]*energy;
    else
      hcal += respcorr[index]*energy;
  }
  return;
}

void DijetRespCorrDatum::GetProbeEnergies(const TArrayD& respcorr, Double_t& ecal, Double_t& hcal, Double_t& hf) const
{
  // calculate the ecal, hcal, and HF energy
  // scale the hcal and hf energy by the response corrections
  ecal=GetProbeEcalE();
  hcal=0.0;
  hf=0.0;

  std::map<Int_t,Double_t> probemap;
  GetProbeHcalE(probemap);
  for(std::map<Int_t, Double_t>::const_iterator mapit=probemap.begin(); mapit!=probemap.end(); ++mapit) {
    int ieta=mapit->first;
    double energy=mapit->second;
    int index=ieta+MAXIETA;
    if(std::abs(ieta)>29)
      hf += respcorr[index]*energy;
    else
      hcal += respcorr[index]*energy;
  }
  return;
}

void DijetRespCorrDatum::GetTrackVariables(const TArrayD& respcorr, const Int_t index_, Double_t& TrackP_, Double_t& EcalE_, Double_t& HcalE_) const
{
  TrackP_ = fCandTrackP[index_];
  EcalE_ = fCandTrackEcalE[index_];

  std::map<Int_t,Double_t> candmap = fCandTrackHcalE[index_];
  for(std::map<Int_t, Double_t>::const_iterator mapit=candmap.begin(); mapit!=candmap.end(); ++mapit) {
    int ieta=mapit->first;
    double energy=mapit->second;
    int index=ieta+MAXIETA;
    HcalE_ += respcorr[index]*energy;
  }
  return;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

DijetRespCorrData::DijetRespCorrData()
{
  fData.clear();
  fPrintLevel=5;
  fParStep=0.10;
  fParMin=0.0;
  fParMax=0.0;//0.0;
  fEcalRes=0.07;
  fHcalRes=1.15;
  fHfRes=1.35;
  fBalanceSigma=1.;

  fDoCandTrackEnergyDiff=0;
}

DijetRespCorrData::~DijetRespCorrData()
{
}

void DijetRespCorrData::push_back(const DijetRespCorrDatum& d)
{
  fData.push_back(d);
  return;
}

const DijetRespCorrDatum& DijetRespCorrData::GetAt(Int_t i) const
{
  return fData[i];
}

Int_t DijetRespCorrData::GetSize(void) const
{
  return fData.size();
}

Double_t DijetRespCorrData::GetLikelihoodDistance(const TArrayD& respcorr) const
{
  Double_t total=0.0;
  //std::cout << "GetLikelihoodDistance" << std::endl;

  if(GetDoCandTrackEnergyDiff()==1){
    // loop over each jet pair
    for(std::vector<DijetRespCorrDatum>::const_iterator it=fData.begin(); it!=fData.end(); ++it) {
      
      // calculate the balance and resolution for each jet pair
      Double_t B, dB;
      GetBalance(*it, respcorr, B, dB);
      
      // this is the total likelihood
      total += 0.5*(std::log(dB*dB)+B*B/dB/dB);
      for(Int_t iCandTrack=0; iCandTrack<(*it).GetCandTrackN(); iCandTrack++){
	// Calculate the difference in energy between rechits and tracks
	Double_t Ediff, dEdiff;
	GetTrackEnergyDiff(*it, iCandTrack, respcorr, Ediff, dEdiff);
	total += 0.5*(std::log(dEdiff*dEdiff) + Ediff*Ediff/dEdiff/dEdiff);
      }
    }
  }
  else if (!GetDoCandTrackEnergyDiff()) {
    // loop over each jet pair
    for(std::vector<DijetRespCorrDatum>::const_iterator it=fData.begin(); it!=fData.end(); ++it) {
      
      // calculate the balance and resolution for each jet pair
      Double_t B, dB;
      GetBalance(*it, respcorr, B, dB);
      
      // this is the total likelihood
      total += 0.5*(std::log(dB*dB)+B*B/dB/dB);
    }
  }
  else if (GetDoCandTrackEnergyDiff()==2) {
    // loop over each jet pair
    for(std::vector<DijetRespCorrDatum>::const_iterator it=fData.begin(); it!=fData.end(); ++it) {

      // calculate the balance and resolution for each jet pair
      Double_t B, dB;
      GetBalance(*it, respcorr, B, dB);

      // this is the total likelihood
      dB= fBalanceSigma;
      total += 0.5*(std::log(dB*dB) + B*B/dB/dB) * pow(it->GetWeight(),2);
    }
    //std::cout << "total=" << total << " dB=" << fBalanceSigma << "\n";
  }
  else if ((GetDoCandTrackEnergyDiff()==3) ||
	   (GetDoCandTrackEnergyDiff()==4)) {
    // loop over each jet pair
    for(std::vector<DijetRespCorrDatum>::const_iterator it=fData.begin(); it!=fData.end(); ++it) {

      // calculate the balance and resolution for each jet pair
      Double_t B, dB, dummy;
      //GetBalance(*it, respcorr, B, dB);
      GetMET(*it, respcorr, B, dummy, dB);

      // this is the total likelihood
      dB= fBalanceSigma;
      total += 0.5*(std::log(dB*dB) + B*B/dB/dB) * pow(it->GetWeight(),2);
    }
    //std::cout << "total=" << total << " dB=" << fBalanceSigma << std::endl;
  }
  return total;
}

TH1D* DijetRespCorrData::doFit(const char* name, const char* title)
{
  TArrayD respcorr, respcorre;
  doFit(respcorr, respcorre);
  TH1D* h=new TH1D(name,title,NUMTOWERS,-MAXIETA-0.5,MAXIETA+0.5);
  for(int i=1; i<=NUMTOWERS; i++) {
    h->SetBinContent(i, respcorr[i-1]);
    h->SetBinError(i, respcorre[i-1]);
  }

  return h;
}

void DijetRespCorrData::doFit(TArrayD& respcorr, TArrayD& respcorre)
{
  // setup the initial response corrections
  const int maxIetaFixed=1;//20;
  //  Double_t array[NUMTOWERS] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.26437, 1.27885, 1.252, 1.26742, 1.32585, 1.27661, 1.29944, 1.451, 1.2652, 1.25045, 1.29709, 1.23643, 1.11458, 1.14601, 1.20513, 1.15064, 1.11951, 1.16586, 1.15791, 1.13728, 1.14483, 1.1412, 1.11142, 0, 1.15833, 1.14589, 1.15, 1.14048, 1.22407, 1.09756, 1.07979, 1.14484, 1.22885, 1.20833, 1.21161, 1.18929, 1.17783, 1.27585, 1.29167, 1.25481, 1.26563, 1.35172, 1.2816, 1.25988, 1.22321, 1.21111, 1.175, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
  Double_t array[NUMTOWERS] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
				//1.151630, 1.148890, 1.144870, 1.161240, 1.195630, 1.190690, 1.162400, 1.130190, 1.128810, 1.114070, 1.109850, 1.098620, 1.099950, 1.096140, 1.087350, 1.076980, 1.095540, 1.092960, 1.091500, 1.087100, 1.0, 1.085640, 1.094110, 1.089310, 1.089260, 1.079420, 1.089310, 1.093900, 1.104690, 1.098890, 1.109600, 1.115400, 1.144240, 1.125160, 1.148010, 1.204180, 1.197330, 1.150100, 1.154010, 1.143610, 1.159000,
				1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
				1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

  TArrayD respCorrInit;
  respCorrInit.Set(NUMTOWERS, array);

  // set the number of parameters to be the number of towers
  xMinuit=new TMinuit(NUMTOWERS);
  xMinuit->SetPrintLevel(fPrintLevel);
  //std::cout << "========================= " << xMinuit->GetMaxIterations() << std::endl;
  //xMinuit->SetMaxIterations(1000);
  xMinuit->SetErrorDef(0.5); // for a log likelihood
  xMinuit->SetFCN(FCN);
  xMinuit->SetObjectFit(this);

  // define the parameters
  for(int i=0; i<respCorrInit.GetSize(); i++) {
    int ieta=-MAXIETA+i;
    std::ostringstream oss;
    oss << "Tower ieta: " << ieta;
    xMinuit->DefineParameter(i, oss.str().c_str(), respCorrInit[i], fParStep, fParMin, fParMax);
    if(ieta>=-maxIetaFixed && ieta<=maxIetaFixed) xMinuit->FixParameter(i);

    // override the HF
    /*if(ieta<=-30 || ieta>=30) {
      xMinuit->DefineParameter(i, oss.str().c_str(), 1.0, fParStep, fParMin, fParMax);
      //xMinuit->DefineParameter(i, oss.str().c_str(), 1.3, fParStep, fParMin, fParMax);
      xMinuit->FixParameter(i);
      }*/
  }

  //gMinuit->Command("SET STR 2");

  // Minimize
  xMinuit->Migrad();

  // get the results
  TArrayD results(NUMTOWERS);
  TArrayD errors(NUMTOWERS);
  for(int i=0; i<results.GetSize(); i++) {
    Double_t val, error;
    xMinuit->GetParameter(i, val, error);
    results[i]=val;
    errors[i]=error;
  }
  respcorr=results;
  respcorre=errors;

  return;
}

TH1D* DijetRespCorrData::doFit_v2(const char* hname, const char* htitle,
				  std::vector<double> &respCorrInit,
				  const std::vector<int> &fixTowers
				  )
{
  int corrArrSize=NUMTOWERS;
  if (fDoCandTrackEnergyDiff==4) corrArrSize++;

  if (int(respCorrInit.size())!=corrArrSize) {
    std::cout << "doFit_v2: respCorrInit.size should be " << corrArrSize
	      << ". NUMTOWERS=" << NUMTOWERS << "\n";
    return NULL;
  }

  if (fDoCandTrackEnergyDiff<2) {
    std::cout << "doFit_v2: fDoCandTrackEnergyDiff expected to be >=2\n";
    std::cout << "  now it is " << fDoCandTrackEnergyDiff << "\n";
    return NULL;
  }

  // determine the balance resolution
  {
    TArrayD noRespCorr(corrArrSize);
    for (int i=0; i<corrArrSize; i++) noRespCorr[i]=1.;
    double sumWeight=0, sumB=0;
    for (unsigned int i=0; i<fData.size(); ++i) {
      double b=0, dummy1=0, dummy2=0;
      if (fDoCandTrackEnergyDiff==2) {
	GetBalance(fData[i],noRespCorr,b,dummy1);
      }
      else {
	GetMET(fData[i],noRespCorr,b,dummy1,dummy2);
      }
      double w= fData[i].GetWeight();
      sumWeight += w;
      sumB += b * w;
    }
    double bAvg= sumB/sumWeight;
    double sumDiffSqr=0;
    for (unsigned int i=0; i<fData.size(); ++i) {
      double b=0, dummy1=0, dummy2=0;
      if (fDoCandTrackEnergyDiff==2) {
	GetBalance(fData[i],noRespCorr,b,dummy1);
      }
      else {
	GetMET(fData[i],noRespCorr,b,dummy1,dummy2);
      }
      double w= fData[i].GetWeight();
      sumWeight += w;
      sumDiffSqr += pow(b-bAvg,2) * w;
    }
    fBalanceSigma = sqrt(sumDiffSqr/sumWeight);
    std::cout << "fBalanceSigma=" << fBalanceSigma << "\n";
  }

  // set the number of parameters to be the number of towers
  xMinuit=new TMinuit(corrArrSize);
  xMinuit->SetPrintLevel(fPrintLevel);
  //std::cout << "========================= " << xMinuit->GetMaxIterations() << std::endl;
  //xMinuit->SetMaxIterations(1000);
  xMinuit->SetErrorDef(0.5); // for a log likelihood
  xMinuit->SetFCN(FCN);
  xMinuit->SetObjectFit(this);

  // define the parameters
  for(int i=0; i<NUMTOWERS; i++) {
    int ieta=-MAXIETA+i;
    std::ostringstream oss;
    oss << "Tower ieta: " << ieta;
    xMinuit->DefineParameter(i, oss.str().c_str(), respCorrInit[i],
			     fParStep, fParMin, fParMax);
    if (std::find(fixTowers.begin(),fixTowers.end(), ieta)!=fixTowers.end()) {
      std::cout << "-- fix parameter\n";
      xMinuit->FixParameter(i);
    }
  }
  if (fDoCandTrackEnergyDiff==4) {
    std::ostringstream oss;
    oss << "HadEn scale: ";
    int idx=NUMTOWERS;
    xMinuit->DefineParameter(idx, oss.str().c_str(),
			     respCorrInit[idx], 0.1, 1e-3, 100);
  }
  //for (unsigned int i=0; i<fixTowers.size(); ++i) {
  //  xMinuit->FixParameter(fixTowers[i]);
  //}

  // Minimize
  xMinuit->mnscan();
  xMinuit->Migrad();

  // return the results
  TH1D* histo=new TH1D(hname,htitle,NUMTOWERS,-MAXIETA-0.5,MAXIETA+0.5);
  for(int i=1; i<=NUMTOWERS; i++) {
    Double_t val, error;
    xMinuit->GetParameter(i-1, val, error);
    histo->SetBinContent(i, val);
    histo->SetBinError(i, error);
  }

  for (unsigned int i=0; i<respCorrInit.size(); ++i) {
    Double_t val,error;
    xMinuit->GetParameter(i+1, val,error);
    respCorrInit[i]=val;
    if (int(i)>=NUMTOWERS) {
      std::cout << "spec param #" << i << " " << val << " +- " << error <<"\n";
    }
  }

  return histo;
}


void DijetRespCorrData::GetBalance(const DijetRespCorrDatum& datum, const TArrayD& respcorr, Double_t& balance_, Double_t& resolution_) const
{
  Double_t te, th, thf;
  Double_t pe, ph, phf;
  datum.GetTagEnergies(respcorr, te, th, thf);
  datum.GetProbeEnergies(respcorr, pe, ph, phf);
  
  // calculate the resolution and balance in E_T, not E
  Double_t tet=(te+th+thf)/std::cosh(datum.GetTagEta());
  Double_t pet=(pe+ph+phf)/std::cosh(datum.GetProbeEta());
  
  // correct the tag/probe E_T's for the "third jet"
  Double_t tpx = tet*std::cos(datum.GetTagPhi());
  Double_t tpy = tet*std::sin(datum.GetTagPhi());
  Double_t ppx = pet*std::cos(datum.GetProbePhi());
  Double_t ppy = pet*std::sin(datum.GetProbePhi());

  tpx += 0.5*datum.GetThirdJetPx();
  tpy += 0.5*datum.GetThirdJetPy();
  ppx += 0.5*datum.GetThirdJetPx();
  ppy += 0.5*datum.GetThirdJetPy();
  
  Double_t tetcorr = std::sqrt(tpx*tpx + tpy*tpy);
  Double_t petcorr = std::sqrt(ppx*ppx + ppy*ppy);

  balance_ = 2*(tetcorr-petcorr)/(tetcorr+petcorr);
  resolution_ = 0.035/datum.GetWeight(); //0.035
  return;
}

void DijetRespCorrData::GetMET(const DijetRespCorrDatum& datum,
			       const TArrayD& respcorr,
			       Double_t& met_,
			       Double_t &metDivAveEt_,
			       Double_t& resolution_) const
{
  Double_t te, th, thf;
  Double_t pe, ph, phf;
  datum.GetTagEnergies(respcorr, te, th, thf);
  datum.GetProbeEnergies(respcorr, pe, ph, phf);

  if (respcorr.GetSize()>NUMTOWERS) {
    th*= respcorr[NUMTOWERS];
    ph*= respcorr[NUMTOWERS];
    int idx=(respcorr.GetSize()==NUMTOWERS+1) ? NUMTOWERS : (NUMTOWERS+1);
    thf*=respcorr[idx];
    phf*=respcorr[idx];
    if (0) {
      double x=respcorr[NUMTOWERS];
      if (x!=1)
	std::cout << " last correction=" << x << std::endl;
    }
  }

  // calculate the resolution and balance in E_T, not E
  Double_t tet=(te+th+thf)/std::cosh(datum.GetTagEta());
  Double_t pet=(pe+ph+phf)/std::cosh(datum.GetProbeEta());

  // correct the tag/probe E_T's for the "third jet"
  Double_t tpx = tet*std::cos(datum.GetTagPhi());
  Double_t tpy = tet*std::sin(datum.GetTagPhi());
  Double_t ppx = pet*std::cos(datum.GetProbePhi());
  Double_t ppy = pet*std::sin(datum.GetProbePhi());

  Double_t sumPx= tpx + ppx + datum.GetThirdJetPx();
  Double_t sumPy= tpy + ppy + datum.GetThirdJetPy();

  met_= sqrt(sumPx*sumPx + sumPy*sumPy);
  metDivAveEt_ = met_*2/(tet+pet);

  resolution_ = 0.035/datum.GetWeight(); //0.035
  return;
}

void DijetRespCorrData::GetTrackEnergyDiff(const DijetRespCorrDatum& datum, const Int_t index, const TArrayD& respcorr, Double_t& Ediff_, Double_t& dEdiff_) const
{
  Double_t TrackP, EcalE, HcalE;
  datum.GetTrackVariables(respcorr, index, TrackP, EcalE, HcalE);
  
  Double_t TrackE = sqrt(TrackP*TrackP + 0.019479835145232396); // sqrt(p^2 + m_pion^2)

  Ediff_ = fabs(TrackE - EcalE - HcalE);
  dEdiff_ = 1.0;
  return;
}

void DijetRespCorrData::FCN(Int_t &npar, Double_t*, Double_t &f, Double_t *par, Int_t)
{
  if (0 && (npar==0)) std::cout << "npar=0\n"; // get rid of compiler msg
  // get the relevant data
  const DijetRespCorrData* data=dynamic_cast<const DijetRespCorrData*>(xMinuit->GetObjectFit());
  TArrayD respcorr;
  int arrSize=NUMTOWERS;
  if (data->GetDoCandTrackEnergyDiff()==4) arrSize++;
  respcorr.Set(arrSize, par);
  f = data->GetLikelihoodDistance(respcorr);
  
  return;
}
