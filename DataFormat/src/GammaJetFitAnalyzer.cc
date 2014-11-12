
#ifdef __localRun
#  include "../interface/GammaJetFitData.h"
#  include "../interface/GammaJetFitAnalyzer.h"
#else
#  include "HcalClosureTest/DataFormat/interface/GammaJetFitData.h"
#  include "HcalClosureTest/DataFormat/interface/GammaJetFitAnalyzer.h"
#endif


// ----------------------------------------------------

TH2D* GammaJetFitAnalyzer_t::plot_EtVsEt(const char *hName, const char *hTitle,
	const TArrayD *hcalCorrCf,
	int nBins, double EtMin, double EtMax) const {
  TString hTitleMdf= hTitle + TString(";E_{T,tag};E_{T,probe}");
  TH2D* h2= new TH2D(hName,hTitleMdf,
		     nBins,EtMin,EtMax,
		     nBins,EtMin,EtMax);
  h2->SetStats(0);
  h2->SetDirectory(0);

  for (unsigned int i=0; i<fData->size(); i++) {
    const GammaJetEvent_t *e= fData->at(i);
    double w= e->GetWeight();
    double eT_tag = (hcalCorrCf) ?
      e->GetTagETtot(*hcalCorrCf) : e->GetTagETtot();
    double eT_probe = (hcalCorrCf) ?
      e->GetProbeETtot(*hcalCorrCf) : e->GetProbeETtot();
    h2->Fill(eT_tag,eT_probe, w);
  }
  return h2;
}

// ----------------------------------------------------

TH2D* GammaJetFitAnalyzer_t::plot_TowerEn(const char *hNameBase,
					  const char *hTitle,
	 unsigned int idxMin, unsigned int idxMax,
	 int plotProbe,
	 const TArrayD *hcalCorrCf) const {
  TString tagStr=(plotProbe) ? "" : "_tag";
  TString hName= TString(hNameBase) +
    TString(Form("_%d_%d",int(idxMin),int(idxMax))) + tagStr;
  TString hTitleMdf= hTitle + TString(";idx;iEta");
  TH2D* h2= new TH2D(hName,hTitleMdf,
		     int(idxMax-idxMin+1),idxMin, idxMax+1,
		     NUMTOWERS+2, -MAXIETA-1, MAXIETA+1);
  h2->SetStats(0);
  h2->SetDirectory(0);

  for (unsigned int i=idxMin; (i<idxMax) && (i<fData->size()); i++) {
    const GammaJetEvent_t *e= fData->at(i);
    const std::map<Int_t,Double_t> *hMap=
      (plotProbe) ? &e->GetProbeHcalE() : &e->GetTagHcalE();
    for (std::map<Int_t,Double_t>::const_iterator it=hMap->begin();
	 it!=hMap->end(); it++) {
      int idx= it->first;
      double val= it->second;
      //std::cout << "idx=" << idx << ", val=" << val << "\n";
      if (hcalCorrCf) val *= (*hcalCorrCf)[idx+MAXIETA];
      h2->Fill(i,idx,val);
    }
  }
  return h2;
}

// ----------------------------------------------------

TH2D* GammaJetFitAnalyzer_t::plot_TowerFitProfile
    (const char *hName, const char *hTitle,
     int normalized,
     int nBins, double cfMin, double cfMax,
     const std::vector<double> *setCfs) const
{

  TString hTitleMdf= hTitle + TString(";Cf_{iEta};iEta");
  TH2D* h2= new TH2D(hName,hTitleMdf,
		     nBins,cfMin,cfMax,
		     NUMTOWERS+1, -MAXIETA-1, MAXIETA+1);
  h2->SetStats(0);
  h2->SetDirectory(0);

  TArrayD cfArr(NUMTOWERS);
  cfArr.Reset(1);
  if (setCfs) {
    for (unsigned int i=0; i<setCfs->size(); ++i) {
      cfArr[i] = setCfs->at(i);
    }
  }

  TAxis *ax= h2->GetXaxis();
  TAxis *ay= h2->GetYaxis();

  for (int jbin=1; jbin<=h2->GetNbinsY(); ++jbin) {
    int iEta= int(ay->GetBinCenter(jbin));
    if (iEta<-MAXIETA) continue;
    if (iEta> MAXIETA) continue;
    int idx= iEta+MAXIETA;
    //std::cout << "iEta=" << iEta << ", idx=" << idx << "\n";

    double minVal=1e9, maxVal=-1e9;
    for (int ibin=1; ibin<=h2->GetNbinsX(); ++ibin) {
      double cfVal= ax->GetBinCenter(ibin);
      double storeVal= cfArr[idx];
      cfArr[idx]=cfVal;
      //std::cout << "cfArr= " << cfArr << "\n";
      double fitVal= fData->GetFitValue(cfArr);
      if (fitVal<minVal) minVal=fitVal;
      if (fitVal>maxVal) maxVal=fitVal;
      h2->SetBinContent(ibin,jbin, fitVal);
      cfArr[idx]=storeVal;
    }
    //std::cout << "GammaJetFitAnalyzer_t::plot_TowerFitProfile: minVal="
    //	      << minVal << ", maxVal=" << maxVal << "\n";
    if (normalized) {
      double norm=fabs(maxVal);
      if (fabs(minVal)>norm) norm=fabs(minVal);
      if (norm<1e-6) {
	std::cout << "WARNING: small norm=" << norm << ". Reset to 1\n";
	norm=1.;
      }
      for (int ibin=1; ibin<=h2->GetNbinsX(); ++ibin) {
	double x= h2->GetBinContent(ibin,jbin);
	h2->SetBinContent(ibin,jbin,x/norm);
      }
    }
  }

  return h2;
}

// ----------------------------------------------------
