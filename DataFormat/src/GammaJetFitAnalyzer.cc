
#ifdef __localRun
#  include "../interface/GammaJetFitData.h"
#  include "../interface/GammaJetFitAnalyzer.h"
#else
#  include "HcalClosureTest/DataFormat/interface/GammaJetFitData.h"
#  include "HcalClosureTest/DataFormat/interface/GammaJetFitAnalyzer.h"
#endif

#include <TPaveStats.h>

// ----------------------------------------------------

TH2D* GammaJetFitAnalyzer_t::plot_EtVsEt(const char *hName, const char *hTitle,
	const TArrayD *hcalCorrCf,
	int nBins, double EtMin, double EtMax) const
{
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

TH1D* GammaJetFitAnalyzer_t::plot_Energy
   (const char *hName, const char *hTitle,
    int tag, int transverse,
    const TArrayD *hcalCorrCf,
    int nBins, double EtMin, double EtMax) const
{
  TString xaxis=Form("E_{%s%s}",
		     (transverse) ? "T," : "",
		     (tag) ? "tag" : "probe");
  TString hTitleMdf= hTitle +TString(Form(";%s;count",xaxis.Data()));
  TH1D* h1= new TH1D(hName,hTitleMdf, nBins,EtMin,EtMax);
  //h1->SetStats("rmne");
  h1->SetDirectory(0);

  for (unsigned int i=0; i<fData->size(); i++) {
    const GammaJetEvent_t *e= fData->at(i);
    double w= e->GetWeight();
    double energy=0.;
    if (tag) {
      if (hcalCorrCf) energy= e->GetTagEtot(*hcalCorrCf);
      else energy= e->GetTagEtot();
      if (transverse) energy /= cosh(e->GetTagEta());
    }
    else {
      if (hcalCorrCf) energy= e->GetProbeEtot(*hcalCorrCf);
      else energy= e->GetProbeEtot();
      if (transverse) energy /= cosh(e->GetProbeEta());
    }
    h1->Fill(energy, w);
  }
  return h1;
}

// ----------------------------------------------------

TH1D* GammaJetFitAnalyzer_t::plot_EnergyOverGenEnergy
   (const char *hName, const char *hTitle,
    int tag, int transverse,
    const TArrayD *hcalCorrCf,
    int nBins, double RatioMin, double RatioMax) const
{
  if (transverse!=-1) {
    std::cout << "plot_EnergyOverGenEnergy: "
	      << "the code is not ready for transverse" << std::endl;
    return NULL;
  }

  TString xaxis=Form("E_{%s}/E_{gen}",
		     (tag) ? "tag" : "probe");
  TString hTitleMdf= hTitle +TString(Form(";%s;count",xaxis.Data()));
  TH1D* h1= new TH1D(hName,hTitleMdf, nBins,RatioMin,RatioMax);
  //h1->SetStats(0);
  h1->SetDirectory(0);

  for (unsigned int i=0; i<fData->size(); i++) {
    const GammaJetEvent_t *e= fData->at(i);
    double w= e->GetWeight();
    double energy=0., genEnergy=0.;
    if (tag) {
      if (hcalCorrCf) energy= e->GetTagEtot(*hcalCorrCf);
      else energy= e->GetTagEtot();
      genEnergy= e->GetAuxInfo().GetTagGenE();
      //if (transverse) energy /= cosh(e->GetTagEta());
    }
    else {
      if (hcalCorrCf) energy= e->GetProbeEtot(*hcalCorrCf);
      else energy= e->GetProbeEtot();
      genEnergy= e->GetAuxInfo().GetProbeGenE();
      //if (transverse) energy /= cosh(e->GetProbeEta());
    }
    if (genEnergy==double(0)) { energy=0; genEnergy=1; }
    h1->Fill(energy/genEnergy, w);
  }
  return h1;
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

TH1D* GammaJetFitAnalyzer_t::plot_TowerEntryCount(
                const char *hName, const char *hTitle,
		int weighted,
		double minWeight, int fractionFromMax,
		std::vector<Double_t> *towers_ptr) const
{
  std::vector<Double_t> towerWeight;
  std::vector<Int_t> emptyTowers;
  if (towers_ptr) towers_ptr->clear();

  if (!GetEmptyTowers(*fData,emptyTowers,weighted,
		      minWeight,fractionFromMax,&towerWeight)) {
    std::cout << "GammaJetFitAnalyzer_t::plot_TowerEntryCount: "
	      << "failed to get weights\n";
    return NULL;
  }

  TH1D *h1= new TH1D(hName,hTitle, NUMTOWERS+1, -MAXIETA-1, MAXIETA+1);
  h1->SetDirectory(0);
  h1->SetStats(0);

  //std::cout << "towerWeight.size=" << towerWeight.size() << "\n";
  for (unsigned int i=0; i<towerWeight.size(); ++i) {
    int idx= i-MAXIETA;
    //int idx_mirror= MAXIETA-i;
    //int i_mirror= NUMTOWERS-i-1;
    //std::cout << "i=" << i << ", idx=" << idx << ", w=" << towerWeight[i] << "\n";
    //std::cout << ", chkSymmetry idx_mirror=" << idx_mirror << ", i_mirror="
    //	      << i_mirror << "   "
    //	      << towerWeight[i] << "\n";
    h1->Fill( idx + 0.1 , towerWeight[i]);
  }

  if (towers_ptr) *towers_ptr= towerWeight;

  return h1;
}

// ----------------------------------------------------
// ----------------------------------------------------

void plotTwoHistos(TCanvas *cx, TString title,
		   TH1D* h1, TString label1, int color1, TString drawOpt1,
		   TH1D* h2, TString label2, int color2, TString drawOpt2)
{
  cx->cd();

  h1->SetLineColor(color1); h1->SetMarkerColor(color1);
  TH1D *h1Tmp= (TH1D*)h1->Clone(h1->GetName() + TString("_tmp"));
  h1Tmp->SetDirectory(0);
  h1Tmp->SetStats(0);
  h1Tmp->SetTitle(title);
  h1Tmp->Draw(drawOpt1);
  h1->Draw(drawOpt1 + TString("sames"));
  cx->Update();

  TPaveStats *stats1= (TPaveStats*)cx->GetPrimitive("stats");
  if (stats1) {
    stats1->SetName(h1->GetName() + TString("_stats"));
    stats1->SetY1NDC(.7);
    stats1->SetY2NDC(.85);
    stats1->SetTextColor(color1);
  }
  else {
    std::cout << "\tplotTwoHistos: stats are not available\n";
  }

  h2->SetLineColor(color2); h2->SetMarkerColor(color2);
  h2->Draw(drawOpt2 + TString("sames"));
  cx->Update();
  TPaveStats *stats2= (TPaveStats*)cx->GetPrimitive("stats");
  if (stats2) {
    stats2->SetName(h2->GetName() + TString("_stats"));
    stats2->SetY1NDC(.5);
    stats2->SetY2NDC(.65);
    stats2->SetTextColor(color2);
  }
  else {
    std::cout << "\tplotTwoHistos: stats are not available\n";
  }

  cx->Modified();
  cx->Update();
}



// ----------------------------------------------------

// ----------------------------------------------------
