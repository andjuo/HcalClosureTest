#include "link_GammaJetFit.h"
#include "link_GammaJetFitAnalyzer.h"
#include <TCanvas.h>

void runGammaJetFitter(const TString fname="gjet_toy1_model2.root",
		       Long64_t maxEntries=-1) {

  GammaJetCuts_t cuts;
  GammaJetFitter_t fitter;

  if (!LoadGammaJetEvents(fname,cuts,fitter,maxEntries)) return;
  HERE("load ok");

  std::vector<Int_t> fixTowers;
  if (!GetEmptyTowers(fitter,fixTowers, 1, 0.1,1)) return;

  printVec("fixTowers",fixTowers);

  if (1) {
    GammaJetFitAnalyzer_t anObj(&fitter);
    TH1D *hTowCount=
      anObj.plot_TowerEntryCount("hTowerCount",
				 "Events in towers;iEta_{tow};count",
				 1,0.1,1);
    TCanvas *cTowCount= new TCanvas("cTowCount","cTowCount",600,600);
    hTowCount->Draw("LPE");
    cTowCount->Update();
    //return;
  }

  std::vector<Double_t> hcalCorrCf(NUMTOWERS,1.);
  if (0) {
    TArrayD cf0= convert(hcalCorrCf);
    std::cout << "init values " << hcalCorrCf << "\n";
    std::cout << " 1st " << fitter.at(0)->CalcDiffEt() << "\n";
    std::cout << " cmp " << fitter.at(0)->CalcDiffEt(cf0) << "\n";

    const GammaJetEvent_t *e = fitter.at(0);
    std::cout << " Et " << e->GetTagETtot() << " vs " << e->GetProbeETtot() << "\n";
    std::cout << " cmp " << e->GetTagETtot(cf0) << " vs " << e->GetProbeETtot(cf0) << "\n";
    return;
  }

  std::cout << "estimated resolution " << fitter.EstimateResolution(1) << "\n";

  if (0) {
    GammaJetFitAnalyzer_t anObj(&fitter);
    TH2D* h2=anObj.plot_TowerFitProfile("h2TowFitProfile","tower fit profile",
					0);
    TCanvas *cx=new TCanvas("cxTFP","cxTFP",600,600);
    h2->Draw("COLZ");
    cx->Update();
    //return;
  }
  if (0) {
    GammaJetFitAnalyzer_t anObj(&fitter);
    TH2D* h2Norm=anObj.plot_TowerFitProfile("h2TowFitProfile",
					    "tower fit profile (normalized)",
					    1);
    TCanvas *cy=new TCanvas("cyTFP","cyTFP",600,600);
    h2Norm->Draw("COLZ");
    cy->Update();
    //return;
  }

  //fixTowers.push_back(31);

  TH1D* h= fitter.doFit("hCorr","hCorr", hcalCorrCf, &fixTowers);
  TArrayD cf= convert(hcalCorrCf);

  //printVec("hcalCorrCf ",hcalCorrCf);

  fitter.PrintFitValueChange(hcalCorrCf);

  if (1 && h) {
    TCanvas *cx= new TCanvas("cx","cx", 600,600);
    h->Draw("hist");
    cx->Update();
  }

  if (1) {
    std::cout << "hcalCorrCf " << hcalCorrCf << "\n";
    std::cout << "cf " << cf << "\n";

    GammaJetFitAnalyzer_t analyzer(&fitter);

    TH2D* h2EtVsEt_unCorr=
      analyzer.plot_EtVsEt("h2EtVsEtUncorr", "Et vs Et : uncorr",
						NULL);
    TH2D* h2EtVsEt_corr=
      analyzer.plot_EtVsEt("h2EtVsEtCorr","Et vs Et : corr",
			   &cf);
    TCanvas *cEt0=new TCanvas("cEt0","cEt0",600,600);
    h2EtVsEt_unCorr->Draw("COLZ");
    cEt0->Update();

    TCanvas *cEt1=new TCanvas("cEt1","cEt1",600,600);
    h2EtVsEt_corr->Draw("COLZ");
    cEt1->Update();
  }

  if (0) {
    GammaJetFitAnalyzer_t analyzer(&fitter);

    TH2D* h2TowerEn_unCorr=
      analyzer.plot_TowerEn("h2TowerEn_uncorr","uncorr tower en",
			    0,20,1,
			    NULL);
    TH2D* h2TowerEn_corr=
      analyzer.plot_TowerEn("h2TowerEn_corr","corr tower en",
			    0,20,1,
			    &cf);

    TCanvas *cEt0=new TCanvas("cTowE0","cTowE0",600,600);
    h2TowerEn_unCorr->Draw("COLZ");
    cEt0->Update();

    TCanvas *cEt1=new TCanvas("cTowE1","cTowE1",600,600);
    h2TowerEn_corr->Draw("COLZ");
    cEt1->Update();
  }

  if (0) {
    GammaJetFitAnalyzer_t anObj(&fitter);
    TH2D* h2=anObj.plot_TowerFitProfile("h2TowFittedProfile",
					"tower fit profile",
					0, 20, 0., 2., &hcalCorrCf);
    TCanvas *cx=new TCanvas("cxTFPfitted","cxTFPfitted",600,600);
    h2->Draw("COLZ");
    cx->Update();
    //return;
  }

  if (1) { // plot jet energy uncorr vs corr
    int transverse=0;
    int tag=0; // jet is the probe
    int nBins=50;
    double etMax=200;
    GammaJetFitAnalyzer_t anObj(&fitter);
    TH1D* h1_uncorr = anObj.plot_Energy("h1JetEnergyUncorr",
			      "Jet energy (uncorr)",tag,transverse,NULL,
					nBins,0,etMax);
    TH1D *h1_corr = anObj.plot_Energy("h1JetEnergyCorr",
				      "Jet energy (corr)",tag,transverse,&cf,
				      nBins,0,etMax);
    TCanvas *cx= new TCanvas("cJetEn","cJetEn",600,600);
    plotTwoHistos(cx,"Jet energy",
		  h1_uncorr,"uncorr",kBlack,"hist",
		  h1_corr,"corr",kRed,"LP");
    //printHisto(h1_uncorr);
    //printHisto(h1_corr);
    //h1_uncorr->Draw("hist");
    //h1_corr->SetLineColor(kRed);
    //h1_corr->Draw("LPsame");
    cx->Update();
  }

}
