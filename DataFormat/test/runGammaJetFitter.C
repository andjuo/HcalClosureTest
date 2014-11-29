#include "link_GammaJetFit.h"
#include "link_GammaJetFitAnalyzer.h"
#include <TCanvas.h>
#include <TLine.h>
#include "../interface/HistoCollector.h"
#include "helper.h"


void runGammaJetFitter(const TString fname="gjet_toy1_model2.root",
		       Long64_t maxEntries=-1,
		       TString outFileNameTag="",
		       double minFraction_user=-1.) {

  gBenchmark->Start("runGammaJetFitter");

  GammaJetCuts_t cuts;
  GammaJetFitter_t fitter;
  HistoCollector_t collector;

  //fitter.SetFittingProcedure(1); // 0 - gammaJet, 1 - dijet

  int nBins=50;
  double etMax=500;

  if (!LoadGammaJetEvents(fname,cuts,fitter,maxEntries)) return;
  HERE("load ok");

  std::vector<Int_t> fixTowers;
  double minFraction=0.1;
  //minFraction=0.00001;
  if (minFraction_user>=0) minFraction=minFraction_user;
  if (!GetEmptyTowers(fitter,fixTowers, 1, minFraction,1)) return;

  printVec("fixTowers",fixTowers);


  if (1) {
    GammaJetFitAnalyzer_t anObj(&fitter);
    std::vector<Double_t> weights;
    TH1D *hTowCount=
      anObj.plot_TowerEntryCount("hTowerCount",
				 "Events in towers;iEta_{tower};count",1,
				 0.,1,&weights);
    TCanvas *cTowCount= new TCanvas("cTowCount","cTowCount",600,600);
    hTowCount->Draw("LPE");
    if (1) {
      cTowCount->Update();
      std::cout << "cTowCount->GetUxmin()= " << cTowCount->GetUxmin() << ".. "<< cTowCount->GetUxmax() << "\n";
      double maxWeight= *std::max_element( weights.begin(), weights.end() );
      double lineAt= minFraction * maxWeight;
      TLine *line= new TLine(cTowCount->GetUxmin(),lineAt,
			     cTowCount->GetUxmax(),lineAt);
      line->Draw();
    }
    cTowCount->Update();

    collector.Add(hTowCount,cTowCount);
    //return;
  }

  std::vector<Double_t> hcalCorrCf(NUMTOWERS,1.);
  TArrayD cfAt1= convert(hcalCorrCf);
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

  if (0) { // plot jet energy uncorr vs corr
    /// This is a debug command set
    TArrayD cf0= convert(hcalCorrCf);
    for (int transverse=0; transverse<1; transverse++) {
      int tag=0; // jet is the probe
      GammaJetFitAnalyzer_t anObj(&fitter);
      TString hName1=(!transverse) ? "h1JetEnergy_Uncorr" : "h1JetET_Uncorr";
      TString hName2=(!transverse) ? "h1JetEnergy_Corr" : "h1JetET_Corr";
      TH1D* h1_uncorr = anObj.plot_Energy(hName1,
		       "Jet energy (uncorr)",tag,transverse,NULL,
					  nBins,0,etMax);
      TH1D *h1_corr = anObj.plot_Energy(hName2,
				      "Jet energy (corr)",tag,transverse,&cf0,
				      nBins,0,etMax);

      TString canvName = (transverse) ? "cJetEt" : "cJetEn";
      TString canvTitle= (transverse) ? "Jet transverse energy" : "Jet energy";
      TCanvas *cx= new TCanvas(canvName,canvName,600,600);
      plotTwoHistos(cx,canvTitle,
		    h1_uncorr,"uncorr",kBlack,"hist",
		    h1_corr,"corr",kBlue,"LP");
      cx->Update();

      collector.Add(h1_uncorr,h1_corr,cx);
    }
    return;
  }

  if (0) { // plot jet energy uncorr vs corr
    // for debug purposes
    TArrayD cf0= convert(hcalCorrCf);
    int tag=0; // jet is the probe
    GammaJetFitAnalyzer_t anObj(&fitter);
    TString hName1="h1JetEnergyRes_Uncorr";
    TString hName2="h1JetEnergyRes_Corr";
    TH1D* h1_uncorr = anObj.plot_EnergyOverGenEnergy(hName1,
		       "Jet energy resolution (uncorr)",tag,-1,NULL,
					  nBins,0,2.);
    TH1D *h1_corr = anObj.plot_EnergyOverGenEnergy(hName2,
		       "Jet energy resolution (corr)",tag,-1,&cf0,
					  nBins,0,2.);

    TString canvName = "cJetEnResol";
    TString canvFigTitle= "Jet energy resolution";
    TCanvas *cx= new TCanvas(canvName,canvName,600,600);
    plotTwoHistos(cx,canvFigTitle,
		  h1_uncorr,"uncorr",kBlack,"hist",
		  h1_corr,"corr",kBlue,"LP");
    cx->Update();

    collector.Add(h1_uncorr,h1_corr,cx);
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
    collector.Add(h2,cx);
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
    collector.Add(h2Norm,cy);
    //return;
  }

  //fixTowers.push_back(31);

  TH1D* h= fitter.doFit("hCorr","hCorr", hcalCorrCf, &fixTowers);
  if (!h) std::cout << "fitting failed\n";
  TArrayD cf= convert(hcalCorrCf);

  //printVec("hcalCorrCf ",hcalCorrCf);

  fitter.PrintFitValueChange(hcalCorrCf);

  if (1 && h) {
    TCanvas *cx= new TCanvas("cx","cx", 600,600);
    TH1D *hCoef= (TH1D*)h->Clone("hCoef");
    hCoef->SetStats(0);
    hCoef->SetDirectory(0);
    hCoef->SetTitle("correction coeffs;iEta_{tower};cf");
    hCoef->Draw("hist");
    cx->Update();

    collector.Add(hCoef,cx);
  }

  if (1) {
    std::cout << "hcalCorrCf " << hcalCorrCf << "\n";
    std::cout << "cf " << cf << "\n";

    GammaJetFitAnalyzer_t analyzer(&fitter);

    TH2D* h2EtVsEt_unCorr=
      analyzer.plot_EtVsEt("h2EtVsEtUncorr", "Et vs Et : uncorr",
			   NULL, nBins,0.,etMax);
    TH2D* h2EtVsEt_corr=
      analyzer.plot_EtVsEt("h2EtVsEtCorr","Et vs Et : corr",
			   &cf, nBins,0.,etMax);
    TCanvas *cEt0=new TCanvas("cEt0","cEt0",600,600);
    h2EtVsEt_unCorr->Draw("COLZ");
    cEt0->Update();

    TCanvas *cEt1=new TCanvas("cEt1","cEt1",600,600);
    h2EtVsEt_corr->Draw("COLZ");
    cEt1->Update();

    collector.Add(h2EtVsEt_unCorr, cEt0);
    collector.Add(h2EtVsEt_corr, cEt1);
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

    collector.Add(h2TowerEn_unCorr,cEt0);
    collector.Add(h2TowerEn_corr,cEt1);
  }

  if (0) {
    GammaJetFitAnalyzer_t anObj(&fitter);
    TH2D* h2=anObj.plot_TowerFitProfile("h2TowFittedProfile",
					"tower fit profile",
					0, 20, 0., 2., &hcalCorrCf);
    TCanvas *cx=new TCanvas("cxTFPfitted","cxTFPfitted",600,600);
    h2->Draw("COLZ");
    cx->Update();

    collector.Add(h2,cx);
    //return;
  }

  if (1) { // plot jet energy uncorr vs corr
    for (int transverse=0; transverse<2; transverse++) {
      int tag=0; // jet is the probe
      GammaJetFitAnalyzer_t anObj(&fitter);
      TString hName1=(!transverse) ? "h1JetEnergy_Uncorr" : "h1JetET_Uncorr";
      TString hName2=(!transverse) ? "h1JetEnergy_Corr" : "h1JetET_Corr";
      TH1D* h1_uncorr = anObj.plot_Energy(hName1,
		       "Jet energy (uncorr)",tag,transverse,NULL,
					  nBins,0,etMax);
      TH1D *h1_corr = anObj.plot_Energy(hName2,
				      "Jet energy (corr)",tag,transverse,&cf,
				      nBins,0,etMax);

      TString canvName = (transverse) ? "cJetEt" : "cJetEn";
      TString canvFigTitle=
	(transverse) ? "Jet transverse energy" : "Jet energy";
      TCanvas *cx= new TCanvas(canvName,canvName,600,600);
      plotTwoHistos(cx,canvFigTitle,
		    h1_uncorr,"uncorr",kBlack,"hist",
		    h1_corr,"corr",kBlue,"LP");
      cx->Update();

      collector.Add(h1_uncorr,h1_corr,cx);
    }
  }

  if (1) { // plot jet energy uncorr vs corr
    int tag=0; // jet is the probe
    GammaJetFitAnalyzer_t anObj(&fitter);
    TString hName1="h1JetEnergyRes_Uncorr";
    TString hName2="h1JetEnergyRes_Corr";
    TH1D* h1_uncorr = anObj.plot_EnergyOverGenEnergy(hName1,
						     "Jet energy resolution (uncorr)",tag,-1,//
						     //NULL,
						     &cfAt1,
					  nBins,0,2.);
    TH1D *h1_corr = anObj.plot_EnergyOverGenEnergy(hName2,
		       "Jet energy resolution (corr)",tag,-1,&cf,
					  nBins,0,2.);

    TString canvName = "cJetEnResol";
    TString canvFigTitle= "Jet energy resolution";
    TCanvas *cx= new TCanvas(canvName,canvName,600,600);
    plotTwoHistos(cx,canvFigTitle,
		  h1_uncorr,"uncorr",kBlack,"hist",
		  h1_corr,"corr",kBlue,"LP");
    cx->Update();

    collector.Add(h1_uncorr,h1_corr,cx);
  }

  if (1) {
    collector.Add(packMessages(3,"producedBy runGammaJetFitter",
			       fname.Data(),
			       Form("maxEntries=%ld",long(maxEntries))));
    if (!fitter.SaveInfoToFile(outFileNameTag,collector)) {
      std::cout << "error saving info\n";
    }
  }

  if (1) {
    collector.SaveCanvases(outFileNameTag);
    /*
    TString destDir="plots_" + outFileNameTag;
    for (unsigned int i=0; i<saveCanvas_vec.size(); ++i) {
      TString figName="fig-" + outFileNameTag + TString("-") +
	TString(saveCanvas_vec[i]->GetName());
      SaveCanvas(saveCanvas_vec[i],figName,destDir);
    }
    */
  }

  if (fitter.GetFittingProcedure()!=0) {
    std::cout << "special fitting procedure: " << fitter.GetFittingProcedure() << "\n";
  }

  gBenchmark->Show("runGammaJetFitter");
}
