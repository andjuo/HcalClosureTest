//#include "../interface/GammaJetFitData.h"
#include "pf_gammajettree.h"
#include <TH2D.h>
#include <TCanvas.h>
#include <TLine.h>
//#include "colorPalettes.hh"
#include "ComparisonPlot.hh"
#include "../interface/HistoCollector.h"
#include "helper.h"

// ------------------------------------------------------

#ifdef helper_HH
  TString plotOutDir="plots";
#endif

// ------------------------------------------------------

inline void prepareHisto(TH1D *h, int setStats=0) {
  h->SetDirectory(0);
  h->SetStats(setStats);
}

inline void prepareHisto(TH2D *h, int setStats=0) {
  h->SetDirectory(0);
  h->SetStats(setStats);
}

// ------------------------------------------------------
// ------------------------------------------------------

HistoCollector_t collector;
int saveCollection=0;

// ----------------


void displayHisto(int show, TH1D* h, TString tag, TString drawOpt="LPE");
void displayHisto(int show,
		  TH2D* h, TString tag, TString drawOpt="COLZ", int drawXeqY=1);

void displayHistoRatio(TH1D *h1, TString label1, TH1D *h2, TString label2,
		       TString xAxisLabel, TString yAxisLabel,
		       TString tag,
		       TString drawOpt1="hist", TString drawOpt2="hist");

// ------------------------------------------------------
// ------------------------------------------------------

void analyzeQCD_PFGammaJetTree(TString inpFileName,
			       int flatQCD=0,
			       Long64_t maxEntries=-1,
			       double extraWeightFactor=-1.,
			       TString plotOutDir_user="",
			       int saveCollection_user=-1) {

  std::cout << "inpFileName=<" << inpFileName << ">\n";

#ifdef helper_HH
  if (plotOutDir_user.Length()>0) plotOutDir= plotOutDir_user;
#endif
  if (saveCollection_user>=0) saveCollection=saveCollection_user;

  double flat_ptMin=0, flat_ptMax=0;
  if (flatQCD==1) { flat_ptMin=80; flat_ptMax=120; }
  else if (flatQCD==2) { flat_ptMin=120; flat_ptMax=170; }
  else if (flatQCD==3) { flat_ptMin=170; flat_ptMax=300; }
  if (flatQCD && (extraWeightFactor<0)) {
    std::cout << "flatQCD requires extraWeightFactor\n";
    return;
  }

  collector.Clear();

  // book histograms
  //double cPI= 4*atan(1);
  int show_evweight=1;
  TH1D *h1_evWeight= new TH1D("h1_evWeight","Event weight;w;count",1000,0,10);
  TH1D *h1_evPtHat= new TH1D("h1_evPtHat","Event #hat p_{T};#hat p_{T};count",100,0,500);
  TH1D *h1_evPtHatEvWeighted= new TH1D("h1_evPtHatEvWeighted","Event #hat p_{T};#hat p_{T};gen-weighted count",100,0,500);
  TH1D *h1_evPtHatWeighted= new TH1D("h1_evPtHatWeighted","Event #hat p_{T};#hat p_{T};weighted count",100,0,500);
  prepareHisto(h1_evWeight);
  prepareHisto(h1_evPtHat);
  prepareHisto(h1_evPtHatEvWeighted);
  prepareHisto(h1_evPtHatWeighted);

  double c_pt_max=500;
  double c_energy_max=1500;


  int show_jet_energy=1;
  TH1D* h1_jet_energyGen= new TH1D("h1_jet_energyGen",
				   "jet energy Gen;E_{gen};count",
				   100,0.,c_energy_max);
  TH1D* h1_jet_energyPF = new TH1D("h1_jet_energyPF",
				   "jet energy PF;E_{PF};count",
				   100,0.,c_energy_max);
  prepareHisto(h1_jet_energyGen,1);
  prepareHisto(h1_jet_energyPF,1);

  int show_jet_pT=1;
  TH1D* h1_jet_pTGen= new TH1D("h1_jet_pTGen",
				   "jet pT Gen;p_{T,gen};count",
				   100,0.,c_pt_max);
  TH1D* h1_jet_pTPF = new TH1D("h1_jet_pTPF",
				   "jet pT PF;p_{T,PF};count",
				   100,0.,c_pt_max);
  prepareHisto(h1_jet_pTGen,1);
  prepareHisto(h1_jet_pTPF,1);

  int show_jet_pT_window=(1 && flatQCD) ? 1:0;
  TH1D* h1_jet_pTGen_window= new TH1D("h1_jet_pTGen_window",
	       Form("jet pT Gen window %2.0lf-%2.0lf;p_{T,gen};count",
		   flat_ptMin,flat_ptMax),
				   100,0.,c_pt_max);
  TH1D* h1_jet_pTPF_window= new TH1D("h1_jet_pTPF_window",
	       Form("jet pT PFreco window %2.0lf-%2.0lf;p_{T,PF};count",
		   flat_ptMin,flat_ptMax),
				   100,0.,c_pt_max);
  prepareHisto(h1_jet_pTGen_window,1);
  prepareHisto(h1_jet_pTPF_window,1);

  int show_jet_energy_PF_over_Gen=1;
  TH1D *h1_jet_energy_PFoverGen= new TH1D("h1_jet_energy_PFoverGen",
		       "Jet energy ratio (PF div gen);E_{PF}/E_{gen};count",
					  100,0.,2.5);
  prepareHisto(h1_jet_energy_PFoverGen,1);

  int show_jet_pT_PF_over_Gen=1;
  TH1D *h1_jet_pT_PFoverGen= new TH1D("h1_jet_pT_PFoverGen",
	       "Jet p_{T} ratio (PF div gen);p_{T}^{PF}/p_{T}^{gen};count",
					  100,0.,2.5);
  prepareHisto(h1_jet_pT_PFoverGen,1);

  int show_jet2D_ptPF=1;
  TH2D *h2_jet_ptPF= new TH2D("h2_jet_ptPF",
			      "Leading jet p_{T};gen p_{T};reco p_{T}^{PF}",
			      100,0.,c_pt_max,
			      100,0.,c_pt_max);
  prepareHisto(h2_jet_ptPF);

  int show_photon_pt=1;
  TH1D* h1_photon_pTreco = new TH1D("h1_photon_pTreco",
				    "photon pT;p_{T}^{#gamma};count",
				   100,0.,c_pt_max);
  prepareHisto(h1_photon_pTreco,1);



  // process the data
  pf_gammajettree inpData(inpFileName);

  inpData.DeactivateBranches();
  //inpData.ActivateBranches_forFitSkim();
  inpData.ActivateBranches(2,"EventWeight","EventPtHat");
  inpData.ActivateBranches(2,"ppfjet_pt","ppfjet_E");
  inpData.ActivateBranches(2,"ppfjet_genpt","ppfjet_genE");
  inpData.ActivateBranches(1,"tagPho_pt");
  //inpData.ActivateBranches_genBasicSet();

  // read in the file
  Long64_t nEntries= inpData.fChain->GetEntries();
  std::cout << "nEntries=" << nEntries << "\n";
  if (nEntries==0) return;
  if (maxEntries<0) maxEntries=nEntries;
  Long64_t nBytes=0, passedCount=0;
  Double_t sum_evWeight_window=0;
  Double_t sum_finalW_window=0;

  for (Long64_t iEntry=0; (iEntry<nEntries) && (iEntry<maxEntries);
       iEntry++) {
    if (inpData.LoadTree(iEntry) < 0) break;
    Long64_t nb = inpData.GetEntry(iEntry);
    nBytes += nb;
    if (iEntry%10000==0) std::cout << " ... reading entry " << iEntry << "\n";
    //std::cout << "ientry=" << iEntry << "\n";

    //if (!inpData.passCuts(1,2)) continue; // recommended use

    passedCount++;

    double w=inpData.EventWeight;
    if (0 && flatQCD) {
      // It seems this factor is already included in the gen.event weight
      //w /= pow(inpData.EventPtHat/double(15.), double(4.5));
    }
    if (extraWeightFactor>0) w*= extraWeightFactor;

    h1_evWeight->Fill( w, 1. );
    h1_evPtHat->Fill( inpData.EventPtHat, 1. );
    h1_evPtHatEvWeighted->Fill( inpData.EventPtHat, inpData.EventWeight );
    h1_evPtHatWeighted->Fill( inpData.EventPtHat, w );
    h1_jet_energyGen->Fill( inpData.ppfjet_genE, w);
    h1_jet_energyPF->Fill( inpData.ppfjet_E, w);
    h1_jet_energy_PFoverGen->Fill( inpData.ppfjet_E/inpData.ppfjet_genE, w);
    h1_jet_pTGen->Fill( inpData.ppfjet_genpt, w);
    h1_jet_pTPF->Fill( inpData.ppfjet_pt, w);
    if ((inpData.EventPtHat>flat_ptMin) &&
	(inpData.EventPtHat<flat_ptMax)) {
      sum_evWeight_window+= inpData.EventWeight;
      sum_finalW_window+= w;
      h1_jet_pTGen_window->Fill( inpData.ppfjet_genpt, w);
      h1_jet_pTPF_window->Fill( inpData.ppfjet_pt, w);
    }
    h1_jet_pT_PFoverGen->Fill( inpData.ppfjet_pt/inpData.ppfjet_genpt, w);
    h2_jet_ptPF->Fill( inpData.ppfjet_genpt, inpData.ppfjet_pt, w);

    h1_photon_pTreco->Fill( inpData.tagPho_pt, w );

    int gen_debug=0;
    if (gen_debug==1) {
      std::cout << "iEntry=" << iEntry << ",    leading jet genE=" << inpData.ppfjet_genE << ", pfE=" << inpData.ppfjet_E << "\n";

      /*
      std::cout << "  ppfjet_twr_ieta " << (inpData.ppfjet_twr_ieta) << "\n";
      std::cout << "  ppfjet_twr_hade " << (inpData.ppfjet_twr_hade) << "\n";
      std::cout << "  ppfjet_twr_frac " << (inpData.ppfjet_twr_frac) << "\n";
      std::cout << "  ppfjet_twr_clusterind " << (inpData.ppfjet_twr_clusterind) << "\n";
      std::cout << "  ppfjet_cluster_dR " << (inpData.ppfjet_cluster_dR) << "\n";
      */
    }
    else if (gen_debug==2) {
      std::cout << "genE=" << inpData.ppfjet_genE << ", pfE=" << inpData.ppfjet_E << "\n";
    }

  }

  std::cout << "nEntries=" << nEntries << ", passedCount=" << passedCount << "\n";

  displayHisto(show_evweight, h1_evWeight,"evWeight","LPE");
  displayHisto(show_evweight, h1_evPtHat, "evPtHat", "LPE");
  displayHisto(show_evweight, h1_evPtHatEvWeighted,"evPtHatEvWeighted", "LPE");
  displayHisto(show_evweight, h1_evPtHatWeighted, "evPtHatWeighted", "LPE");

  displayHisto(show_jet_energy, h1_jet_energyGen, "jet_E_gen","LPE");
  displayHisto(show_jet_energy, h1_jet_energyPF, "jet_E_recoPF","LPE");
  displayHisto(show_jet_energy_PF_over_Gen, h1_jet_energy_PFoverGen,"e_PFoverGen","hist");

  displayHisto(show_jet_pT, h1_jet_pTGen, "jet_pT_gen","LPE");
  displayHisto(show_jet_pT, h1_jet_pTPF, "jet_pT_recoPF","LPE");
  displayHisto(show_jet_pT_PF_over_Gen, h1_jet_pT_PFoverGen,"pT_PFoverGen","hist");
  //if (sum_window!=Double_t(0)) h1_jet_pTGen_window->Scale(1/sum_window);
  displayHisto(show_jet_pT_window, h1_jet_pTGen_window, "jet_pT_gen_window","LPE");
  displayHisto(show_jet_pT_window, h1_jet_pTPF_window, "jet_pT_reco_window","LPE");

  displayHisto(show_jet2D_ptPF, h2_jet_ptPF, "jet2D_ptPF","COLZ");

  displayHisto(show_photon_pt, h1_photon_pTreco, "photon_pTreco","LPE");

  if (saveCollection) {
    TString outFName=TString("saved_") + plotOutDir + TString(".root");
    collector.Add("producedBy_analyzePFGammaJetTree");
    collector.SaveToFile(outFName);
  }

  if (flatQCD) {
    std::cout << "flatQCD pt window : " << flat_ptMin << " .. "
	      << flat_ptMax << "\n";
    std::cout << "sum_window. evWeight=" << sum_evWeight_window
	      << ", finalW=" << sum_finalW_window << "\n";
  }

  return;
}

// ------------------------------------------------------

void displayHisto(int show, TH1D* h, TString tag, TString drawOpt)
{
  bool isBatch= gROOT->IsBatch();
  if (!show && !saveCollection) return;
  if (!show &&  saveCollection) gROOT->SetBatch(true);

  TString canvName="c_" + tag;
  TCanvas *c=new TCanvas(canvName,canvName,600,600);
  h->Draw(drawOpt);
  c->Update();
  TPaveStats *stats= (TPaveStats*)c->GetPrimitive("stats");
  if (stats) {
    stats->SetY1NDC(0.7);
    stats->SetY2NDC(0.85);
    c->Update();
  }
  if (show) {
    double effSigma= calc_effSigma(h);
    std::cout << "histogram " << h->GetName()<<", effSigma="<<effSigma<<"\n";
  }

#ifdef helper_HH
  TString figName="fig-" + tag;
  SaveCanvas(c,figName,plotOutDir);
#endif

  if (saveCollection) collector.Add(h,c);
  if (!show && saveCollection) gROOT->SetBatch(isBatch);
}

// ------------------------------------------------------

void displayHisto(int show,
		  TH2D* h, TString tag, TString drawOpt, int drawXeqY)
{
  bool isBatch= gROOT->IsBatch();
  if (!show && !saveCollection) return;
  if (!show &&  saveCollection) gROOT->SetBatch(true);

  TString canvName="c_" + tag;
  TCanvas *c=new TCanvas(canvName,canvName,600,600);
#ifdef ColorPalettes_HH
  AdjustFor2DplotWithHeight(c,0.18);
  SetSideSpaces(c,0.05,0.0,0.0,0.0);
  h->GetYaxis()->SetTitleOffset(1.8);
#endif
  h->Draw(drawOpt);
  if (drawXeqY) {
    c->Update();
    Double_t xmin= c->GetUxmin();
    Double_t xmax= c->GetUxmax();
    Double_t ymin= c->GetUymin();
    Double_t ymax= c->GetUymax();
    if (xmin>ymin) ymin=xmin; else xmin=ymin;
    if (xmax<ymax) xmax=ymax; else ymax=xmax;
    if (xmin<xmax) {
      TLine *line= new TLine(xmin,ymin,xmax,ymax);
      line->Draw();
    }
    else {
      std::cout << "\tcould not draw x=y line\n";
    }
  }
  c->Update();
  TPaveStats *stats= (TPaveStats*)c->GetPrimitive("stats");
  if (stats) {
    stats->SetY1NDC(0.7);
    stats->SetY2NDC(0.85);
    c->Update();
  }
#ifdef helper_HH
  TString figName="fig-" + tag;
  SaveCanvas(c,figName,plotOutDir);
#endif

  if (saveCollection) collector.Add(h,c);
  if (!show &&  saveCollection) gROOT->SetBatch(isBatch);
}

// ------------------------------------------------------

void displayHistoRatio(TH1D *h1, TString label1, TH1D *h2, TString label2,
		       TString xAxisLabel, TString yAxisLabel, TString tag,
		       TString drawOpt1, TString drawOpt2)
{
#ifndef ComparisonPlot_HH
  std::cout << "displayHistoRatio: ComparisonPlot.hh is not included\n";
  return;
#endif
#ifdef ComparisonPlot_HH
  ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"cp","",
		      xAxisLabel,yAxisLabel,"ratio");

  cp.AddHist1D(h1,label1,drawOpt1, kBlack,1,0,1);
  cp.AddHist1D(h2,label2,drawOpt2, kBlue, 1,0,1);
  TString canvName="c_" + tag;
  TCanvas *c= new TCanvas(canvName,canvName,600,700);
  cp.Prepare2Pads(c);
  cp.Draw(c);
  c->Update();
#endif
}


// ------------------------------------------------------
