#include "../interface/GammaJetFitData.h"
#include "link_PFDijetTree.h"
#include <TH2D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TText.h>
//#include "colorPalettes.hh"
#include "ComparisonPlot.hh"
#include "../interface/HistoCollector.h"
#include "helper.h"

// ------------------------------------------------------

#ifdef helper_HH
  TString plotOutDir="plots";
#endif

// ------------------------------------------------------

inline void prepareHisto(TH1D *h) {
  h->SetDirectory(0);
  h->SetStats(0);
}

inline void prepareHisto(TH2D *h) {
  h->SetDirectory(0);
  h->SetStats(0);
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

void analyzePFDiJetTree(TString inpFileName,
			Long64_t maxEntries=-1,
			double extraWeightFactor=-1.,
			TString plotOutDir_user="") {

#ifdef helper_HH
  if (plotOutDir_user.Length()>0) plotOutDir= plotOutDir_user;
#endif

  // book histograms
  double cPI= 4*atan(1);
  int show_dPhi=0;
  TH1D *h1_dPhi= new TH1D("h1_dPhi","#Delta#Phi;#Delta#Phi;count",100,-2*cPI,2*cPI);
  prepareHisto(h1_dPhi);

  double c_pt_max=500;
  //double c_energy_max=1500;

  int show_dPt_PF=1;
  TH1D *h1_dPt_PF= new TH1D("h1_dPt_PF",
			    "#Deltap_{T}^{PF};#Deltap_{T}^{PF};count",
			    100,-c_pt_max,c_pt_max);
  prepareHisto(h1_dPt_PF);

  /*
  int show_jet_tightID=1;
  TH1D* h1_tjet_tightID= new TH1D("h1_tjet_tightID",
				  "Tag jet tightID flag;flag;count",
				  3,0.,2.);
  TH1D* h1_pjet_tightID= new TH1D("h1_pjet_tightID",
				  "Probe jet tightID flag;flag;count",
				  3,0.,2.);
  prepareHisto(h1_tjet_tightID);
  prepareHisto(h1_pjet_tightID);

  int show_jet_energy_PF_over_Gen_sideBySide=0;
  TH1D* h1_jet_energyGen= new TH1D("h1_jet_energyGen",
				   "jet energy Gen;E_{gen};count",
				   100,0.,c_energy_max);
  TH1D* h1_jet_energyPF = new TH1D("h1_jet_energyPF",
				   "jet energy PF;E_{PF};count",
				   100,0.,c_energy_max);
  TH1D* h1_jet_energyRH = new TH1D("h1_jet_energyRH",
				   "jet energy RH;E_{RecHits};count",
				   100,0.,c_energy_max);
  prepareHisto(h1_jet_energyGen);
  prepareHisto(h1_jet_energyPF);
  prepareHisto(h1_jet_energyRH);
  */


  int show_jet_energy_PF_over_Gen=1;
  TH1D *h1_tjet_energy_PFoverGen= new TH1D("h1_tjet_energy_PFoverGen",
	  "Tag jet energy ratio (PF div gen);E_{PF}^{j1}/E_{gen}^{j1};count",
					  100,0.,2.5);
  TH1D *h1_tjet_energy_RHoverGen= new TH1D("h1_tjet_energy_RHoverGen",
	  "Tag jet energy ratio (RH div gen);E_{RH}^{j1}/E_{gen}^{j1};count",
					  100,0.,2.5);
  TH1D *h1_pjet_energy_PFoverGen= new TH1D("h1_pjet_energy_PFoverGen",
	  "Probe jet energy ratio (PF div gen);E_{PF}^{j2}/E_{gen}^{j2};count",
					  100,0.,2.5);
  TH1D *h1_pjet_energy_RHoverGen= new TH1D("h1_pjet_energy_RHoverGen",
	  "Probe jet energy ratio (RH div gen);E_{RH}^{j2}/E_{gen}^{j2};count",
					  100,0.,2.5);
  prepareHisto(h1_tjet_energy_PFoverGen);
  prepareHisto(h1_tjet_energy_RHoverGen);
  h1_tjet_energy_PFoverGen->SetStats(1);
  h1_tjet_energy_RHoverGen->SetStats(1);
  prepareHisto(h1_pjet_energy_PFoverGen);
  prepareHisto(h1_pjet_energy_RHoverGen);
  h1_pjet_energy_PFoverGen->SetStats(1);
  h1_pjet_energy_RHoverGen->SetStats(1);

  /*
  int show_jet_pT_PF_over_Gen=1;
  TH1D *h1_jet_pT_PFoverGen= new TH1D("h1_jet_pT_PFoverGen",
	       "Jet E_{T} ratio (PF div gen);p_{T}^{PF}/p_{T}^{gen};count",
					  100,0.,2.5);
  TH1D *h1_jet_pT_RHoverGen= new TH1D("h1_jet_pT_RHoverGen",
	       "Jet E_{T} ratio (RH div gen);p_{T}^{RH}/p_{T}^{gen};count",
					  100,0.,2.5);
  prepareHisto(h1_jet_pT_PFoverGen);
  prepareHisto(h1_jet_pT_RHoverGen);
  h1_jet_pT_PFoverGen->SetStats(1);
  h1_jet_pT_RHoverGen->SetStats(1);
  */

  int show_jet2D_ptPF=1;
  TH2D *h2_tjet_ptPF= new TH2D("h2_tjet_ptPF",
			      "Leading jet p_{T};gen p_{T};reco p_{T}^{PF}",
			      100,0.,c_pt_max,
			      100,0.,c_pt_max);
  TH2D *h2_pjet_ptPF= new TH2D("h2_pjet_ptPF",
			      "Subleading jet p_{T};gen p_{T};reco p_{T}^{PF}",
			      100,0.,c_pt_max,
			      100,0.,c_pt_max);
  prepareHisto(h2_tjet_ptPF);
  prepareHisto(h2_pjet_ptPF);

  int show_jet2D_ptRH=1;
  TH2D *h2_tjet_ptRH= new TH2D("h2_tjet_ptRH",
		      "Leading jet p_{T};gen p_{T};reco p_{T}^{RecHits}",
			      100,0.,c_pt_max,
			      100,0.,c_pt_max);
  TH2D *h2_pjet_ptRH= new TH2D("h2_pjet_ptRH",
		     "Subleading jet p_{T};gen p_{T};reco p_{T}^{RecHits}",
			      100,0.,c_pt_max,
			      100,0.,c_pt_max);
  prepareHisto(h2_tjet_ptRH);
  prepareHisto(h2_pjet_ptRH);

  int show_jet_vs_jet_genPt=1;
  TH2D *h2_jet_vs_jet_genPt= new TH2D("h2_jet_vs_jet_genPt",
				   "Gen p_{T};jet1 gen p_{T};jet2 gen p_{T}",
				   100,0.,c_pt_max,
				   100,0.,c_pt_max);
  prepareHisto(h2_jet_vs_jet_genPt);

  /*
  int show_pho_vs_jet_recoPtPF=1;
  TH2D *h2_pho_vs_jet_recoPtPF= new TH2D("h2_pho_vs_jet_recoPtPF",
			   "Reco p_{T};#gamma p_{T};jet p_{T}^{PF}",
				   100,0.,c_pt_max,
				   100,0.,c_pt_max);
  prepareHisto(h2_pho_vs_jet_recoPtPF);

  int show_pho_vs_jet_recoPtRH=1;
  TH2D *h2_pho_vs_jet_recoPtRH= new TH2D("h2_pho_vs_jet_recoPtRH",
			   "Reco p_{T};#gamma p_{T};jet p_{T}^{RecHits}",
				   100,0.,c_pt_max,
				   100,0.,c_pt_max);
  prepareHisto(h2_pho_vs_jet_recoPtRH);
  */

  int show_jet_vs_jet_recoEta=1;
  TH2D *h2_jet_vs_jet_recoEta= new TH2D("h2_jet_vs_jet_recoEta",
					"Reco #eta;jet1 #eta;jet2 #eta",
					100,-2.5*cPI, 2.5*cPI,
					100,-2.5*cPI, 2.5*cPI);
  prepareHisto(h2_jet_vs_jet_recoEta);

  // process the data
  dijet_PFNtuple inpData(inpFileName);

  //inpData.DeactivateBranches();

  GammaJetEvent_t *dt= new GammaJetEvent_t();
  GammaJetEventAuxInfo_t aux;

  // read in the file
  Long64_t nEntries= inpData.fChain->GetEntries();
  if (maxEntries<0) maxEntries=nEntries;
  Long64_t nBytes=0, passedCount=0;

  for (Long64_t iEntry=0; (iEntry<nEntries) && (iEntry<maxEntries);
       iEntry++) {
    if (inpData.LoadTree(iEntry) < 0) break;
    Long64_t nb = inpData.GetEntry(iEntry);
    nBytes += nb;
    if (iEntry%1000==0) std::cout << " ... reading entry " << iEntry << "\n";
    //std::cout << "ientry=" << iEntry << "\n";

    //if (!inpData.passCuts(1,2)) continue; // recommended use

    passedCount++;

     // Note: the dijet event has the tag jet (idx=1) and the probe jet (0)

    double tagjet_ecalE= inpData.getSumEcalE(1,1);
    double tagjet_hcalE_noRecHits= inpData.getSumHcalE_trackDiffEcal(1);

    double ecalE= inpData.getSumEcalE(0,1);
    double hcalE_noRecHits= inpData.getSumHcalE_trackDiffEcal(0);

    aux.SetEventNo(-1); //inpData.EventNumber);
    aux.SetRunNo(-1); //inpData.RunNumber);
    aux.SetGenE(inpData.tpfjet_genE,inpData.ppfjet_genE);
    aux.SetProbeHcalENoRecHits(hcalE_noRecHits);
    dt->SetAuxInfo(aux);

    double w=inpData.pf_weight;
    if (extraWeightFactor>0) w*= extraWeightFactor;
    dt->SetWeight(w);

    dt->SetTagEtaPhiEn(inpData.tpfjet_eta, inpData.tpfjet_phi,
		       tagjet_ecalE + tagjet_hcalE_noRecHits,
		       inpData.getHcalEMap(1,1e-4));

    dt->SetProbeEtaPhiEn(inpData.ppfjet_eta, inpData.ppfjet_phi,
			 ecalE+hcalE_noRecHits,
			 inpData.getHcalEMap(0,1e-4));


    h1_dPhi->Fill( dt->GetTagPhi() - dt->GetProbePhi() , w);
    h1_dPt_PF->Fill( inpData.tpfjet_pt - inpData.ppfjet_pt , w);
    //h1_tjet_tightID->Fill( inpData.passTightJetID(1), w);
    //h1_pjet_tightID->Fill( inpData.passTightJetID(0), w);
    //h1_jet_energyGen->Fill( inpData.ppfjet_genE, w);
    //h1_jet_energyPF->Fill( inpData.ppfjet_E, w);
    //h1_jet_energyRH->Fill( dt->GetProbeEtot(), w);

    h1_tjet_energy_PFoverGen->Fill( inpData.tpfjet_E/inpData.tpfjet_genE, w);
    h1_tjet_energy_RHoverGen->Fill( dt->GetTagEtot()/inpData.tpfjet_genE, w);
    h1_pjet_energy_PFoverGen->Fill( inpData.ppfjet_E/inpData.ppfjet_genE, w);
    h1_pjet_energy_RHoverGen->Fill( dt->GetProbeEtot()/inpData.ppfjet_genE, w);

    //h1_jet_pT_PFoverGen->Fill( inpData.ppfjet_pt/inpData.ppfjet_genpt, w);
    //h1_jet_pT_RHoverGen->Fill( dt->GetProbeETtot()/inpData.ppfjet_genpt, w);

    h2_tjet_ptPF->Fill( inpData.tpfjet_genpt, inpData.tpfjet_pt, w);
    h2_pjet_ptPF->Fill( inpData.ppfjet_genpt, inpData.ppfjet_pt, w);
    h2_tjet_ptRH->Fill( inpData.tpfjet_genpt, dt->GetTagETtot(), w );
    h2_pjet_ptRH->Fill( inpData.ppfjet_genpt, dt->GetProbeETtot(), w );

    h2_jet_vs_jet_genPt->Fill( inpData.tpfjet_genpt, inpData.ppfjet_genpt, w);
    //h2_pho_vs_jet_recoPtPF->Fill( inpData.tagPho_pt, inpData.ppfjet_pt, w);
    //h2_pho_vs_jet_recoPtRH->Fill( inpData.tagPho_pt, dt->GetProbeETtot(), w);
    h2_jet_vs_jet_recoEta->Fill( inpData.tpfjet_eta, inpData.ppfjet_eta, w);

    /*
    int gen_debug=0;
    if (gen_debug==1) {
      std::cout << "iEntry=" << iEntry << ",    leading jet genE=" << inpData.ppfjet_genE << ", pfE=" << inpData.ppfjet_E << ", rhE=" << dt->GetProbeEtot() << "\n";
      std::cout << " -- " << dt->CalcProbeHcalEtot() << " + " << dt->GetProbeEcalE() << "\n";
      std::cout << "  ... " << dt->GetProbeHcalE() << "\n";
      std::cout << "  chk " << inpData.getHcalEMap(1,0.) << "\n";

      std::cout << "  ppfjet_twr_ieta " << (inpData.ppfjet_twr_ieta) << "\n";
      std::cout << "  ppfjet_twr_hade " << (inpData.ppfjet_twr_hade) << "\n";
      std::cout << "  ppfjet_twr_frac " << (inpData.ppfjet_twr_frac) << "\n";
      std::cout << "  ppfjet_twr_clusterind " << (inpData.ppfjet_twr_clusterind) << "\n";
      std::cout << "  ppfjet_cluster_dR " << (inpData.ppfjet_cluster_dR) << "\n";
    }
    else if (gen_debug==2) {
      std::cout << "genE=" << inpData.ppfjet_genE << ", pfE=" << inpData.ppfjet_E << ", rhE=" << dt->GetProbeEtot() << "\n";
    }


    if (gen_debug) {
      GammaJetEvent_t dt_sublead;
      GammaJetEventAuxInfo_t aux2(aux);
      double ecalE2=inpData.getSumEcalE_id(_subleadJet,1);
      double hcalE_noRecHits2=inpData.getSumHcalE_trackDiffEcal_id(_subleadJet);
      aux2.SetProbeHcalENoRecHits(hcalE_noRecHits2);
      dt_sublead.SetAuxInfo(aux2);
      dt_sublead.SetWeight(w);
      dt_sublead.SetTagEEtaPhi(0.,0.,0.); // irrelevant
      dt_sublead.SetProbeEtaPhiEn(inpData.pfjet2_eta, inpData.pfjet2_phi,
				  ecalE2+hcalE_noRecHits2,
				  inpData.getHcalEMap_id(_subleadJet,1e-4));

      if (gen_debug==1) {
	std::cout << "iEntry=" << iEntry << ", subleading jet genE=" << inpData.pfjet2_genE << ", pfE=" << inpData.pfjet2_E << ", rhE=" << dt_sublead.GetProbeEtot() << "\n";
      }
      else {
	std::cout << "genE=" << inpData.pfjet2_genE << ", pfE=" << inpData.pfjet2_E << ", rhE=" << dt_sublead.GetProbeEtot() << "\n";
      }
    }
    */

  }

  std::cout << "nEntries=" << nEntries << ", passedCount=" << passedCount << "\n";

  displayHisto(show_dPhi,h1_dPhi,"dPhi","LPE");
  displayHisto(show_dPt_PF,h1_dPt_PF, "dPt_PF","LPE");
  //displayHisto(show_jet_tightID,h1_jet_tightID,"jet_tightID","LPE");
  displayHisto(show_jet_energy_PF_over_Gen,h1_tjet_energy_PFoverGen,"te_PFoverGen","hist");
  displayHisto(show_jet_energy_PF_over_Gen,h1_tjet_energy_RHoverGen,"te_RHoverGen","hist");
  displayHisto(show_jet_energy_PF_over_Gen,h1_pjet_energy_PFoverGen,"pe_PFoverGen","hist");
  displayHisto(show_jet_energy_PF_over_Gen,h1_pjet_energy_RHoverGen,"pe_RHoverGen","hist");

  /*
  if (show_jet_pT_PF_over_Gen) {
    displayHisto(h1_jet_pT_PFoverGen,"pT_PFoverGen","hist");
    displayHisto(h1_jet_pT_RHoverGen,"pT_RHoverGen","hist");
  }
  */

  displayHisto(show_jet2D_ptPF,h2_tjet_ptPF,"1stjet2D_ptPF","COLZ");
  displayHisto(show_jet2D_ptPF,h2_pjet_ptPF,"2ndjet2D_ptPF","COLZ");

  displayHisto(show_jet2D_ptRH,h2_tjet_ptRH,"1stjet2D_ptRH","COLZ");
  displayHisto(show_jet2D_ptRH,h2_pjet_ptRH,"2ndjet2D_ptRH","COLZ");

  displayHisto(show_jet_vs_jet_genPt,h2_jet_vs_jet_genPt,"jet_vs_jet_genPt","COLZ");
  //if (show_pho_vs_jet_recoPtPF) displayHisto(h2_pho_vs_jet_recoPtPF,"pho_vs_jet_recoPtPF","COLZ");
  //if (show_pho_vs_jet_recoPtRH) displayHisto(h2_pho_vs_jet_recoPtRH,"pho_vs_jet_recoPtRH","COLZ");
  displayHisto(show_jet_vs_jet_recoEta,h2_jet_vs_jet_recoEta,"jet_vs_jet_recoEta","COLZ");

  /*
  if (show_jet_energy_PF_over_Gen_sideBySide) {
    displayHistoRatio(h1_jet_energyGen,"gen",h1_jet_energyPF,"PF",
		      "E [GeV]","count", "jet_en_gen_over_pf",
		      "hist","hist");
    displayHistoRatio(h1_jet_energyGen,"gen",h1_jet_energyRH,"RecHits",
		      "E [GeV]","count", "jet_en_gen_over_rh",
		      "hist","hist");
  }
  */

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
		       TString drawOpt1, TString drawOpt2) {
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
