#include "../interface/GammaJetFitData.h"
#include "pf_gammajettree.h"
#include <TH2D.h>
#include <TCanvas.h>
//#include "colorPalettes.hh"

// ------------------------------------------------------

inline void prepareHisto(TH1D *h) {
  h->SetDirectory(0);
}

inline void prepareHisto(TH2D *h) {
  h->SetDirectory(0);
}

// ------------------------------------------------------

inline
int passCuts(const pf_gammajettree &e, int theSet=0) {
  if (theSet==0) return 1;
  if (e.tagPho_idLoose==0) return 0;
  if (theSet==1) return 1;
  if (e.tagPho_idTight==0) return 0;
  return 1;
}

// ------------------------------------------------------

void displayHisto(TH1D* h, TString tag, TString drawOpt="LPE");
void displayHisto(TH2D* h, TString tag, TString drawOpt="COLZ");

// ------------------------------------------------------
// ------------------------------------------------------

void analyzePFGammaJetTree(TString inpFileName,
			   Long64_t maxEntries=-1,
			   double extraWeightFactor=-1.) {

  // book histograms
  double cPI= 4*atan(1);
  int show_dPhi=0;
  TH1D *h1_dPhi= new TH1D("h1_dPhi","#Delta#Phi",100,-2*cPI,2*cPI);
  prepareHisto(h1_dPhi);

  int show_dPt_PF=0;
  TH1D *h1_dPt_PF= new TH1D("h1_dPt_PF","#Deltap_{T}^{PF}",100,-1000,1000);
  prepareHisto(h1_dPt_PF);

  int show_pho2D_pt=1;
  TH2D *h2_pho_pt= new TH2D("h2_pho_pt","Photon p_{T};gen p_{T};reco p_{T}",
			    100,0.,1000.,
			    100,0.,1000.);
  prepareHisto(h2_pho_pt);

  int show_pho2D_ptMap=1;
  TH2D *h2_pho_ptMap= (TH2D*)h2_pho_pt->Clone("h2_pho_ptMap");
  h2_pho_ptMap->SetTitle("Photon average quality");
  TH2D *h2_pho_ptMapCount= (TH2D*)h2_pho_pt->Clone("h2_pho_ptMapCount");
  prepareHisto(h2_pho_ptMap);
  prepareHisto(h2_pho_ptMapCount);

  int show_jet2D_ptPF=1;
  TH2D *h2_jet_ptPF= new TH2D("h2_jet_pt",
			      "Leading jet p_{T};gen p_{T};reco p_{T}^{PF}",
			      100,0.,1000.,
			      100,0.,1000.);
  prepareHisto(h2_jet_ptPF);

  int show_pho_vs_jet_genPt=1;
  TH2D *h2_pho_vs_jet_genPt= new TH2D("h2_pho_vs_jet_genPt",
				   "Gen p_{T};#gamma gen p_{T};jet gen p_{T}",
				   100,0.,1000.,
				   100,0.,1000.);
  prepareHisto(h2_pho_vs_jet_genPt);

  int show_pho_vs_jet_recoPtPF=1;
  TH2D *h2_pho_vs_jet_recoPtPF= new TH2D("h2_pho_vs_jet_recoPtPF",
			   "Reco p_{T};#gamma p_{T};jet p_{T}^{PF}",
				   100,0.,1000.,
				   100,0.,1000.);
  prepareHisto(h2_pho_vs_jet_recoPtPF);


  // process the data
  pf_gammajettree inpData(inpFileName);

  inpData.DeactivateBranches();
  inpData.ActivateBranches_forFitSkim();
  inpData.ActivateBranches(1,"ppfjet_pt");
  inpData.ActivateBranches_genBasicSet();

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
    if (iEntry%10000==0) std::cout << " ... reading entry " << iEntry << "\n";
    std::cout << "ientry=" << iEntry << "\n";

    if (!passCuts(inpData,0)) continue;

    passedCount++;

    double ecalE= inpData.getSumEcalE(0,1);
    double hcalE_noRecHits= inpData.getSumHcalE_trackDiffEcal(1);

    aux.SetEventNo(inpData.EventNumber);
    aux.SetRunNo(inpData.RunNumber);
    aux.SetProbeHcalENoRecHits(hcalE_noRecHits);
    dt->SetAuxInfo(aux);

    double w=inpData.EventWeight;
    if (extraWeightFactor>0) w*= extraWeightFactor;
    dt->SetWeight(w);

    dt->SetTagEEtaPhi(inpData.tagPho_energy,
		      inpData.tagPho_eta, inpData.tagPho_phi);

    dt->SetProbeEtaPhiEn(inpData.ppfjet_eta, inpData.ppfjet_phi,
			 ecalE+hcalE_noRecHits, inpData.getHcalEMap(1,1e-4));

    h1_dPhi->Fill( dt->GetTagPhi() - dt->GetProbePhi() , w);
    h1_dPt_PF->Fill( inpData.tagPho_pt - inpData.ppfjet_pt , w);
    h2_pho_pt->Fill( inpData.tagPho_genPt, inpData.tagPho_pt, w);
    h2_jet_ptPF->Fill( inpData.ppfjet_genpt, inpData.ppfjet_pt, w);
    h2_pho_vs_jet_genPt->Fill( inpData.tagPho_genPt, inpData.ppfjet_genpt, w);
    h2_pho_vs_jet_recoPtPF->Fill( inpData.tagPho_pt, inpData.ppfjet_pt, w);

    int cathegory=0;
    if (inpData.tagPho_idTight) cathegory=2;
    else if (inpData.tagPho_idLoose) cathegory=1;
    h2_pho_ptMap->Fill( inpData.tagPho_genPt, inpData.tagPho_pt, cathegory);
    h2_pho_ptMapCount->Fill( inpData.tagPho_genPt, inpData.tagPho_pt, 1);
  }

  std::cout << "nEntries=" << nEntries << ", passedCount=" << passedCount << "\n";

  if (show_dPhi) displayHisto(h1_dPhi,"dPhi","LPE");
  if (show_dPt_PF) displayHisto(h1_dPt_PF, "dPt_PF","LPE");
  if (show_pho2D_pt) displayHisto(h2_pho_pt,"pho2D_pt","COLZ");
  if (show_jet2D_ptPF) displayHisto(h2_jet_ptPF,"jet2D_ptPF","COLZ");
  if (show_pho_vs_jet_genPt) displayHisto(h2_pho_vs_jet_genPt,"pho_vs_jet_genPt","COLZ");
  if (show_pho_vs_jet_recoPtPF) displayHisto(h2_pho_vs_jet_recoPtPF,"pho_vs_jet_recoPtPF","COLZ");

  if (show_pho2D_ptMap) {
    for (int ibin=1; ibin<=h2_pho_ptMapCount->GetNbinsX(); ++ibin) {
      for (int jbin=1; jbin<=h2_pho_ptMapCount->GetNbinsY(); ++jbin) {
	double cnt= h2_pho_ptMapCount->GetBinContent(ibin,jbin);
	if (cnt==double(0)) continue;
	double sum= h2_pho_ptMap->GetBinContent(ibin,jbin);
	h2_pho_ptMap->SetBinContent(ibin,jbin, sum/cnt);
      }
    }
    displayHisto(h2_pho_ptMap,"pho_ptMap","COLZ");
  }

  return;
}

// ------------------------------------------------------

void displayHisto(TH1D* h, TString tag, TString drawOpt) {
  TString canvName="c_" + tag;
  TCanvas *c=new TCanvas(canvName,canvName,600,600);
  h->Draw(drawOpt);
  c->Update();

}

// ------------------------------------------------------

void displayHisto(TH2D* h, TString tag, TString drawOpt) {
  TString canvName="c_" + tag;
  TCanvas *c=new TCanvas(canvName,canvName,600,600);
#ifdef ColorPalettes_HH
  AdjustFor2DplotWithHeight(c);
#endif
  h->Draw(drawOpt);
  c->Update();

}

// ------------------------------------------------------
