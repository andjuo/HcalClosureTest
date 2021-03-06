#include "../interface/GammaJetFitData.h"
#include "pf_gammajettree.h"
#include <TH2D.h>
#include <TCanvas.h>
#include <TLine.h>
//#include "colorPalettes.hh"
#include "ComparisonPlot.hh"
#include "../interface/HistoCollector.h"
#include <TLorentzVector.h>
#include "helper.h"

// ------------------------------------------------------

#ifdef helper_HH
  TString plotOutDir="plots";
#endif

// ------------------------------------------------------

inline void prepareHisto(TH1D *h, int showStats=0) {
  h->SetDirectory(0);
  h->SetStats(showStats);
}

inline void prepareHisto(TH2D *h, int showStats=0) {
  h->SetDirectory(0);
  h->SetStats(showStats);
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

void analyzePFGammaJetTree(TString inpFileName,
			   Long64_t maxEntries=-1,
			   double extraWeightFactor=-1.,
			   TString plotOutDir_user="",
			   int saveCollection_user=-1,
			   TString restrictTriggers="|") {

#ifdef helper_HH
  if (plotOutDir_user.Length()>0) plotOutDir= plotOutDir_user;
#endif
  if (saveCollection_user>=0) saveCollection=saveCollection_user;

  collector.Clear();

  std::vector<int> chkPhoTriggers,chkJetTriggers;
  if (!identifyTriggerIndices(restrictTriggers,
			      chkPhoTriggers,chkJetTriggers,0)) return;

  // book histograms
  double cPI= 4*atan(1);
  int show_dPhi=1;
  TH1D *h1_dPhi= new TH1D("h1_dPhi","#Delta#Phi(#gamma,j);#Delta#Phi;count",100,-2*cPI,2*cPI);
  TH1D *h1_dPhi_Abs= new TH1D("h1_dPhi_Abs","Absolute value of |#Delta#Phi(#gamma,j)|;|#Delta#Phi|;count",100,0,2*cPI);
  prepareHisto(h1_dPhi);
  prepareHisto(h1_dPhi_Abs);

  int show_pho_vs_jet_phi=1;
  TH2D *h2_pho_vs_jet_phi= new TH2D("h2_pho_jet_phi",
		  "#Phi correlation after selection;#phi_{#gamma};#phi_{jet}",
				 13,-1.5*cPI,1.5*cPI,
				 13,-1.5*cPI,1.5*cPI);
  prepareHisto(h2_pho_vs_jet_phi);

  double c_pt_max=500;
  double c_energy_max=1500;

  int show_dPt_PF=1;
  TH1D *h1_dPt_PF= new TH1D("h1_dPt_PF",
			    "#Deltap_{T}^{PF};#Deltap_{T}^{PF};count",
			    100,-c_pt_max,c_pt_max);
  prepareHisto(h1_dPt_PF);

  int show_jet_tightID=1;
  TH1D* h1_jet_tightID= new TH1D("h1_jet_tightID",
				 "Jet tightID flag;flag;count",
				 3,0.,2.);
  prepareHisto(h1_jet_tightID);

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


  int show_jet_energy_PF_over_Gen=1;
  TH1D *h1_jet_energy_PFoverGen= new TH1D("h1_jet_energy_PFoverGen",
		       "Jet energy ratio (PF div gen);E_{PF}/E_{gen};count",
					  100,0.,2.5);
  TH1D *h1_jet_energy_RHoverGen= new TH1D("h1_jet_energy_RHoverGen",
		       "Jet energy ratio (RH div gen);E_{RH}/E_{gen};count",
					  100,0.,2.5);
  prepareHisto(h1_jet_energy_PFoverGen,1);
  prepareHisto(h1_jet_energy_RHoverGen,1);

  int show_jet_pT_PF_over_Gen=1;
  TH1D *h1_jet_pT_PFoverGen= new TH1D("h1_jet_pT_PFoverGen",
	       "Jet E_{T} ratio (PF div gen);p_{T}^{PF}/p_{T}^{gen};count",
					  100,0.,2.5);
  TH1D *h1_jet_pT_RHoverGen= new TH1D("h1_jet_pT_RHoverGen",
	       "Jet E_{T} ratio (RH div gen);p_{T}^{RH}/p_{T}^{gen};count",
					  100,0.,2.5);
  prepareHisto(h1_jet_pT_PFoverGen,1);
  prepareHisto(h1_jet_pT_RHoverGen,1);

  int show_pho2D_pt=1;
  TH2D *h2_pho_pt= new TH2D("h2_pho_pt","Photon p_{T};gen p_{T};reco p_{T}",
			    100,0.,c_pt_max,
			    100,0.,c_pt_max);
  prepareHisto(h2_pho_pt);

  int show_pho2D_ptMap=0;
  TH2D *h2_pho_ptMap= (TH2D*)h2_pho_pt->Clone("h2_pho_ptMap");
  h2_pho_ptMap->SetTitle("Photon average quality");
  TH2D *h2_pho_ptMapCount= (TH2D*)h2_pho_pt->Clone("h2_pho_ptMapCount");
  prepareHisto(h2_pho_ptMap);
  prepareHisto(h2_pho_ptMapCount);

  int show_jet2D_ptPF=1;
  TH2D *h2_jet_ptPF= new TH2D("h2_jet_ptPF",
			      "Leading jet p_{T};gen p_{T};reco p_{T}^{PF}",
			      100,0.,c_pt_max,
			      100,0.,c_pt_max);
  prepareHisto(h2_jet_ptPF);

  int show_jet2D_ptRH=1;
  TH2D *h2_jet_ptRH= new TH2D("h2_jet_ptRH",
			      "Leading jet p_{T};gen p_{T};reco p_{T}^{RecHits}",
			      100,0.,c_pt_max,
			      100,0.,c_pt_max);
  prepareHisto(h2_jet_ptRH);

  int show_pho_pt=1;
  TH1D *h1_pho_genPt= new TH1D("h1_pho_genPt",
			       "Photon gen p_{T};gen p_{#gamma,T};count",
			       100,0.,c_pt_max);
  TH1D *h1_pho_recoPt= new TH1D("h1_pho_recoPt",
			       "Photon reco p_{T};reco p_{#gamma,T};count",
			       100,0.,c_pt_max);
  prepareHisto(h1_pho_genPt,1);
  prepareHisto(h1_pho_recoPt,1);

  int show_pho_vs_jet_genPt=1;
  TH2D *h2_pho_vs_jet_genPt= new TH2D("h2_pho_vs_jet_genPt",
				   "Gen p_{T};#gamma gen p_{T};jet gen p_{T}",
				   100,0.,c_pt_max,
				   100,0.,c_pt_max);
  prepareHisto(h2_pho_vs_jet_genPt);

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

  int show_pho_vs_jet_recoEta=1;
  TH2D *h2_pho_vs_jet_recoEta= new TH2D("h2_pho_vs_jet_recoEta",
					"Reco #eta;#gamma #eta;jet #eta",
					100,-2.5*cPI, 2.5*cPI,
					100,-2.5*cPI, 2.5*cPI);
  prepareHisto(h2_pho_vs_jet_recoEta);

  int show_jet_div_gamma_pT_Reco=1;
  int show_jet_div_gamma_pT_Reco_chk=0;
  TH1D *h1_jet_div_gamma_pTPF_Reco= new TH1D("h1_jet_div_gamma_pTPF_Reco",
      "Jet E_{T} ratio to #gamma p_{T} (PF);p_{jet,T}^{PF}/p_{#gamma,T};count",
					     100,0.,2.5);
  TH1D *h1_jet_div_gamma_pTRH_Reco= new TH1D("h1_jet_div_gamma_pTRH_Reco",
      "Jet E_{T} ratio to #gamma p_{T} (RH);p_{jet,T}^{RH}/p_{#gamma,T};count",
					     100,0.,2.5);
  TH1D *h1_jet_div_gamma_pTRH_Reco_chk= new TH1D("h1_jet_div_gamma_pTRH_Reco_chk",
      "Jet E_{T} ratio to #gamma p_{T} (RH) control plot;p_{jet,T}^{RH}/p_{#gamma,T};count",
					     100,0.,2.5);
  prepareHisto(h1_jet_div_gamma_pTPF_Reco,1);
  prepareHisto(h1_jet_div_gamma_pTRH_Reco,1);
  prepareHisto(h1_jet_div_gamma_pTRH_Reco_chk,1);

  // jet vs gamma+addedJet2
  int show_jet_div_gammaJet2_pT_Reco=1;
  TH1D *h1_jet_div_gammaJet2_pTPF_Reco=
    new TH1D("h1_jet_div_gammaJet2_pTPF_Reco",
      "Jet E_{T} ratio to (#gamma+jet2) p_{T} (PF);p_{jet,T}^{PF}/(p_{#gamma,T}+p_{jet2,T});count",
	     100,0.,2.5);
  TH1D *h1_jet_div_gammaJet2_pTRH_Reco=
    new TH1D("h1_jet_div_gammaJet2_pTRH_Reco",
      "Jet E_{T} ratio to (#gamma+jet2) p_{T} (RH);p_{jet,T}^{RH}/(p_{#gamma,T}+p_{jet2,T});count",
	     100,0.,2.5);
  TH1D *h1_jetJet2_div_gamma_pTPF_Reco=
    new TH1D("h1_jetJet2_div_gamma_pTPF_Reco",
      "(Jet+jet2) E_{T} ratio to #gamma p_{T} (PF);(p_{jet,T}^{PF}+p_{jet2,T})/p_{#gamma,T};count",
	     100,0.,2.5);

  prepareHisto(h1_jet_div_gammaJet2_pTPF_Reco,1);
  prepareHisto(h1_jet_div_gammaJet2_pTRH_Reco,1);
  prepareHisto(h1_jetJet2_div_gamma_pTPF_Reco,1);

  int show_jet_div_gamma_Hemi=1;
  TH1D *h1_jet_div_gamma_pTPF_jetHemi_Reco= new TH1D("h1_jet_div_gamma_pTPF_jetHemi_Reco",
      "Jet E_{T} ratio to #gamma p_{T} (PF), 2nd jet in jet Hemi;p_{jet,T}^{PF}/p_{#gamma,T} | (2nd jet jetHemi);count",
					     100,0.,2.5);
  TH1D *h1_jet_div_gamma_pTPF_phoHemi_Reco= new TH1D("h1_jet_div_gamma_pTPF_phoHemi_Reco",
      "Jet E_{T} ratio to #gamma p_{T} (PF), 2nd jet in pho Hemi;p_{jet,T}^{PF}/p_{#gamma,T} | (2nd jet phoHemi);count",
					     100,0.,2.5);

  TH1D *h1_jetJet2_div_gamma_pTPF_jetHemi_Reco=
    new TH1D("h1_jetJet2_div_gamma_pTPF_jetHemi_Reco",
	     "(Jet+jet2) E_{T} ratio to #gamma p_{T} (PF,jetHemi);(p_{jet,T}^{PF}+p_{jet2,T})/p_{#gamma,T};count",
	     100,0.,2.5);
  TH1D *h1_jet_div_gammaJet2_pTPF_phoHemi_Reco=
    new TH1D("h1_jet_div_gammaJet2_pTPF_phoHemi_Reco",
	     "Jet E_{T} ratio to (#gamma+jet2) p_{T} (PF,phoHemi);p_{jet,T}^{PF}/(p_{#gamma,T}+p_{jet2,T}^{PF});count",
	     100,0.,2.5);

  prepareHisto(h1_jet_div_gamma_pTPF_jetHemi_Reco,1);
  prepareHisto(h1_jet_div_gamma_pTPF_phoHemi_Reco,1);
  prepareHisto(h1_jetJet2_div_gamma_pTPF_jetHemi_Reco,1);
  prepareHisto(h1_jet_div_gammaJet2_pTPF_phoHemi_Reco,1);

  // process the data
  pf_gammajettree inpData(inpFileName);

  if (1) { // deactivation speeds-up file reading
    inpData.DeactivateBranches();
    inpData.ActivateBranches_forFitSkim();
    inpData.ActivateBranches_genBasicSet();
  }


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

    if (!inpData.passCuts(20.,int(_phoTightID),1)) continue; // recommended use
    //if ( inpData.pfjet2_pt/inpData.tagPho_pt > 0.05) continue; // stricter alpha cut

    if (1) {
      if ((fabs(inpData.tagPho_eta)>2.4) || (fabs(inpData.ppfjet_eta)>2.4))
	continue;

      if ((restrictTriggers.Length()>1) &&
	  !inpData.hltTriggerFired(chkPhoTriggers,chkJetTriggers)) continue;
    }


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

    double dphi_abs= fabs(dt->GetTagPhi() - dt->GetProbePhi());
    if (dphi_abs>cPI) dphi_abs= 2*cPI - dphi_abs;

    h1_dPhi->Fill( dt->GetTagPhi() - dt->GetProbePhi() , w);
    h1_dPhi_Abs->Fill ( dphi_abs, w );
    h2_pho_vs_jet_phi->Fill( inpData.tagPho_phi, inpData.ppfjet_phi, w);
    h1_dPt_PF->Fill( inpData.tagPho_pt - inpData.ppfjet_pt , w);
    h1_jet_tightID->Fill( inpData.passTightJetID(1), w);
    h1_jet_energyGen->Fill( inpData.ppfjet_genE, w);
    h1_jet_energyPF->Fill( inpData.ppfjet_E, w);
    h1_jet_energyRH->Fill( dt->GetProbeEtot(), w);
    h1_jet_energy_PFoverGen->Fill( inpData.ppfjet_E/inpData.ppfjet_genE, w);
    h1_jet_energy_RHoverGen->Fill( dt->GetProbeEtot()/inpData.ppfjet_genE, w);
    h1_jet_pT_PFoverGen->Fill( inpData.ppfjet_pt/inpData.ppfjet_genpt, w);
    h1_jet_pT_RHoverGen->Fill( dt->GetProbeETtot()/inpData.ppfjet_genpt, w);
    h2_pho_pt->Fill( inpData.tagPho_genPt, inpData.tagPho_pt, w);
    h2_jet_ptPF->Fill( inpData.ppfjet_genpt, inpData.ppfjet_pt, w);
    //h2_jet_ptPF->Fill( inpData.pfjet2_genpt, inpData.pfjet2_pt, w);
    h2_jet_ptRH->Fill( inpData.ppfjet_genpt, dt->GetProbeETtot(), w );
    h1_pho_genPt->Fill (inpData.tagPho_genPt, w );
    h1_pho_recoPt->Fill(inpData.tagPho_pt, w );
    h2_pho_vs_jet_genPt->Fill( inpData.tagPho_genPt, inpData.ppfjet_genpt, w);
    h2_pho_vs_jet_recoPtPF->Fill( inpData.tagPho_pt, inpData.ppfjet_pt, w);
    h2_pho_vs_jet_recoPtRH->Fill( inpData.tagPho_pt, dt->GetProbeETtot(), w);
    h2_pho_vs_jet_recoEta->Fill( inpData.tagPho_eta, inpData.ppfjet_eta, w);

    h1_jet_div_gamma_pTPF_Reco->Fill( inpData.ppfjet_pt/inpData.tagPho_pt, w );
    h1_jet_div_gamma_pTRH_Reco->Fill( dt->GetProbeETtot()/inpData.tagPho_pt, w );
    h1_jet_div_gamma_pTRH_Reco_chk->Fill( dt->GetProbeETtot() / dt->GetTagETtot(), w );

    if (1) {
      // add gamma and 2nd jet
      GammaJetEvent_t dt2RH_gammaJet2(*dt);
      // add the jet contribution to the transverse momentum
      if (0) {
	// pfjet PF energy is different from pT
	if ( inpData.pfjet2_E / cosh(inpData.pfjet2_eta) != inpData.pfjet2_pt )
	  {
	    std::cout << "\tpfjet2_E check failed: "
		      << inpData.pfjet2_E/cosh(inpData.pfjet2_eta)
		      << " vs pt=" << inpData.pfjet2_pt << "\n";
	  }
      }
      dt2RH_gammaJet2.AddTagEEtaPhi( inpData.pfjet2_E,
			   inpData.pfjet2_eta, inpData.pfjet2_phi );
      double pt_gammaJet2= dt2RH_gammaJet2.GetTagETtot();

      h1_jet_div_gammaJet2_pTPF_Reco->Fill( inpData.ppfjet_pt/ pt_gammaJet2, w );
      h1_jet_div_gammaJet2_pTRH_Reco->Fill( dt->GetProbeETtot() / pt_gammaJet2, w );

      TLorentzVector jet,jet2,sum;
      jet .SetPtEtaPhiE(inpData.ppfjet_pt,inpData.ppfjet_eta,inpData.ppfjet_phi,inpData.ppfjet_E);
      jet2.SetPtEtaPhiE(inpData.pfjet2_pt,inpData.pfjet2_eta,inpData.pfjet2_phi,inpData.pfjet2_E);
      sum= jet+jet2;
      // take into account only the momentum change
      double jetJet2_pt=sqrt(pow(sum.Px(),2)+pow(sum.Py(),2));
      h1_jetJet2_div_gamma_pTPF_Reco->Fill( jetJet2_pt / inpData.tagPho_pt, w );

      double dPhi1 = fabs(inpData.tagPho_phi - inpData.pfjet2_phi);
      double dPhi2 = fabs(inpData.ppfjet_phi - inpData.pfjet2_phi);
      if (dPhi1 > cPI) dPhi1= 2*cPI - dPhi1;
      if (dPhi2 > cPI) dPhi2= 2*cPI - dPhi2;
      int photonHemi=(dPhi1<dPhi2) ? 1:0;

      //double phoPhi=inpData.tagPho_phi;
      //double jetPhi=inpData.ppfjet_phi;
      //double jet2Phi=inpData.pfjet2_phi;
      //std::cout << "phi: " << phoPhi << "," << jetPhi << "," << jet2Phi << "  ";
      //std::cout << "  dPhi:" << dPhi1 << ", " << dPhi2 << ": photonSide=" << photonHemi << "\n";
      if (photonHemi) {
	h1_jet_div_gamma_pTPF_phoHemi_Reco->Fill( inpData.ppfjet_pt / inpData.tagPho_pt, w );
	h1_jet_div_gammaJet2_pTPF_phoHemi_Reco->Fill( inpData.ppfjet_pt / pt_gammaJet2, w );
      }
      else {
	h1_jet_div_gamma_pTPF_jetHemi_Reco->Fill ( inpData.ppfjet_pt / inpData.tagPho_pt, w );
	h1_jetJet2_div_gamma_pTPF_jetHemi_Reco->Fill (jetJet2_pt / inpData.tagPho_pt, w );
      }
    }

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

    int cathegory=0;
    if (inpData.tagPho_idTight) cathegory=2;
    else if (inpData.tagPho_idLoose) cathegory=1;
    h2_pho_ptMap->Fill( inpData.tagPho_genPt, inpData.tagPho_pt, cathegory);
    h2_pho_ptMapCount->Fill( inpData.tagPho_genPt, inpData.tagPho_pt, 1);
  }

  std::cout << "nEntries=" << nEntries << ", passedCount=" << passedCount << "\n";
  //return;

  displayHisto(show_dPhi, h1_dPhi,"dPhi","LPE");
  displayHisto(show_dPhi, h1_dPhi_Abs,"dPhiAbs","LPE");
  displayHisto(show_pho_vs_jet_phi, h2_pho_vs_jet_phi,"pho_vs_jet_phi","COLZ");
  displayHisto(show_dPt_PF, h1_dPt_PF, "dPt_PF","LPE");
  displayHisto(show_jet_tightID, h1_jet_tightID,"jet_tightID","LPE");
  displayHisto(show_jet_energy_PF_over_Gen, h1_jet_energy_PFoverGen,"e_PFoverGen","hist");
  displayHisto(show_jet_energy_PF_over_Gen, h1_jet_energy_RHoverGen,"e_RHoverGen","hist");

  displayHisto(show_jet_pT_PF_over_Gen, h1_jet_pT_PFoverGen,"pT_PFoverGen","hist");
  displayHisto(show_jet_pT_PF_over_Gen, h1_jet_pT_RHoverGen,"pT_RHoverGen","hist");

  displayHisto(show_pho2D_pt, h2_pho_pt,"pho2D_pt","COLZ");
  displayHisto(show_jet2D_ptPF, h2_jet_ptPF,"jet2D_ptPF","COLZ");
  displayHisto(show_jet2D_ptRH, h2_jet_ptRH,"jet2D_ptRH","COLZ");
  displayHisto(show_pho_pt, h1_pho_genPt, "pho_genPt","LPE");
  displayHisto(show_pho_pt, h1_pho_recoPt, "pho_recoPt","LPE");
  displayHisto(show_pho_vs_jet_genPt, h2_pho_vs_jet_genPt,"pho_vs_jet_genPt","COLZ");
  displayHisto(show_pho_vs_jet_recoPtPF, h2_pho_vs_jet_recoPtPF,"pho_vs_jet_recoPtPF","COLZ");
  displayHisto(show_pho_vs_jet_recoPtRH, h2_pho_vs_jet_recoPtRH,"pho_vs_jet_recoPtRH","COLZ");
  displayHisto(show_pho_vs_jet_recoEta,h2_pho_vs_jet_recoEta,"pho_vs_jet_recoEta","COLZ");

  if (show_pho2D_ptMap) {
    for (int ibin=1; ibin<=h2_pho_ptMapCount->GetNbinsX(); ++ibin) {
      for (int jbin=1; jbin<=h2_pho_ptMapCount->GetNbinsY(); ++jbin) {
	double cnt= h2_pho_ptMapCount->GetBinContent(ibin,jbin);
	if (cnt==double(0)) continue;
	double sum= h2_pho_ptMap->GetBinContent(ibin,jbin);
	h2_pho_ptMap->SetBinContent(ibin,jbin, sum/cnt);
      }
    }
    displayHisto(show_pho2D_ptMap, h2_pho_ptMap,"pho_ptMap","COLZ");
  }

  if (show_jet_energy_PF_over_Gen_sideBySide) {
    displayHistoRatio(h1_jet_energyGen,"gen",h1_jet_energyPF,"PF",
		      "E [GeV]","count", "jet_en_gen_over_pf",
		      "hist","hist");
    displayHistoRatio(h1_jet_energyGen,"gen",h1_jet_energyRH,"RecHits",
		      "E [GeV]","count", "jet_en_gen_over_rh",
		      "hist","hist");
  }

  displayHisto(show_jet_div_gamma_pT_Reco, h1_jet_div_gamma_pTPF_Reco, "jet_div_pho_pT_PF","LPE");
  displayHisto(show_jet_div_gamma_pT_Reco, h1_jet_div_gamma_pTRH_Reco, "jet_div_pho_pT_RH","LPE");
  displayHisto(show_jet_div_gamma_pT_Reco_chk, h1_jet_div_gamma_pTRH_Reco_chk, "jet_div_pho_pT_RH_chk","LPE");

  displayHisto(show_jet_div_gammaJet2_pT_Reco, h1_jet_div_gammaJet2_pTPF_Reco, "jet_div_phoJet2_pT_PF","LPE");
  displayHisto(show_jet_div_gammaJet2_pT_Reco, h1_jet_div_gammaJet2_pTRH_Reco, "jet_div_phoJet2_pT_RH","LPE");
  displayHisto(show_jet_div_gammaJet2_pT_Reco, h1_jetJet2_div_gamma_pTPF_Reco, "jetJet2_div_pho_pT_PF","LPE");

  displayHisto(show_jet_div_gamma_Hemi, h1_jet_div_gamma_pTPF_jetHemi_Reco, "jet_div_gamma_pTPF_jetHemi","LPE");
  displayHisto(show_jet_div_gamma_Hemi, h1_jet_div_gamma_pTPF_phoHemi_Reco, "jet_div_gamma_pTPF_phoHemi","LPE");
  displayHisto(show_jet_div_gamma_Hemi, h1_jetJet2_div_gamma_pTPF_jetHemi_Reco, "jetJet2_div_pho_pTPF_jetHemi","LPE");
  displayHisto(show_jet_div_gamma_Hemi, h1_jet_div_gammaJet2_pTPF_phoHemi_Reco, "jet_div_phoJet2_pTPF_phoHemi","LPE");



  if (saveCollection) {
    TString outFName=TString("saved_") + plotOutDir + TString(".root");
    collector.Add("producedBy_analyzePFGammaJetTree");
    collector.SaveToFile(outFName);
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
