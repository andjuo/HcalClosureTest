#include "DijetRespCorrDataExtended.h"
#include "helper.hh"
#include "Analyzer.hh"
#include "colorPalettes.hh"
#include <TFile.h>
#include <TTree.h>
#include <assert.h>

// ----------------------------------------------------------------

void analyzeSkim(Long64_t maxEntries=10) {
  DijetRespCorrDatumExtended_t *d= new DijetRespCorrDatumExtended_t();
  TFile inpFile("skim.root","READ");
  TTree *inpTree = (TTree*)inpFile.Get("dijet_data");
  assert(inpTree);
  inpTree->SetBranchAddress("dijet_data",&d);

  const int plotBalanceVsCf=0;
  const int balanceKind=1; //  0 - original, 1 - met, 2 - metDivAvePt

  std::vector<Cuts_t*> cutV;
  std::vector<TH1D*> hTowV; // (iTow,Count)
  unsigned int balanceTargetCutIdx=0;
  std::vector<TH1D*> hTowBalance_iTow_V;
  std::vector<TH1D*> hTowBalance_iTow_SumW_V; // iTow: (Ci,Balance)
  std::vector<TH1D*> hTowBalance_iTow_SumXW_V; // iTow: (Ci,Balance)
  std::vector<TH1D*> hTowBalance_iTow_SumXXW_V; // iTow: (Ci,Balance)
  const int colorV[]= { kBlack, kRed, kGreen+1, kBlue, kViolet, 43 };
  std::vector<TH2D*> hBalanceEtDiffV;
  std::vector<TH2D*> hDiffRatioVsDeltaEotherV;

  cutV.reserve(5);
  if (0) {
    cutV.push_back(new Cuts_t("EtDiff_0_10"));
    cutV.back()->setEtDiff(0.,10., -1.,-1.);
    cutV.push_back(new Cuts_t("EtDiff_10_20"));
    cutV.back()->setEtDiff(10.,20., -1.,-1.);
    cutV.push_back(new Cuts_t("EtDiff_20_50"));
    cutV.back()->setEtDiff(20.,50., -1.,-1.);
    cutV.push_back(new Cuts_t("EtDiff_50_100"));
    cutV.back()->setEtDiff(50.,100., -1.,-1.);
    cutV.push_back(new Cuts_t("EtDiff_100_200"));
    cutV.back()->setEtDiff(100.,200., -1.,-1.);
    cutV.push_back(new Cuts_t("EtDiff_gt_200"));
    cutV.back()->setEtDiff(200.,10000, -1.,-1.);
  }

  balanceTargetCutIdx=cutV.size();  
  cutV.push_back(new Cuts_t("EtDiff_0_200"));
  //cutV.back()->setEtDiff(0,200.,-1.,-1.);
  cutV.back()->setEtDiffRatio(-2.,0.);


  if (1)
  for (unsigned int i=0; i<cutV.size(); ++i) {
    cutV[i]->swapEtFields("EtDiff","RHEtDiff");
  }

  for (unsigned int i=0; i<cutV.size(); ++i) {
    const double c_PI= 4*atan(1);
    cutV[i]->setDijetDiPhi(c_PI-0.2,c_PI+0.2);
  }

  // prepare histogram of count vs iTower
  for (unsigned int i=0; i<cutV.size(); ++i) {
    TString hname="hTow_" + cutV[i]->name();
    TString htitle= hname + TString(";ieta;count");
    TH1D *h= new TH1D(hname,htitle, 2*MAXIETA+2, -MAXIETA-1, MAXIETA+1);
    //printHisto(h); return;
    h->SetDirectory(0);
    hTowV.push_back(h);
  }

  // prepare histogram of balance vs iTower
  if (plotBalanceVsCf) {
    TString balanceStr;
    switch(balanceKind) {
    case 0: balanceStr="Balance"; break;
    case 1: balanceStr="MET"; break;
    case 2: balanceStr="METdivEt"; break;
    default: std::cout << "not ready for balanceKind=" << balanceKind << "\n";
      return;
    }
  for (int iTow=-MAXIETA; iTow<=MAXIETA; iTow++) {
    TString hName=Form("hTow%s_%d_",balanceStr.Data(),iTow)
      + cutV[balanceTargetCutIdx]->name();
    TString hTitle=hName +
      TString(Form(";C_{iTow=%d};%s [GeV]",iTow,balanceStr.Data()));
    int locCount=40;
    double cfMin=0.;
    double cfMax=2.;
    TH1D* h= new TH1D(hName,hTitle, locCount,cfMin,cfMax);
    h->SetDirectory(0);
    hTowBalance_iTow_V.push_back(h);

    hName.ReplaceAll(balanceStr,balanceStr+TString("SumW"));
    hTitle.ReplaceAll(balanceStr,balanceStr+TString("SumW"));
    TH1D* hSumW= new TH1D(hName,hTitle, locCount,cfMin,cfMax);
    hSumW->SetDirectory(0);
    hTowBalance_iTow_SumW_V.push_back(hSumW);

    hName.ReplaceAll("SumW","SumXW");
    hTitle.ReplaceAll("SumW","SumXW");
    TH1D* hSumXW= new TH1D(hName,hTitle, locCount,cfMin,cfMax);
    hSumXW->SetDirectory(0);
    hTowBalance_iTow_SumXW_V.push_back(hSumXW);

    hName.ReplaceAll("SumXW","SumXXW");
    hTitle.ReplaceAll("SumXW","SumXXW");
    TH1D* hSumXXW= new TH1D(hName,hTitle, locCount,cfMin,cfMax);
    hSumXXW->SetDirectory(0);
    hTowBalance_iTow_SumXXW_V.push_back(hSumXXW);
  }
  }
  // prepare the coefficients
  TArrayD unitArr(NUMTOWERS);
  for (int i=0; i<NUMTOWERS; ++i) unitArr[i]=1.;


  double maxEtDiff=300.;
  // prepare histograms for HadEn balance vs otherEn balance
  for (unsigned int i=0; i<cutV.size(); i++) {
    const int divCount=100;
    TString hname="hBalanceEtDiff_" + cutV[i]->name();
    TString htitle=hname + TString(";#DeltahadEt;#DeltanonHadEt");
    TH2D* h2= new TH2D(hname,htitle, divCount, -maxEtDiff, maxEtDiff,
		       divCount, -maxEtDiff, maxEtDiff);
    h2->SetDirectory(0);
    hBalanceEtDiffV.push_back(h2);
  }

  double maxDiffRatio=3.5;
  for (unsigned int i=0; i<cutV.size(); ++i) {
    const int divCount=50;
    TString hname="hDiffRatioVsDeltaEother_" + cutV[i]->name();
    TString htitle= hname +
          TString(";#DeltaE_{other};#DeltaE_{had}/#DeltaE_{other}");
    TH2D* h2= new TH2D(hname,htitle, divCount,-maxEtDiff,maxEtDiff,
		       divCount, -maxDiffRatio,maxDiffRatio);
    h2->SetDirectory(0);
    h2->SetStats(0);
    hDiffRatioVsDeltaEotherV.push_back(h2);
  }

  TH1D *h1TagJet_pfE_over_genE = new TH1D("h1TPFJet_pfE_over_genE",
		 "Tag jet pfE/genE;E_{jet1}^{PF};E_{jet1}^{gen}",100,0.,2.);
  TH1D *h1ProbeJet_pfE_over_genE = new TH1D("h1PPFJet_pfE_over_genE",
		 "Probe jet pfE/genE;E_{jet1}^{PF};E_{jet1}^{gen}",100,0.,2.);
  h1TagJet_pfE_over_genE->SetDirectory(0);
  h1ProbeJet_pfE_over_genE->SetDirectory(0);

  TH1D *h1TagJet_rhE_over_genE = new TH1D("h1TPFJet_rhE_over_genE",
		 "Tag jet rhE/genE;E_{jet1}^{RH};E_{jet1}^{gen}",100,0.,2.);
  TH1D *h1ProbeJet_rhE_over_genE = new TH1D("h1PPFJet_rhE_over_genE",
		 "Probe jet rhE/genE;E_{jet1}^{RH};E_{jet1}^{gen}",100,0.,2.);
  h1TagJet_rhE_over_genE->SetDirectory(0);
  h1ProbeJet_rhE_over_genE->SetDirectory(0);

  // --------------------------------------

  DijetRespCorrData DRCD;

  Long64_t nEntries= inpTree->GetEntries();
  if (maxEntries<0) maxEntries=nEntries;

  for (Long64_t iEntry=0; (iEntry<nEntries) && (iEntry<maxEntries);
       ++iEntry) {
    inpTree->GetEntry(iEntry);
    if (iEntry%10000==0) {
      std::cout << "iEntry=" << iEntry << "\n";
    }
    if (0) {
      std::cout << "print iEntry=" << iEntry << "\n";
      if (0) {
	std::cout << (*d) << "\n";
      }
      else d->PrintBaseClass(1);
    }

    h1TagJet_pfE_over_genE->Fill(d->fTPFJetE/d->fTPFJetGenE, d->GetWeight());
    h1ProbeJet_pfE_over_genE->Fill(d->fPPFJetE/d->fPPFJetGenE, d->GetWeight());
    h1TagJet_rhE_over_genE->Fill(d->CalcRecHits_TagJetEn()/d->fTPFJetGenE, d->GetWeight());
    h1ProbeJet_rhE_over_genE->Fill(d->CalcRecHits_ProbeJetEn()/d->fPPFJetGenE, d->GetWeight());

    std::cout << "iEntry=" << iEntry << ", tag   jet genE=" << d->fTPFJetGenE << ", pfE=" << d->fTPFJetE << ", rhE=" << d->CalcRecHits_TagJetEn() << "\n";
    std::cout << " -- " << d->fTJet_SumHadE << " + " << d->fTJet_SumHadOtherE << " + " << d->fTJet_SumNonHadOtherE << "\n";
    std::map<Int_t,Double_t> hm;
    hm.clear();
    d->GetTagHcalE(hm);
    for (std::map<Int_t,Double_t>::const_iterator it=hm.begin(); it!=hm.end(); it++) {
      std::cout << " (" << it->first << "," << it->second << ")";
    }
    std::cout << "\n";
    std::cout << "iEntry=" << iEntry << ", probe jet genE=" << d->fPPFJetGenE << ", pfE=" << d->fPPFJetE << ", rhE=" << d->CalcRecHits_ProbeJetEn() << "\n";
    std::cout << " -- " << d->fPJet_SumHadE << " + " << d->fPJet_SumHadOtherE << " + " << d->fPJet_SumNonHadOtherE << "\n";
    hm.clear();
    d->GetProbeHcalE(hm);
    for (std::map<Int_t,Double_t>::const_iterator it=hm.begin(); it!=hm.end(); it++) {
      std::cout << " (" << it->first << "," << it->second << ")";
    }
    std::cout << "\n";


    for (unsigned int iCut=0; iCut<cutV.size(); ++iCut) {
      if (cutV[iCut]->passCuts(*d)) {
	FillTowerCount(hTowV[iCut],*d,-1);

	double hadEtDiff= d->CalcRecHits_EnergyDiff(_sumHadEt);
	double nonHadEtDiff= d->CalcRecHits_EnergyDiff(_sumNonHadEt);
	double ratio=hadEtDiff/nonHadEtDiff;
	if (ratio<-maxDiffRatio) ratio=-maxDiffRatio+1e-3;
	if (ratio> maxDiffRatio) ratio= maxDiffRatio-1e-3;
	hDiffRatioVsDeltaEotherV[iCut]->Fill(nonHadEtDiff,ratio,d->GetWeight());

	if (hadEtDiff<-maxEtDiff) hadEtDiff = -maxEtDiff+1;
	if (hadEtDiff> maxEtDiff) hadEtDiff =  maxEtDiff-1;
	if (nonHadEtDiff<-maxEtDiff) nonHadEtDiff = -maxEtDiff+1;
	if (nonHadEtDiff> maxEtDiff) nonHadEtDiff =  maxEtDiff-1;
	hBalanceEtDiffV[iCut]->Fill(hadEtDiff,nonHadEtDiff, d->GetWeight());
      }


      if (plotBalanceVsCf)
      if ((cutV[iCut]->passCuts(*d)) && (iCut==balanceTargetCutIdx)) {
	for (int iTow=-MAXIETA; iTow<=MAXIETA; iTow++) {
	  //HERE("iTow=",iTow);
	  int idx=iTow+MAXIETA;
	  TH1D* hW= hTowBalance_iTow_SumW_V[idx];
	  TH1D* hXW= hTowBalance_iTow_SumXW_V[idx];
	  TH1D* hXXW= hTowBalance_iTow_SumXXW_V[idx];
	  double balance=0, balanceV2=0, resolution=0;
	  for (int ibin=1; ibin<=hW->GetNbinsX(); ++ibin) {
	    double cf= hW->GetBinCenter(ibin);
	    unitArr[idx]=cf;
	    switch(balanceKind) {
	    case 0: DRCD.GetBalance(*d,unitArr,balance,resolution); break;
	    case 1: DRCD.GetMET(*d,unitArr,balance,balanceV2,resolution);break;
	    case 2: DRCD.GetMET(*d,unitArr,balanceV2,balance,resolution);break;
	    default:
	      std::cout << "DRCD not ready for balanceKind="
			<< balanceKind << "\n";
	      return;
	    }
	    hW->Fill(cf, d->GetWeight());
	    hXW->Fill(cf, balance * d->GetWeight());
	    hXXW->Fill(cf, balance*balance * d->GetWeight());
	    unitArr[idx]=1.;
	  }
	}
      }
    }
  }
  inpFile.Close();
  HERE("file closed: ",inpFile.GetName());

  for (unsigned int i=0; i<hTowBalance_iTow_SumW_V.size(); ++i) {
    TH1D* h= hTowBalance_iTow_V[i];
    if (!h) HERE("h is null");
    printHisto(h);
    const TH1D* hW= hTowBalance_iTow_SumW_V[i];
    const TH1D* hXW= hTowBalance_iTow_SumXW_V[i];
    const TH1D* hXXW= hTowBalance_iTow_SumXXW_V[i];
    h->Reset();
    for (int ibin=1; ibin<= h->GetNbinsX(); ++ibin) {
      double sumW= hW->GetBinContent(ibin);
      double avg= (sumW!=0) ? hXW->GetBinContent(ibin)/sumW : 0.;
      double sigma=(sumW!=0) ? hXXW->GetBinContent(ibin)/sumW : 0.;
      if (sumW==0) sumW=1;
      h->SetBinContent(ibin, avg);
      h->SetBinError(ibin, 0.);
      //double errSqr=sigma - avg*avg;
      //h->SetBinError  (ibin, sqrt(errSqr));
    }
  }
  HERE("plot");

  // Plot the distributions
  
  // Tower counts
  if (1) {
    CPlot cpTowCount("towCount","Tower count","Tower No.","Count");
    for (unsigned int i=0; i<hTowV.size(); ++i) {
      cpTowCount.AddToStack(hTowV[i],cutV[i]->name(),colorV[i]);
    }
    TCanvas *cx=new TCanvas("cTC","cTC",600,600);
    cpTowCount.Draw(cx);
    cx->Update();
  }

  // Balancing
  if (1) {
    HERE("a");
    for (unsigned int i=0; i<hBalanceEtDiffV.size(); ++i) {
      TString canvName=TString("canvBalance_") + hBalanceEtDiffV[i]->GetName();
      TCanvas *ct= new TCanvas(canvName,canvName,600,600);
      AdjustFor2DplotWithHeight(ct);
      hBalanceEtDiffV[i]->Draw("COLZ");
      ct->Update();
    }
  }

  // Ratio vs diffEt
  if (1) {
    HERE("a");
    for (unsigned int i=0; i<hDiffRatioVsDeltaEotherV.size(); ++i) {
      TString canvName=TString("canvRatio_")
	+ hDiffRatioVsDeltaEotherV[i]->GetName();
      TCanvas *ct= new TCanvas(canvName,canvName,600,600);
      AdjustFor2DplotWithHeight(ct);
      hDiffRatioVsDeltaEotherV[i]->Draw("COLZ");
      ct->Update();
    }
  }

  // Balance vs cf
  if (plotBalanceVsCf) {
    for (unsigned int i=0; i<hTowBalance_iTow_V.size(); ++i) {
      TString canvName= Form("canvBalance_iTow%d",int(i)-MAXIETA);
      TCanvas *ct=new TCanvas(canvName,canvName,600,600);
      //printHisto(hTowBalance_iTow_V[i]);
      hTowBalance_iTow_V[i]->Draw("LPE");
      ct->Update();
    }
  }

  if (1) {
    TString canvName="canv_tag_pfE_over_genE";
    TCanvas *ctag=new TCanvas(canvName,canvName,600,600);
    h1TagJet_pfE_over_genE->Draw("hist");
    ctag->Update();

    canvName="canv_probe_pfE_over_genE";
    TCanvas *cprobe=new TCanvas(canvName,canvName,600,600);
    h1ProbeJet_pfE_over_genE->Draw("hist");
    cprobe->Update();
  }

  if (1) {
    TString canvName="canv_tag_rhE_over_genE";
    TCanvas *ctag=new TCanvas(canvName,canvName,600,600);
    h1TagJet_rhE_over_genE->Draw("hist");
    ctag->Update();

    canvName="canv_probe_rhE_over_genE";
    TCanvas *cprobe=new TCanvas(canvName,canvName,600,600);
    h1ProbeJet_rhE_over_genE->Draw("hist");
    cprobe->Update();
  }
}

// ----------------------------------------------------------------
