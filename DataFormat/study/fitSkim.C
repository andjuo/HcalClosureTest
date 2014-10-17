#include "DijetRespCorrDataExtended.h"
#include "helper.hh"
#include "Analyzer.hh"
#include <TFile.h>
#include <TTree.h>
#include <assert.h>

// ----------------------------------------------------------------

void fitSkim(Long64_t maxEntries=10, int fitStrategy_User=-1) {
  DijetRespCorrDatumExtended_t *d= new DijetRespCorrDatumExtended_t();
  TFile inpFile("skim.root","READ");
  TTree *inpTree = (TTree*)inpFile.Get("dijet_data");
  assert(inpTree);
  inpTree->SetBranchAddress("dijet_data",&d);

  Cuts_t applyCut;
  //applyCut.setEtDiff(0.,10., -1,-1);
  //applyCut.setEtDiff(0.,10., 0, 10);
  //applyCut.setEtDiffRatio(-2.,0.);

  // The third jet is a sum of the remaining jets
  // It's better to apply a cut rather than completely remove
  // the contribution
  int removeThirdJetContribution=0;

  TString hname="hTowCount";
  TH1D *hTowCount= new TH1D(hname,hname, 2*MAXIETA+3, -MAXIETA-1, MAXIETA+2);
  hTowCount->SetDirectory(0);

  DijetRespCorrData data;
  data.SetDoCandTrackEnergyDiff(2); // original balance
  //data.SetDoCandTrackEnergyDiff(3); // MET
  if (fitStrategy_User>0) data.SetDoCandTrackEnergyDiff(fitStrategy_User);
  data.SetParStep(1e-4);
  data.SetParMax(2.);
  data.SetParMin(1e-6);

  Long64_t nEntries= inpTree->GetEntriesFast();
  if (maxEntries<0) maxEntries=nEntries;

  double minWeight0=1e9, maxWeight0=-1e9;
  double thrWeight=1e-5;
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

    if (applyCut.passCuts(*d)) {
      // Passed the selection
      FillTowerCount(hTowCount,*d,-1,1);

      if (removeThirdJetContribution==1) {
	d->SetThirdJetPx(0.);
	d->SetThirdJetPy(0.);
      }

      if (minWeight0>d->GetWeight()) {
	//std::cout << "new minWeight0=" << d->GetWeight() << "\n";
	minWeight0= d->GetWeight();
      }
      if (maxWeight0<d->GetWeight()) {
	//std::cout << "new maxWeight0=" << d->GetWeight() << "\n";
	maxWeight0= d->GetWeight();
      }

      if (d->GetWeight()>thrWeight) {
	data.push_back(*d);
      }
    }
  }
  inpFile.Close();

  std::cout << data.GetSize() << " data entries" << std::endl;
  std::cout << "minWeight0=" << minWeight0 << ", maxWeight0=" << maxWeight0 << "\n";

  printHisto(hTowCount);

  int iniCfArrSize= NUMTOWERS;
  if (data.GetDoCandTrackEnergyDiff()==4) iniCfArrSize++;
  std::vector<double> iniCfVals(iniCfArrSize,1.); // no correction
  std::vector<int> fixParameters;
  double maxWeight= 0.01*hTowCount->GetMaximum();
  fixParameters.reserve(NUMTOWERS);
  if (1)
  for (int ibin=1; ibin<=hTowCount->GetNbinsX(); ++ibin) {
    int iTow= ibin - 2 - MAXIETA;
    std::cout << "iTow=" << iTow << "\n";
    if (hTowCount->GetBinContent(ibin) < 0.05*maxWeight) {
      std::cout << "fixing iTow=" << iTow << "\n";
      fixParameters.push_back(iTow);
    }
  }

  TH1D* hCoef= data.doFit_v2("hCoef","hCoef", iniCfVals, fixParameters);

  //xMinuit->mnscan();
  printHisto(hCoef);

  if (1) {
    TCanvas *cx=new TCanvas("cC","cC",600,600);
    hCoef->SetMarkerStyle(24);
    hCoef->Draw("LPE");
    cx->Update();
  }

  /*
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
  */
}

// ----------------------------------------------------------------
