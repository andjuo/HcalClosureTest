#include "dijet_PFNtuple.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>


void plotThirdJetPt(int imax=100) {
  dijet_PFNtuple inpData;

  // deactivate all branches, except for pTs of the jets
  if (0) {
    inpData.fChain->SetBranchStatus("*",0);
    inpData.fChain->SetBranchStatus("tpfjet_pt",1);
    inpData.fChain->SetBranchStatus("ppfjet_pt",1);
    inpData.fChain->SetBranchStatus("pf_thirdjet_px",1);
    inpData.fChain->SetBranchStatus("pf_thirdjet_py",1);
  }
  else {
    inpData.DeactivateBranches();
    inpData.ActivateBranches(5,"tpfjet_pt","ppfjet_pt","pf_thirdjet_px","pf_thirdjet_py","pf_weight");
  }

  TH1D *hJetTagPt= new TH1D("hJetTagPt","pT of the tag jet;p_{T} [GeV];count", 200,0.,1000.);
  hJetTagPt->SetDirectory(0);
  hJetTagPt->SetLineColor(kRed);
  TH1D *hJetProbePt= new TH1D("hJetProbePt","pT of the probe jet;p_{T} [GeV];count", 200,0.,1000.);
  hJetProbePt->SetDirectory(0);
  hJetProbePt->SetLineColor(kBlue);
  TH1D *hJet3Pt= new TH1D("hJet3Pt","pT of third jet;p_{T} [GeV];count", 100,0.,100.);
  hJet3Pt->SetDirectory(0);
  TH1D *hJetPtAve= new TH1D("hJetPtAve","average pT of 2 leading jets;0.5*(p_{T,j1}+p_{T,j2}) [GeV];count", 200,0.,1000.);
  hJetPtAve->SetDirectory(0);

  TH1D *hJet3Pt_div_ave= new TH1D("hJet3Pt_div_ave","pT_{j3}*2/(pT_{j1}+pT_{j2});p_{T,j3}/[0.5*(p_{T,j1}+p_{T,j2})];count", 50,0.,5.);
  hJet3Pt_div_ave->SetDirectory(0);

  TH2D *h2TagVsProbeJetPt= new TH2D("h2TagVsProbeJetPt","leading jet pT comparison;tag jet p_{T};probe jet p_{T}",100,0.,1000.,100,0.,1000.);
  h2TagVsProbeJetPt->SetDirectory(0);
  h2TagVsProbeJetPt->SetStats(0);

  Long64_t nentries = inpData.fChain->GetEntriesFast();
  if (nentries<=0) { std::cout << "nentries=" << nentries << "\n"; return; }
  if (imax<0) imax=int(nentries);

  Long64_t nbytes = 0, nb = 0;
  int badRatioCount=0;

  for (Long64_t jentry=0; (jentry<nentries) && (jentry<imax); jentry++) {
    Long64_t ientry = inpData.LoadTree(jentry);
    if (ientry < 0) break;
    nb = inpData.GetEntry(jentry);   nbytes += nb;
    if (jentry%100000==0) std::cout << "... reading entry " << jentry << "\n";
    // if (Cut(ientry) < 0) continue;

    double jet3pt= sqrt(pow(inpData.pf_thirdjet_px,2) + pow(inpData.pf_thirdjet_py,2));
    //if (jet3pt>15) continue;
    double jet12ptAve= 0.5*(inpData.tpfjet_pt + inpData.ppfjet_pt);
    double ratio=(jet12ptAve==0) ? -9.99 : jet3pt/jet12ptAve;
    if (ratio<0) badRatioCount++;
    //if (ratio>0.2) continue;

    hJetTagPt->Fill(inpData.tpfjet_pt, inpData.pf_weight);
    hJetProbePt->Fill(inpData.ppfjet_pt, inpData.pf_weight);
    hJet3Pt->Fill(jet3pt, inpData.pf_weight);
    hJetPtAve->Fill(jet12ptAve, inpData.pf_weight);
    hJet3Pt_div_ave->Fill(ratio, inpData.pf_weight);
    h2TagVsProbeJetPt->Fill(inpData.tpfjet_pt, inpData.ppfjet_pt,
			    inpData.pf_weight);
  }

  TCanvas *cx= new TCanvas("cx","pT3",600,600);
  hJet3Pt->Draw("hist");
  cx->Update();

  TCanvas *cy= new TCanvas("cy","pT12ave",600,600);
  hJetPtAve->Draw("hist");
  hJetTagPt->Draw("histsame");
  hJetProbePt->Draw("histsame");
  cy->Update();

  TCanvas *cz= new TCanvas("cz","pT3_div_pT12ave",600,600);
  hJet3Pt_div_ave->Draw("hist");
  cz->Update();

  TCanvas *ct= new TCanvas("ct","tag_vs_probe_jet_pt",600,600);
  h2TagVsProbeJetPt->Draw("COLZ");
  ct->Update();

  if (badRatioCount) {
    std::cout << "there were " << badRatioCount << "bad ratios\n";
  }
  else {
    std::cout << "there were no bad ratios (div 0)\n";
  }
}
