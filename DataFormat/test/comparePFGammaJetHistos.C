#define __localRun
#include "../interface/HistoCollector.h"
#include "CPlot.hh"
#include "ComparisonPlot.hh"

// -----------------------------------------------------------------

void CompareHistos(TH1D *h1, TH1D *h2,
		   double weight1, double weight2,
		   TString label1, TString label2,
		   int normalize, int stack,
		   TString drawOpt="LPE");


// -----------------------------------------------------------------

void comparePFGammaJetHistos(TString fname1, TString fname2,
			     double weight1=1., double weight2=1.,
			     TString label1="file1", TString label2="file2",
			     int compareFlag= 1 /*(1|2|4)*/,
			     int normalize=0 )
{
  HistoCollector_t hc1, hc2;

  int loadTH1D = ((compareFlag & 1) != 0) ? 1:0;
  int loadTH2D = ((compareFlag & 2) != 0) ? 1:0;
  int loadCanvas= ((compareFlag & 4) != 0) ? 1:0;

  // load all histos
  if (!hc1.LoadFromFile(fname1,loadTH1D,loadTH2D,loadCanvas,"file1_")) {
    std::cout << "failed to load from file <" << fname1 << ">\n";
    return;
  }
  if (!hc2.LoadFromFile(fname2,loadTH1D,loadTH2D,loadCanvas)) {
    std::cout << "failed to load from file <" << fname2 << ">\n";
    return;
  }

  if (loadTH1D) {
    // comapre several histograms
    for (unsigned int i1=0; i1<hc1.CountTH1D(); i1++) {
      //if (i1>2) break;
      TH1D *h1= hc1.GetH1D(i1);
      TString searchName=h1->GetName();
      searchName.Remove(0,TString("file1_").Length());
      std::cout << "i1=" << i1 << ", searchName=" << searchName << "\n";
      // try to locate the histo in the second set
      TH1D *h2= hc2.GetH1D(searchName);

      int stack=0;
      CompareHistos(h1,h2,weight1,weight2,label1,label2,normalize,stack);
    }
  }

  return;
}

// -----------------------------------------------------------------

void CompareHistos(TH1D *h1_inp, TH1D *h2_inp,
		   double weight1, double weight2,
		   TString label1, TString label2,
		   int normalize, int stack,
		   TString drawOpt)
{

  if (!h1_inp) {
    std::cout << "CompareHistos1D: h1_inp=NULL is not allowed\n";
    return;
  }

  TH1D *h1= (TH1D*)h1_inp->Clone(h1_inp->GetName() + TString("_tmp"));
  h1->SetDirectory(0);
  h1->Sumw2();
  TH1D *h2= NULL;
  if (h2_inp) {
    h2= (TH1D*)h2_inp->Clone(h2_inp->GetName() + TString("_tmp"));
    h2->SetDirectory(0);
    h2->Sumw2();
  }

  if (normalize) {
    // 1) draw the histogram in a hidden mode to determine the x range
    // 2) normalize histograms
    bool isBatch= gROOT->IsBatch();
    gROOT->SetBatch(true);
    TCanvas *cTmp= new TCanvas("cTmp","cTmp",600,600);
    h1->Draw("LPE");
    cTmp->Update();
    Double_t xmin= cTmp->GetUxmin();
    Double_t xmax= cTmp->GetUxmax();
    int binMin= h1->FindBin(xmin);
    int binMax= h1->FindBin(xmax);
    std::cout << "binMin=" << binMin << ", binMax=" << binMax << "\n";
    double norm1= h1->Integral(binMin,binMax);
    if (norm1!=double(0)) h1->Scale(1/norm1);
    if (h2) {
      double norm2= h2->Integral(binMin,binMax);
      if (norm2!=double(0)) h2->Scale(1/norm2);
    }
    delete cTmp;
    gROOT->SetBatch(isBatch);
  }

  if (h1 && (weight1!=double(1))) h1->Scale(weight1);
  if (h2 && (weight2!=double(1))) h2->Scale(weight2);

  TString canvName= TString("canv_") + h1_inp->GetName();
  canvName.ReplaceAll("h1_","");
  TCanvas *cx= new TCanvas(canvName,canvName,700,700+(1-normalize)*100);

  if (stack) {
    // Plot a stack plot
    CPlot cp("cp","",h1->GetXaxis()->GetTitle(),h1->GetYaxis()->GetTitle());
    cp.AddToStack(h1,label1,kBlack);
    if (h2) cp.AddToStack(h2,label2,kBlue);
    cp.Draw(cx);
  }
  else {
    ComparisonPlot_t cp("cp","",h1->GetXaxis()->GetTitle(),
			h1->GetYaxis()->GetTitle(),"ratio");
    cp.Prepare2Pads(cx);
    cp.AddHist1D(h1,label1,drawOpt,kBlack);
    if (h2) cp.AddHist1D(h2,label2,drawOpt,kBlue);
    cp.Draw(cx);
  }
  cx->Update();
  return;
}

// -----------------------------------------------------------------
