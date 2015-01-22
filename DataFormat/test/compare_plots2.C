#define __localRun
#include "../interface/HistoCollector.h"
#include "CPlot.hh"
#include "ComparisonPlot.hh"
#include "helper.h"

// -----------------------------------------------------------------

typedef enum { _rangeFirst=0,
               _50to80, _80to120, _120to170, _170to300, _300to470, 
               _470to600, _600to800, _rangeLast } TRanges_t;

TString getName(TRanges_t r);

// -----------------------------------------------------------------

/*
void CompareHistos(TH1D *h1, TH1D *h2,
		   double weight1, double weight2,
		   TString label1, TString label2,
		   int normalize, int stack,
		   TString drawOpt="LPE");
*/

// -----------------------------------------------------------------

void DrawHistos(std::vector<TH1D*> &histosV,
		std::vector<TString> &labelsV,
		TString drawOpt="LPE",
		int normalizeHistos=0);


// -----------------------------------------------------------------

void compare_plots2(int the_case=1, int plotIdx=0)
{
  std::vector<HistoCollector_t*> hcV;
  std::vector<TString> fnameV,labelV;
  std::vector<double> weightV;

  if (the_case==1) {
    TString fname1="dir-plots/saved_chkplot_QCDflat-80to120.root";
    TString fname2="dir-plots/saved_chkplot_QCDflat-120to170.root";
    TString fname3="dir-plots/saved_chkplot_QCDflat-170to300.root";
    labelV.push_back( "80to120"); fnameV.push_back(fname1);
    labelV.push_back("120to170"); fnameV.push_back(fname2);
    labelV.push_back("170to300"); fnameV.push_back(fname3);
    for (unsigned int i=0; i<fnameV.size(); i++) weightV.push_back(1.);
  }
  else if ((the_case>=2) && (the_case<=5)) {
    double QCDweight[7]= { 1.35839, 0.172428, 0.026111, 0.00587131, 0.000294313, 2.85065e-05, 6.75332e-06 };
    if ((the_case==3) || (the_case==4) || (the_case==5)) {
      labelV.push_back("QCD flat (Z2star)");
      fnameV.push_back("dir-plots/saved_aodplot-QCD_Z2starFlat.root");
      weightV.push_back(1.);
    }
    if ((the_case==4) || (the_case==5)) {
      //weightV[0]=3000;
      labelV.push_back("QCD flat (DT6)");
      fnameV.push_back("dir-plots/saved_recoplot-QCD_D6T.root");
      //weightV.push_back(3000.);
      weightV.push_back(1.);
    }
    if (the_case!=5)
    for (int r=int(_50to80); r<int(_470to600); r++) {
      TString name=getName(TRanges_t(r));
      labelV.push_back(name);
      fnameV.push_back(Form("dir-plots/saved_aodplot-QCD_pt%s.root",
			    name.Data()));
      weightV.push_back(QCDweight[r-int(_50to80)]);
    }
  }
  else if (the_case==6) {
    labelV.push_back("QCD flat (Z2*,AOD)");
    fnameV.push_back("dir-plots/saved_aodplot-QCD_Z2starFlat.root");
    weightV.push_back(1.);

    TString fname1="dir-plots/saved_recoplot-QCD_Z2starReco_pt80to120.root";
    TString fname2="dir-plots/saved_recoplot-QCD_Z2starReco_pt120to170.root";
    TString fname3="dir-plots/saved_recoplot-QCD_Z2starReco_pt170to300.root";
    fnameV.push_back(fname1); labelV.push_back("80to120 (Z2*,RECO)");
    weightV.push_back(0.172428);
    fnameV.push_back(fname2); labelV.push_back("120to170 (Z2*,RECO)");
    weightV.push_back(0.026111 * 5985732./5443023.);
    fnameV.push_back(fname3); labelV.push_back("170to300 (Z2*,RECO)");
    weightV.push_back(0.00587131);
  }
  //else if (the_case==2) {
  //  fname1="dir-plots/saved_chkplot_relValQCD-80to120.root";
  //  fname2="dir-plots/saved_chkplot_relValQCDflat-80to120.root";
  //}
  else {
    std::cout << "not ready for the_case=" << the_case << "\n";
    return;
  }

  if (labelV.size()!=fnameV.size()) {
    std::cout << "sizes: labelV " << labelV.size()
	      << ", fnameV " << fnameV.size() << "\n";
    return;
  }


  int loadTH1D = 1; //((compareFlag & 1) != 0) ? 1:0;
  int loadTH2D = 0; //((compareFlag & 2) != 0) ? 1:0;
  int loadCanvas= 0; //((compareFlag & 4) != 0) ? 1:0;

  // load all histos
  hcV.reserve(fnameV.size());  
  for (unsigned int i=0; i<fnameV.size(); ++i) {
    HistoCollector_t *hc= new HistoCollector_t();
    if (!hc->LoadFromFile(fnameV[i],loadTH1D,loadTH2D,loadCanvas,
			  Form("file%d_",i+1))) {
      std::cout << "failed to load from file <" << fnameV[i] << ">\n";
      return;
    }
    hcV.push_back(hc);
  }

  TString searchHisto="h1_evPtHat";
  switch(plotIdx) {
  case 0: break;
  case 1: searchHisto= "h1_jet_pTPF"; break;
  case 2: searchHisto="h1_jet_pTGen"; break;
  case 3: searchHisto="h1_jet_energy_PFoverGen"; break;
  case 4: searchHisto="h1_photon_pTreco"; break;
  default:
    std::cout << "plotIdx=" << plotIdx << " is not ready\n";
    return;
  }

  TH1D* hFlat=NULL;
  TH1D* hTot=NULL;


  std::vector<TH1D*> histosV;
  histosV.reserve(hcV.size());
  for (unsigned int i=0; i<hcV.size(); ++i) {
    TString do_search= searchHisto;
    if (searchHisto=="h1_evPtHat") {
      if (((the_case==3) || (the_case==6)) && 
	  (i==0)) do_search.Append("Weighted");
      if (((the_case==4) || (the_case==5)) && (i<2)) {
	if (1) do_search.Append("Weighted");
	else do_search.Append("EvWeighted");
      }
    }
    TH1D *h1= hcV[i]->GetH1D(Form("file%d_%s",i+1,do_search.Data()));
    if (!h1) {
      std::cout << "failed to find the histo\n";
      return;
    }
    h1->Scale(weightV[i]);
    histosV.push_back(h1);
    std::cout << "i=" << i << " " << h1->GetTitle()
	      << ". coded extra weight=" << weightV[i] << "\n";

    for (int ibin=1; ibin<=h1->GetNbinsX(); ++ibin) {
	h1->SetBinError(ibin,sqrt(h1->GetBinContent(ibin)));
    }

    if (plotIdx!=3) {
      if ((the_case==3) || (the_case==4)) {
	if (i==0) {
	  hFlat=h1;
	  hTot=(TH1D*)hFlat->Clone("hTot");
	  hTot->SetDirectory(0);
	  hTot->Reset();
	}
	else {
	  if ((the_case!=4) ||
	      ((the_case==4) && (i>1))) {
	    hTot->Add(h1);
	  }
	}
      }
    }
  }

  if (hTot) {
    histosV.push_back(hTot);
    labelV.push_back("hTot");

    TH1D *hR=NULL;
    if (1) {
      hR=(TH1D*)hTot->Clone("hR");
      hR->SetDirectory(0);
      hR->Reset();
    }
    for (int ibin=1; ibin<=hTot->GetNbinsX(); ++ibin) {
      double wtot= hTot->GetBinContent(ibin);
      double wQCD= hFlat->GetBinContent(ibin);
      if (wtot!=0) {
	double pThat= hTot->GetBinLowEdge(ibin) + 0.5*hTot->GetBinWidth(ibin);
	std::cout << "ibin=" << ibin << ", pThat=" << pThat << ", wtot=" << wtot << ", wQCD=" << wQCD << ", ratio=" << wtot/wQCD << "\n";
	std::cout << "log(ratio)=" << log(wtot/wQCD) << "\n";
	if (hR) hR->SetBinContent(ibin,log(wtot));
      }
    }

    if (hR) {
      TCanvas *canvR = new TCanvas("canvR","canvR",600,600);
      hR->Draw("LPE");
      canvR->Update();
      //return;
    }
  }

  int normalizeHistos=0;
  TString drawOpt="LPE";
  if (plotIdx==3) {
    TAxis *ay= histosV[0]->GetYaxis();
    ay->SetTitle(TString(ay->GetTitle()) + TString(" (normalized)"));
    normalizeHistos=1;
    drawOpt="hist";
  }

  DrawHistos(histosV,labelV,drawOpt,normalizeHistos);

  return;
}

// -----------------------------------------------------------------

/*
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
*/

// -----------------------------------------------------------------

void DrawHistos(std::vector<TH1D*> &histosV,
		std::vector<TString> &labelV,
		TString drawOpt,
		int normalizeHistos)
{
  if (histosV.size()!=labelV.size()) {
    std::cout << "histosV.size=" << histosV.size()
	      << ", labelV.size=" << labelV.size() << "\n";
    return;
  }
  if (histosV.size()==0) {
    std::cout << "histosV.size=0\n";
    return;
  }

  const unsigned int colorCount=8;
  const int autoColor[colorCount] = { kBlack, kBlue, kGreen+1, kOrange+1, kRed+1, kViolet, 43, 37 };

  TString canvName= TString("canv_");
  TCanvas *cx= new TCanvas(canvName,canvName,700,700+100);

  TH1D *h1=histosV[0];
  ComparisonPlot_t cp("cp",h1->GetTitle(),h1->GetXaxis()->GetTitle(),
		      h1->GetYaxis()->GetTitle(),"ratio");
  if (!normalizeHistos) {
    cp.SetYRange(1e-1,1e9);
    cp.SetLogy(1);
  }
  else {
    cp.SetYRange(0,0.5);
  }
  cp.Prepare2Pads(cx);
  cp.ErrorsOnRatios(0);
  //cp.SetRatioYRange(0,1.2);
  //cp.SetRatioYRange(0.9,1.1);
 
  for (unsigned int i=0; i<histosV.size(); ++i) {
    TH1D *h1_inp= histosV[i];
    TH1D *h1a= (TH1D*)h1_inp->Clone(h1_inp->GetName() + TString("_tmp"));
    if (i==0) h1a->GetYaxis()->SetNdivisions(510);
    h1a->SetStats(0);
    h1a->SetDirectory(0);
    h1a->Sumw2();
    h1a->SetLineColor(autoColor[i%colorCount]);
    h1a->SetMarkerColor(autoColor[i%colorCount]);

    if (normalizeHistos) {
      h1a->Scale(1/h1a->Integral());
    }

    cp.AddHist1D(h1a,labelV[i],drawOpt,autoColor[i%colorCount]);
  }

  cp.Draw(cx);

  cx->Update();
  return;
}

// -----------------------------------------------------------------

// ----------------------------------------------

TString getName(TRanges_t r) {
  TString name;
  switch(r) {
  case _50to80: name="50to80"; break;
  case _80to120: name="80to120"; break;
  case _120to170: name="120to170"; break;
  case _170to300: name="170to300"; break;
  case _300to470: name="300to470"; break;
  case _470to600: name="470to600"; break;
  case _600to800: name="600to800"; break;
  default: name="unknown";
  }
  return name;
}

// ----------------------------------------------
