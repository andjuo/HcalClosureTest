// plotTreeField.C
// Macro produces histograms for the branches in a specified tree
//
// For large chains macro crashes. Try to use plotTreeFields_v2.C instead.
//
// Known features:
//   ROOT cannot determine the min/max values of the std::vector<>, so
//   some ranges are not very good. These variables could be added manually
//   either in function fillNBinsXMinXMax or getNBinsXMinMax
//
// Jan 27, 2015

#include <TROOT.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TString.h>
#include <TList.h>
#include <TKey.h>
#include <TClass.h>
#include <TLeaf.h>
#include <map>
#include <iostream>
#include <string>

// ----------------------------------------------------------------------

std::map<std::string,int> mapNBins;
std::map<std::string,double> mapMinVal, mapMaxVal;

// ----------------------------------------------------------------------

int fillNBinsXMinXMax(int theSet=0);

int getNBinsXMinMax(std::string branchName,
		    int &nBins, double &xmin, double &xmax);

int adjustRange(double &xmin, double &xmax, int verbose=0);

// ----------------------------------------------------------------------


void plotTreeFields(TString fname, TString treeName, TString outFname,
		    int listBranchesOnly=0,
		    int saveCanvases=0,
		    int keepCanvasesOpen=0)
{
  TChain chain(treeName);
  chain.Add(fname);

  // list files
  if (1) {
    TObjArray* files= chain.GetListOfFiles();
    if (!files) {
      std::cout << "chain.GetListOfFiles is null\n";
      return;
    }

    TIter next(files);
    TChainElement *chEl=0;
    int fileCount=0;
    while ( (chEl=(TChainElement*)next()) ) {
      fileCount++;
      std::cout << "file #" << fileCount << " is " << chEl->GetTitle() << "\n";
      // chEl->GetName would return the tree name
    }
  }

  // activate the chain
  chain.GetEntry(0);
  TTree *ttree= (TTree*)chain.GetTree();
  if (!ttree) std::cout << "ttree is null\n";

  /*
  TString br="tagPho_pt";
  std::cout << "range=" << chain.GetMinimum(br) << " .. " << chain.GetMaximum(br) << "\n";
  return;
  */

  // prepare the ranges
  if (!fillNBinsXMinXMax(1)) {
    std::cout << "failed to prepare the ranges\n";
    return;
  }

  TFile fout(outFname,"RECREATE");
  if (!fout.IsOpen()) {
    std::cout << "failed to create a file <" << fout.GetName() << ">\n";
    return;
  }

  // plot the branches
  TObjArray* branchList= ttree->GetListOfBranches();
  TIter branchIter(branchList);
  TChainElement *brEl=0;
  int brCount=0;
  int nBins=0;
  double xmin=0,xmax=0;
  //TCanvas *cx=new TCanvas("cx","cx",600,600);

  while ( (brEl = (TChainElement*)branchIter()) ) {
    if (!listBranchesOnly) {
      std::cout << "branch name=" << brEl->GetTitle() << "\n";
    }
    TString brNameOrig= brEl->GetName();
    TString brName= brNameOrig;
    int idx=brName.Index("/");
    if (idx>0) brName.Remove(idx,brName.Length());
    if (!listBranchesOnly) std::cout << "work with branch=" << brName << "\n";
    std::cout << brName << "\n";
    brCount++;
    if (listBranchesOnly) continue;
    //if ((brCount>40) && !listBranchesOnly) break;

    nBins=100; xmin=0; xmax=0;
    if (!getNBinsXMinMax(brName.Data(),nBins,xmin,xmax)) {
      // if not listed in the table, try to use auto-range
      nBins=100;
      xmin=chain.GetMinimum(brNameOrig);
      xmax=chain.GetMaximum(brNameOrig);
      adjustRange(xmin,xmax);
    }
    TString canvName= Form("cv_%s",brName.Data());
    TCanvas *cx=new TCanvas(canvName,canvName,600,600);
    std::cout << "ranges: " << nBins << " " << xmin << " .. " << xmax << "\n";
    TString histoName="h1D_" + brName;
    TH1D *h1= new TH1D(histoName,brName,nBins,xmin,xmax);
    TString command= brName + TString(">>") + histoName;
    //ttree->Draw(command); // will draw only currect tree
    chain.Draw(command);
    h1->SetDirectory(0);
    //h1->DrawClone();
    cx->Update();

    fout.cd();
    h1->Write(brName);
    if (saveCanvases) cx->Write();

    if (!keepCanvasesOpen) delete cx;
  }

  return;
}



// ----------------------------------------------------------------------

int getNBinsXMinMax(std::string branchName,
		    int &nBins, double &xmin, double &xmax)
{
  nBins=mapNBins[branchName]; xmin=0; xmax=1;
  if (nBins==0) {
    if ((branchName.find("_px")!=std::string::npos) ||
	(branchName.find("_py")!=std::string::npos) ||
	(branchName.find("_pz")!=std::string::npos)) {
      nBins=100; xmin=-1000; xmax=1000;
      return 1;
    }
    std::cout << "could not locate branch <" << branchName << ">\n";
    return 0;
  }
  xmin= mapMinVal[branchName];
  xmax= mapMaxVal[branchName];
  return 1;
}

// ----------------------------------------------------------------------

int fillNBinsXMinXMax(int theSet)
{
  std::string s;
  if (theSet==1) {
    std::cout << "here\n";
    mapNBins[s]= 10; mapMinVal[s]= 0; mapMaxVal[s]= 1e3;
    s= "photonTrig_fired";
    mapNBins[s]= 2;  mapMinVal[s]= 0; mapMaxVal[s]= 2;
    s= "photonTrig_prescale";
    mapNBins[s]= 10; mapMinVal[s]= 0; mapMaxVal[s]= 1e3;
    s= "jetTrig_fired";
    mapNBins[s]= 2;  mapMinVal[s]= 0; mapMaxVal[s]= 2;
    s= "jetTrig_prescale";
    mapNBins[s]= 10; mapMinVal[s]= 0; mapMaxVal[s]= 1e3;
    s= "tagPho_idLoose";
    mapNBins[s]= 2; mapMinVal[s]= 0; mapMaxVal[s]= 2;
    s= "tagPho_idTight";
    mapNBins[s]= 2; mapMinVal[s]= 0; mapMaxVal[s]= 2;
    s= "tagPho_pixelSeed";
    mapNBins[s]= 2; mapMinVal[s]= 0; mapMaxVal[s]= 2;
  }
  else {
    std::cout << "fillNBinsXMinXMax: theSet=" << theSet << " is not ready\n";
  }
  return 1;
}

// ----------------------------------------------------------------------

int adjustRange(double &xmin, double &xmax, int verbose) {
  double n1= trunc(log(fabs(xmin))/log(10.));
  double n2= trunc(log(fabs(xmax))/log(10.));
  double y1= trunc(xmin/pow(10,n1));
  double y2= trunc(xmax/pow(10,n2));

  if (xmin==0) { n1=0; y1=0; }
  if (xmax==0) { n2=0; y2=0; }

  if (verbose) {
    std::cout << "xmin=" << xmin << ", xmax=" << xmax << "\n";
    std::cout << " n: " << n1 << " " << n2
	      << ", y: " << y1 << " " << y2 << "\n";
  }
  if (n1==0) {
    if (verbose) std::cout << "case n1=0\n";
    xmin= (y1==0) ? 0 : (y1-1)*pow(10,n1);
    xmax= (y2==0) ? 0 : (y2+1)*pow(10,n2);
  }
  else if (n2==0) {
    if (verbose) std::cout << "case n2=0\n";
    if (n1==0) xmin= (xmin>0) ? 10. : -10;
    xmax= (xmax>0) ? 10 : 0;
    if (xmax==xmin) xmax+=1;
  }
  else if (n1<0) {
    // not checked
    if (verbose) std::cout << "case n1<0\n";
    xmin= (xmin<0) ? ((y1-1)*pow(10,n1)) : 0;
    xmax= (y2+1)*pow(10,n2);
  }
  else if ((n1>0) && (n2>0)) {
    if (verbose) std::cout << "case (n1>0) && (n2>0)\n";
    if ((xmin>0) && (xmax>=0)) {
      if (verbose) std::cout << "subcase >0\n";
      xmin=pow(10,n1-1);
      xmax=(y2+1)*pow(10,n2);
    }
    else if ((xmin<0) && (xmax>=0)) {
      if (verbose) std::cout << "subcase <0, >0\n";
      xmin=(y1-1)*pow(10,n1);
      xmax=(y2+1)*pow(10,n2);
    }
    else if ((xmin<0) && (xmax<0)) {
      if (verbose) std::cout << "subcase <0\n";
      xmin=(y1-1)*pow(10,n1);
      xmax=(y2+1)*pow(10,n2);
    }
  }

  if (xmin==xmax) xmax+=1;
  if (verbose) std::cout << " new : " << xmin << " " << xmax << "\n";

  return 1;
}

// ----------------------------------------------------------------------
