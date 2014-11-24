// Small macro to check the influence of the bug

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TLine.h>
#include "MitStyleRemix.hh"
#include <iostream>

const float cPi= 4*atan(1);

// ------------------------------------------------------

inline
double calc_dPhi(double phi1, double phi2, int version) {
  double dphi=fabs(phi1-phi2);
  //std::cout << "phi1=" << phi1 << ", phi2=" << phi2 << "\n";
  if (version==0) {
    // earlier version
    while (fabs(dphi-cPi)>cPi) {
      //std::cout << "ver0: dphi=" << dphi << "\n";
      dphi = (2*cPi - dphi);
    }
  }
  else {
    // corrected version
    while (dphi>cPi) {
      dphi = fabs(2*cPi - dphi);
      //std::cout << "ver1: dphi=" << dphi << "\n";
    }
    //std::cout << "ver1 final: dphi=" << dphi << "\n";
  }
  return dphi;
}

// ------------------------------------------------------

void mapDPhi() {
  gStyle->SetPalette(1);

  double phiMin= -cPi;
  double phiMax=  cPi;
  int phiCount=13;
  double phiStep=(phiMax-phiMin)/double(phiCount);
  int cCorr=0;
  std::cout << "phiStep=" << phiStep << ", phiMin+0.5phiStep=" << (phiMin+phiStep/2) << "\n";

  TH2D *h2_dPhi_ver0= new TH2D("h2_dPhi_ver0",
			       "#Delta#phi ver.0 (bugged);#phi_{1};#phi_{2}",
			       phiCount+cCorr, phiMin,phiMax,
			       phiCount+cCorr, phiMin,phiMax);
  TH2D *h2_dPhi_ver1= new TH2D("h2_dPhi_ver1",
			       "#Delta#phi ver.1;#phi_{1};#phi_{2}",
			       phiCount+cCorr, phiMin,phiMax,
			       phiCount+cCorr, phiMin,phiMax);

  TAxis *ax= h2_dPhi_ver0->GetXaxis();
  std::cout << "histo 1st bin: " << ax->GetBinLowEdge(1) << " " << (ax->GetBinLowEdge(1)+ax->GetBinWidth(1)) << "\n";


  h2_dPhi_ver0->SetDirectory(0);
  h2_dPhi_ver1->SetDirectory(0);
  h2_dPhi_ver0->SetStats(0);
  h2_dPhi_ver1->SetStats(0);

  double phi1=0, phi2=0;
  for (int iPhi1=0; iPhi1<=phiCount; ++iPhi1) {
    phi1= phiMin + (iPhi1+0.5)*phiStep;
    std::cout << "phi1=" << phi1 << "\n";
    for (int iPhi2=0; iPhi2<=phiCount; ++iPhi2) {
      phi2= phiMin + (iPhi2+0.5)*phiStep;
      if (0 && (iPhi1==iPhi2)) {
	std::cout << "iphi1=iphi2=" << iPhi1 << ", dPhi=" << calc_dPhi(phi1,phi2,0) << " " << calc_dPhi(phi1,phi2,1) << "\n";
      }
      h2_dPhi_ver0->Fill(phi1,phi2, calc_dPhi(phi1,phi2,0));
      h2_dPhi_ver1->Fill(phi1,phi2, calc_dPhi(phi1,phi2,1));
    }
  }

  TCanvas *cx= new TCanvas("cx","cx",1200,600);
  cx->Divide(2,1);
  AdjustFor2DplotWithHeight(cx);

  h2_dPhi_ver0->GetZaxis()->SetRangeUser(-1e-6,2*cPi);
  h2_dPhi_ver1->GetZaxis()->SetRangeUser(-1e-6,2*cPi);

  cx->cd(1);
  h2_dPhi_ver0->Draw("COLZ");
  TLine *lineA1= new TLine(phiMin,0, 0, cPi);
  TLine *lineA2= new TLine(0,-cPi, phiMax,0);
  lineA1->Draw();
  lineA2->Draw();
  cx->cd(2);
  h2_dPhi_ver1->Draw("COLZ");
  TLine *lineB1= new TLine(phiMin,0, 0, cPi);
  TLine *lineB2= new TLine(0,-cPi, phiMax,0);
  lineB1->Draw();
  lineB2->Draw();
  cx->Update();
}
