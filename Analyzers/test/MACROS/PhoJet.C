#define PhoJet_cxx
#include "PhoJet.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "JetUtilMC.h"
#include "TVector2.h"

using namespace std;
void PhoJet::Loop()
{
   if (fChain == 0) return;
   TFile* file=new TFile("Hist_gammaJet.root","recreate");
   TH1D* h_response = new TH1D("h_response","h_response",100,0.,5.);
   TH1D* h_phoPt0=new TH1D("h_phoPt0","h_phoPt0",100,0.,300.);
   TH1D* h_jetPt0=new TH1D("h_jetPt0","h_jetPt0",100,0.,300.);
   TH1D* h_response1 = new TH1D("h_response1","h_response1",100,0.,5.);
   TH1D* h_dphi=new TH1D("h_dphi","h_dphi",100,0,3.2);
   TH1D* h_alpha=new TH1D("h_alpha","h_alpha",100,0.,2);
   TH1D* h_response2 = new TH1D("h_response2","h_response2",100,0.,5.);
   TH1D* h_2ndjetPt0=new TH1D("h_2ndjetPt0","h_2ndjetPt0",100,-2.,300.);


   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (ientry%100 == 0 ){
        float perc = float(ientry)*1000./float(nentries);
        printf("Processed Events          %i       out of %i      %.0f %%\n",int(ientry), int (nentries), perc);
      }
      
      //kine.
      h_phoPt0->Fill(tagPho_et);
      h_jetPt0->Fill(ppfjet_pt);
      h_2ndjetPt0->Fill(pf_2ndjet_pt);
      double response = ppfjet_pt/tagPho_et;
      h_response->Fill(response); 


      //selections
      bool PhoSel = 0;
      if(fabs(tagPho_eta)<1.4442 && tagPho_sieie<0.013 && tagPho_HoE<0.05 && 
	 (!tagPho_pixelSeed))PhoSel=1;
      if(ppfjet_pt<15 || tagPho_et<15)continue;
      if(!PhoSel)continue;
      if(PhoSel==1)h_response1->Fill(response); 

      //require back-to-back pair & control extra jet activity
      double ddphi =dPhi(ppfjet_phi,tagPho_phi);
      h_dphi->Fill(ddphi);
      
      double alpha=pf_2ndjet_pt/tagPho_et;
      h_alpha->Fill(alpha);

      if(ddphi>2.95 && alpha<0.20)h_response2->Fill(response);
   }
   file->Write();
   file->Close();
}
