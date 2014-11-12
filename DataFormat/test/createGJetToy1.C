#include "link_GammaJetFit.h"
#include <TRandom3.h>

void createGJetToy1(int nEvents=10, int model=1) {
  double cPi= 4*atan(1);
  double ETgammaTrue=100.;
  double ETgammaRes=1.;

  double gammaPhiRes= 0.1*cPi;
  double gammaEtaRes= 0.1*cPi;

  double jetPhiRes= 0.3*cPi;
  double jetEtaRes= 0.3*cPi;
  double jetETRes = 0.08*ETgammaTrue;
  if (model==4) {
    jetETRes=0.5*ETgammaTrue;
    jetEtaRes=0.5*cPi;
  }
  double jetEcalFrac= 0.3;
  double jetEcalFracRes= 0.1;

  GammaJetEvent_t::Class()->IgnoreTObjectStreamer();
  GammaJetEventAuxInfo_t::Class()->IgnoreTObjectStreamer();

  GammaJetEvent_t *dt= new GammaJetEvent_t();
  GammaJetEventAuxInfo_t aux;

  TFile fout(Form("gjet_toy1_model%d.root",model),"recreate");
  TTree *tree= new TTree("gjet_data","gjet_data");
  tree->Branch("gjet_data","GammaJetEvent_t",&dt);

  for (int iEv=0; iEv<nEvents; ++iEv) {
    aux.SetEventNo(iEv);
    dt->SetAuxInfo(aux);

    double gammaPhiTrue= gRandom->Uniform(2*cPi);
    double gammaPhi = gRandom->Gaus(gammaPhiTrue,gammaPhiRes);
    double gammaEtaTrue= gRandom->Uniform(2*cPi) - cPi;
    double gammaEta = gRandom->Gaus(gammaEtaTrue,gammaEtaRes);
    double gammaEt = gRandom->Gaus(ETgammaTrue,ETgammaRes);

    double gammaE= gammaEt * cosh(gammaEta);
    dt->SetTagEEtaPhi(gammaE,gammaEta,gammaPhi);

    double jetPhiTrue= 2*cPi - gammaPhiTrue;
    double jetPhi= gRandom->Gaus(jetPhiTrue,jetPhiRes);
    double jetEtaTrue = gRandom->Uniform(2*cPi) - cPi;
    double jetEta = gRandom->Gaus(jetEtaTrue,jetEtaRes);
    double jetETtot = gRandom->Gaus(ETgammaTrue,jetETRes);
    double jetE= jetETtot * cosh(jetEta);

    std::map<Int_t,Double_t> jetHCal;
    double jetEcalE=0;
    double jetHcalEtot=0;

    switch(model) {
    case 1: jetEcalE=jetE; break;
    case 2: {
      jetHCal[0+MAXIETA]= jetHcalEtot;
    }
      break;
    case 3:
    case 4: {
      jetEcalE= jetE* gRandom->Gaus(jetEcalFrac,jetEcalFracRes);
      jetHcalEtot= jetE - jetEcalE;
      const int distrSize=10;
      int idxCenter= int(gRandom->Uniform(distrSize) - 10);
      double en=jetHcalEtot;
      std::cout << "center=" << idxCenter << ", jetHcalEtot=" << jetHcalEtot << "\n";
      while (en>1e-3) {
	int idx= int(gRandom->Gaus(idxCenter,0.75*distrSize));
	if (idx>MAXIETA) idx=MAXIETA;
	else if (idx<-MAXIETA) idx=-MAXIETA;
	double w= fabs(0.1*gRandom->Gaus(1.,0.5));
	if (w>1.) w=1.;
	double portion= w*en;
	if ((portion>en) || (en<1)) portion=en;
	jetHCal[idx] += portion;
	//std::cout << "distribute energy " << en << " put into " << idx << " dE=" << portion << "\n";
	en-= portion;
      }
    }
      break;
    default:
      std::cout << "not ready for model=" << model << "\n";
      return;
    }
    dt->SetProbeEtaPhiEn(jetEta,jetPhi,jetEcalE,jetHCal);
    tree->Fill();
  }

  tree->Write();
  fout.Close();
  std::cout << "file <" << fout.GetName() << "> created\n";
  std::cout << "nEvents=" << nEvents << "\n";
}
