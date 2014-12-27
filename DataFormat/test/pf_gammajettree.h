//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 28 11:21:11 2014 by ROOT version 5.34/09
// from TTree pf_gammajettree/tree for gamma+jet balancing using PFJets
// found on file: /media/ssd/ntuple-data-20141028/PhoJet_tree_Summer12DR53X_GPt170to300.root
//////////////////////////////////////////////////////////

#ifndef pf_gammajettree_h
#define pf_gammajettree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Other classes
#include <iostream>
#include <stdarg.h>

//extern const int MAXIETA = 41;

typedef enum { _notJet=0, _leadJet=1, _subleadJet=2 } TJetId_t;

// Fixed size dimensions of array or collections stored in the TTree if any.

class pf_gammajettree {
public :
  //TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<int>     *photonTrig_fired;
   vector<int>     *photonTrig_prescale;
   vector<int>     *jetTrig_fired;
   vector<int>     *jetTrig_prescale;
   Int_t           RunNumber;
   Int_t           LumiBlock;
   Int_t           EventNumber;
   Float_t         EventWeight;
   Float_t         rho2012;
   Float_t         tagPho_pt;
   Float_t         pho_2nd_pt;
   Float_t         tagPho_energy;
   Float_t         tagPho_eta;
   Float_t         tagPho_phi;
   Float_t         tagPho_sieie;
   Float_t         tagPho_HoE;
   Float_t         tagPho_r9;
   Float_t         tagPho_EcalIsoDR04;
   Float_t         tagPho_HcalIsoDR04;
   Float_t         tagPho_HcalIsoDR0412;
   Float_t         tagPho_TrkIsoHollowDR04;
   Float_t         tagPho_pfiso_myphoton03;
   Float_t         tagPho_pfiso_myneutral03;
   vector<vector<float> > *tagPho_pfiso_mycharged03;
   Int_t           tagPho_pixelSeed;
   Int_t           tagPho_ConvSafeEleVeto;
   Int_t           tagPho_idTight;
   Int_t           tagPho_idLoose;
   Float_t         tagPho_genPt;
   Float_t         tagPho_genEnergy;
   Float_t         tagPho_genEta;
   Float_t         tagPho_genPhi;
   Float_t         tagPho_genDeltaR;
   Int_t           nPhotons;
   Int_t           nGenJets;
   Int_t           nPFJets;
   Float_t         ppfjet_pt;
   Float_t         ppfjet_p;
   Float_t         ppfjet_E;
   Float_t         ppfjet_E_NPVcorr;
   Float_t         ppfjet_area;
   Float_t         ppfjet_eta;
   Float_t         ppfjet_phi;
   Float_t         ppfjet_scale;
   Float_t         ppfjet_NeutralHadronFrac;
   Float_t         ppfjet_NeutralEMFrac;
   Int_t           ppfjet_nConstituents;
   Float_t         ppfjet_ChargedHadronFrac;
   Float_t         ppfjet_ChargedMultiplicity;
   Float_t         ppfjet_ChargedEMFrac;
   Float_t         ppfjet_genpt;
   Float_t         ppfjet_genp;
   Float_t         ppfjet_genE;
   Float_t         ppfjet_gendr;
   Float_t         ppfjet_unkown_E;
   Float_t         ppfjet_electron_E;
   Float_t         ppfjet_muon_E;
   Float_t         ppfjet_photon_E;
   Float_t         ppfjet_unkown_px;
   Float_t         ppfjet_electron_px;
   Float_t         ppfjet_muon_px;
   Float_t         ppfjet_photon_px;
   Float_t         ppfjet_unkown_py;
   Float_t         ppfjet_electron_py;
   Float_t         ppfjet_muon_py;
   Float_t         ppfjet_photon_py;
   Float_t         ppfjet_unkown_pz;
   Float_t         ppfjet_electron_pz;
   Float_t         ppfjet_muon_pz;
   Float_t         ppfjet_photon_pz;
   Float_t         ppfjet_unkown_EcalE;
   Float_t         ppfjet_electron_EcalE;
   Float_t         ppfjet_muon_EcalE;
   Float_t         ppfjet_photon_EcalE;
   Int_t           ppfjet_unkown_n;
   Int_t           ppfjet_electron_n;
   Int_t           ppfjet_muon_n;
   Int_t           ppfjet_photon_n;
   Int_t           ppfjet_had_n;
   vector<float>   *ppfjet_had_E;
   vector<float>   *ppfjet_had_px;
   vector<float>   *ppfjet_had_py;
   vector<float>   *ppfjet_had_pz;
   vector<float>   *ppfjet_had_EcalE;
   vector<float>   *ppfjet_had_rawHcalE;
   vector<float>   *ppfjet_had_emf;
   vector<int>     *ppfjet_had_id;
   vector<int>     *ppfjet_had_candtrackind;
   vector<float>   *ppfjet_had_E_mctruth;
   vector<int>     *ppfjet_had_mcpdgId;
   vector<int>     *ppfjet_had_ntwrs;
   Int_t           ppfjet_ntwrs;
   vector<int>     *ppfjet_twr_ieta;
   vector<int>     *ppfjet_twr_iphi;
   vector<int>     *ppfjet_twr_depth;
   vector<int>     *ppfjet_twr_subdet;
   vector<float>   *ppfjet_twr_hade;
   vector<float>   *ppfjet_twr_frac;
   vector<int>     *ppfjet_twr_candtrackind;
   vector<int>     *ppfjet_twr_hadind;
   vector<int>     *ppfjet_twr_elmttype;
   vector<float>   *ppfjet_twr_dR;
   vector<int>     *ppfjet_twr_clusterind;
   Int_t           ppfjet_cluster_n;
   vector<float>   *ppfjet_cluster_eta;
   vector<float>   *ppfjet_cluster_phi;
   vector<float>   *ppfjet_cluster_dR;
   Int_t           ppfjet_ncandtracks;
   vector<float>   *ppfjet_candtrack_px;
   vector<float>   *ppfjet_candtrack_py;
   vector<float>   *ppfjet_candtrack_pz;
   vector<float>   *ppfjet_candtrack_EcalE;
   Float_t         pfjet2_pt;
   Float_t         pfjet2_p;
   Float_t         pfjet2_E;
   Float_t         pfjet2_E_NPVcorr;
   Float_t         pfjet2_area;
   Float_t         pfjet2_eta;
   Float_t         pfjet2_phi;
   Float_t         pfjet2_scale;
   Float_t         pfjet2_NeutralHadronFrac;
   Float_t         pfjet2_NeutralEMFrac;
   Int_t           pfjet2_nConstituents;
   Float_t         pfjet2_ChargedHadronFrac;
   Float_t         pfjet2_ChargedMultiplicity;
   Float_t         pfjet2_ChargedEMFrac;
   Float_t         pfjet2_genpt;
   Float_t         pfjet2_genp;
   Float_t         pfjet2_genE;
   Float_t         pfjet2_gendr;
   Float_t         pfjet2_unkown_E;
   Float_t         pfjet2_electron_E;
   Float_t         pfjet2_muon_E;
   Float_t         pfjet2_photon_E;
   Float_t         pfjet2_unkown_px;
   Float_t         pfjet2_electron_px;
   Float_t         pfjet2_muon_px;
   Float_t         pfjet2_photon_px;
   Float_t         pfjet2_unkown_py;
   Float_t         pfjet2_electron_py;
   Float_t         pfjet2_muon_py;
   Float_t         pfjet2_photon_py;
   Float_t         pfjet2_unkown_pz;
   Float_t         pfjet2_electron_pz;
   Float_t         pfjet2_muon_pz;
   Float_t         pfjet2_photon_pz;
   Float_t         pfjet2_unkown_EcalE;
   Float_t         pfjet2_electron_EcalE;
   Float_t         pfjet2_muon_EcalE;
   Float_t         pfjet2_photon_EcalE;
   Int_t           pfjet2_unkown_n;
   Int_t           pfjet2_electron_n;
   Int_t           pfjet2_muon_n;
   Int_t           pfjet2_photon_n;
   Int_t           pfjet2_had_n;
   vector<float>   *pfjet2_had_E;
   vector<float>   *pfjet2_had_px;
   vector<float>   *pfjet2_had_py;
   vector<float>   *pfjet2_had_pz;
   vector<float>   *pfjet2_had_EcalE;
   vector<float>   *pfjet2_had_rawHcalE;
   vector<float>   *pfjet2_had_emf;
   vector<int>     *pfjet2_had_id;
   vector<int>     *pfjet2_had_candtrackind;
   vector<float>   *pfjet2_had_E_mctruth;
   vector<int>     *pfjet2_had_mcpdgId;
   vector<int>     *pfjet2_had_ntwrs;
   Int_t           pfjet2_ntwrs;
   vector<int>     *pfjet2_twr_ieta;
   vector<int>     *pfjet2_twr_iphi;
   vector<int>     *pfjet2_twr_depth;
   vector<int>     *pfjet2_twr_subdet;
   vector<float>   *pfjet2_twr_hade;
   vector<float>   *pfjet2_twr_frac;
   vector<int>     *pfjet2_twr_candtrackind;
   vector<int>     *pfjet2_twr_hadind;
   vector<int>     *pfjet2_twr_elmttype;
   vector<float>   *pfjet2_twr_dR;
   vector<int>     *pfjet2_twr_clusterind;
   Int_t           pfjet2_cluster_n;
   vector<float>   *pfjet2_cluster_eta;
   vector<float>   *pfjet2_cluster_phi;
   vector<float>   *pfjet2_cluster_dR;
   Int_t           pfjet2_ncandtracks;
   vector<float>   *pfjet2_candtrack_px;
   vector<float>   *pfjet2_candtrack_py;
   vector<float>   *pfjet2_candtrack_pz;
   vector<float>   *pfjet2_candtrack_EcalE;
   Float_t         pf_thirdjet_et;
   Float_t         pf_thirdjet_pt;
   Float_t         pf_thirdjet_p;
   Float_t         pf_thirdjet_px;
   Float_t         pf_thirdjet_py;
   Float_t         pf_thirdjet_E;
   Float_t         pf_thirdjet_eta;
   Float_t         pf_thirdjet_phi;
   Float_t         pf_thirdjet_scale;
   Float_t         met_value;
   Float_t         met_phi;
   Float_t         met_sumEt;
   Float_t         metType1_value;
   Float_t         metType1_phi;
   Float_t         metType1_sumEt;
   Int_t           pf_NPV;

   // List of branches
   TBranch        *b_photonTrig_fired;   //!
   TBranch        *b_photonTrig_prescale;   //!
   TBranch        *b_jetTrig_fired;   //!
   TBranch        *b_jetTrig_prescale;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_EventWeight;   //!
   TBranch        *b_rho2012;   //!
   TBranch        *b_tagPho_pt;   //!
   TBranch        *b_pho_2nd_pt;   //!
   TBranch        *b_tagPho_energy;   //!
   TBranch        *b_tagPho_eta;   //!
   TBranch        *b_tagPho_phi;   //!
   TBranch        *b_tagPho_sieie;   //!
   TBranch        *b_tagPho_HoE;   //!
   TBranch        *b_tagPho_r9;   //!
   TBranch        *b_tagPho_EcalIsoDR04;   //!
   TBranch        *b_tagPho_HcalIsoDR04;   //!
   TBranch        *b_tagPho_HcalIsoDR0412;   //!
   TBranch        *b_tagPho_TrkIsoHollowDR04;   //!
   TBranch        *b_tagPho_pfiso_myphoton03;   //!
   TBranch        *b_tagPho_pfiso_myneutral03;   //!
   TBranch        *b_tagPho_pfiso_mycharged03;   //!
   TBranch        *b_tagPho_pixelSeed;   //!
   TBranch        *b_tagPho_ConvSafeEleVeto;   //!
   TBranch        *b_tagPho_idTight;   //!
   TBranch        *b_tagPho_idLoose;   //!
   TBranch        *b_tagPho_genPt;   //!
   TBranch        *b_tagPho_genEnergy;   //!
   TBranch        *b_tagPho_genEta;   //!
   TBranch        *b_tagPho_genPhi;   //!
   TBranch        *b_tagPho_genDeltaR;   //!
   TBranch        *b_nPhotons;   //!
   TBranch        *b_nGenJets;   //!
   TBranch        *b_nPFJets;   //!
   TBranch        *b_ppfjet_pt;   //!
   TBranch        *b_ppfjet_p;   //!
   TBranch        *b_ppfjet_E;   //!
   TBranch        *b_ppfjet_E_NPVcorr;   //!
   TBranch        *b_ppfjet_area;   //!
   TBranch        *b_ppfjet_eta;   //!
   TBranch        *b_ppfjet_phi;   //!
   TBranch        *b_ppfjet_scale;   //!
   TBranch        *b_ppfjet_NeutralHadronFrac;   //!
   TBranch        *b_ppfjet_NeutralEMFrac;   //!
   TBranch        *b_ppfjet_nConstituents;   //!
   TBranch        *b_ppfjet_ChargedHadronFrac;   //!
   TBranch        *b_ppfjet_ChargedMultiplicity;   //!
   TBranch        *b_ppfjet_ChargedEMFrac;   //!
   TBranch        *b_ppfjet_genpt;   //!
   TBranch        *b_ppfjet_genp;   //!
   TBranch        *b_ppfjet_genE;   //!
   TBranch        *b_ppfjet_gendr;   //!
   TBranch        *b_ppfjet_unkown_E;   //!
   TBranch        *b_ppfjet_electron_E;   //!
   TBranch        *b_ppfjet_muon_E;   //!
   TBranch        *b_ppfjet_photon_E;   //!
   TBranch        *b_ppfjet_unkown_px;   //!
   TBranch        *b_ppfjet_electron_px;   //!
   TBranch        *b_ppfjet_muon_px;   //!
   TBranch        *b_ppfjet_photon_px;   //!
   TBranch        *b_ppfjet_unkown_py;   //!
   TBranch        *b_ppfjet_electron_py;   //!
   TBranch        *b_ppfjet_muon_py;   //!
   TBranch        *b_ppfjet_photon_py;   //!
   TBranch        *b_ppfjet_unkown_pz;   //!
   TBranch        *b_ppfjet_electron_pz;   //!
   TBranch        *b_ppfjet_muon_pz;   //!
   TBranch        *b_ppfjet_photon_pz;   //!
   TBranch        *b_ppfjet_unkown_EcalE;   //!
   TBranch        *b_ppfjet_electron_EcalE;   //!
   TBranch        *b_ppfjet_muon_EcalE;   //!
   TBranch        *b_ppfjet_photon_EcalE;   //!
   TBranch        *b_ppfjet_unkown_n;   //!
   TBranch        *b_ppfjet_electron_n;   //!
   TBranch        *b_ppfjet_muon_n;   //!
   TBranch        *b_ppfjet_photon_n;   //!
   TBranch        *b_ppfjet_had_n;   //!
   TBranch        *b_ppfjet_had_E;   //!
   TBranch        *b_ppfjet_had_px;   //!
   TBranch        *b_ppfjet_had_py;   //!
   TBranch        *b_ppfjet_had_pz;   //!
   TBranch        *b_ppfjet_had_EcalE;   //!
   TBranch        *b_ppfjet_had_rawHcalE;   //!
   TBranch        *b_ppfjet_had_emf;   //!
   TBranch        *b_ppfjet_had_id;   //!
   TBranch        *b_ppfjet_had_candtrackind;   //!
   TBranch        *b_ppfjet_had_E_mctruth;   //!
   TBranch        *b_ppfjet_had_mcpdgId;   //!
   TBranch        *b_ppfjet_had_ntwrs;   //!
   TBranch        *b_ppfjet_ntwrs;   //!
   TBranch        *b_ppfjet_twr_ieta;   //!
   TBranch        *b_ppfjet_twr_iphi;   //!
   TBranch        *b_ppfjet_twr_depth;   //!
   TBranch        *b_ppfjet_twr_subdet;   //!
   TBranch        *b_ppfjet_twr_hade;   //!
   TBranch        *b_ppfjet_twr_frac;   //!
   TBranch        *b_ppfjet_twr_candtrackind;   //!
   TBranch        *b_ppfjet_twr_hadind;   //!
   TBranch        *b_ppfjet_twr_elmttype;   //!
   TBranch        *b_ppfjet_twr_dR;   //!
   TBranch        *b_ppfjet_twr_clusterind;   //!
   TBranch        *b_ppfjet_cluster_n;   //!
   TBranch        *b_ppfjet_cluster_eta;   //!
   TBranch        *b_ppfjet_cluster_phi;   //!
   TBranch        *b_ppfjet_cluster_dR;   //!
   TBranch        *b_ppfjet_ncandtracks;   //!
   TBranch        *b_ppfjet_candtrack_px;   //!
   TBranch        *b_ppfjet_candtrack_py;   //!
   TBranch        *b_ppfjet_candtrack_pz;   //!
   TBranch        *b_ppfjet_candtrack_EcalE;   //!
   TBranch        *b_pfjet2_pt;   //!
   TBranch        *b_pfjet2_p;   //!
   TBranch        *b_pfjet2_E;   //!
   TBranch        *b_pfjet2_E_NPVcorr;   //!
   TBranch        *b_pfjet2_area;   //!
   TBranch        *b_pfjet2_eta;   //!
   TBranch        *b_pfjet2_phi;   //!
   TBranch        *b_pfjet2_scale;   //!
   TBranch        *b_pfjet2_NeutralHadronFrac;   //!
   TBranch        *b_pfjet2_NeutralEMFrac;   //!
   TBranch        *b_pfjet2_nConstituents;   //!
   TBranch        *b_pfjet2_ChargedHadronFrac;   //!
   TBranch        *b_pfjet2_ChargedMultiplicity;   //!
   TBranch        *b_pfjet2_ChargedEMFrac;   //!
   TBranch        *b_pfjet2_genpt;   //!
   TBranch        *b_pfjet2_genp;   //!
   TBranch        *b_pfjet2_genE;   //!
   TBranch        *b_pfjet2_gendr;   //!
   TBranch        *b_pfjet2_unkown_E;   //!
   TBranch        *b_pfjet2_electron_E;   //!
   TBranch        *b_pfjet2_muon_E;   //!
   TBranch        *b_pfjet2_photon_E;   //!
   TBranch        *b_pfjet2_unkown_px;   //!
   TBranch        *b_pfjet2_electron_px;   //!
   TBranch        *b_pfjet2_muon_px;   //!
   TBranch        *b_pfjet2_photon_px;   //!
   TBranch        *b_pfjet2_unkown_py;   //!
   TBranch        *b_pfjet2_electron_py;   //!
   TBranch        *b_pfjet2_muon_py;   //!
   TBranch        *b_pfjet2_photon_py;   //!
   TBranch        *b_pfjet2_unkown_pz;   //!
   TBranch        *b_pfjet2_electron_pz;   //!
   TBranch        *b_pfjet2_muon_pz;   //!
   TBranch        *b_pfjet2_photon_pz;   //!
   TBranch        *b_pfjet2_unkown_EcalE;   //!
   TBranch        *b_pfjet2_electron_EcalE;   //!
   TBranch        *b_pfjet2_muon_EcalE;   //!
   TBranch        *b_pfjet2_photon_EcalE;   //!
   TBranch        *b_pfjet2_unkown_n;   //!
   TBranch        *b_pfjet2_electron_n;   //!
   TBranch        *b_pfjet2_muon_n;   //!
   TBranch        *b_pfjet2_photon_n;   //!
   TBranch        *b_pfjet2_had_n;   //!
   TBranch        *b_pfjet2_had_E;   //!
   TBranch        *b_pfjet2_had_px;   //!
   TBranch        *b_pfjet2_had_py;   //!
   TBranch        *b_pfjet2_had_pz;   //!
   TBranch        *b_pfjet2_had_EcalE;   //!
   TBranch        *b_pfjet2_had_rawHcalE;   //!
   TBranch        *b_pfjet2_had_emf;   //!
   TBranch        *b_pfjet2_had_id;   //!
   TBranch        *b_pfjet2_had_candtrackind;   //!
   TBranch        *b_pfjet2_had_E_mctruth;   //!
   TBranch        *b_pfjet2_had_mcpdgId;   //!
   TBranch        *b_pfjet2_had_ntwrs;   //!
   TBranch        *b_pfjet2_ntwrs;   //!
   TBranch        *b_pfjet2_twr_ieta;   //!
   TBranch        *b_pfjet2_twr_iphi;   //!
   TBranch        *b_pfjet2_twr_depth;   //!
   TBranch        *b_pfjet2_twr_subdet;   //!
   TBranch        *b_pfjet2_twr_hade;   //!
   TBranch        *b_pfjet2_twr_frac;   //!
   TBranch        *b_pfjet2_twr_candtrackind;   //!
   TBranch        *b_pfjet2_twr_hadind;   //!
   TBranch        *b_pfjet2_twr_elmttype;   //!
   TBranch        *b_pfjet2_twr_dR;   //!
   TBranch        *b_pfjet2_twr_clusterind;   //!
   TBranch        *b_pfjet2_cluster_n;   //!
   TBranch        *b_pfjet2_cluster_eta;   //!
   TBranch        *b_pfjet2_cluster_phi;   //!
   TBranch        *b_pfjet2_cluster_dR;   //!
   TBranch        *b_pfjet2_ncandtracks;   //!
   TBranch        *b_pfjet2_candtrack_px;   //!
   TBranch        *b_pfjet2_candtrack_py;   //!
   TBranch        *b_pfjet2_candtrack_pz;   //!
   TBranch        *b_pfjet2_candtrack_EcalE;   //!
   TBranch        *b_pf_thirdjet_et;   //!
   TBranch        *b_pf_thirdjet_pt;   //!
   TBranch        *b_pf_thirdjet_p;   //!
   TBranch        *b_pf_thirdjet_px;   //!
   TBranch        *b_pf_thirdjet_py;   //!
   TBranch        *b_pf_thirdjet_E;   //!
   TBranch        *b_pf_thirdjet_eta;   //!
   TBranch        *b_pf_thirdjet_phi;   //!
   TBranch        *b_pf_thirdjet_scale;   //!
   TBranch        *b_met_value;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_sumEt;   //!
   TBranch        *b_metType1_value;   //!
   TBranch        *b_metType1_phi;   //!
   TBranch        *b_metType1_sumEt;   //!
   TBranch        *b_pf_NPV;   //!

   pf_gammajettree(TString fname="");
   virtual ~pf_gammajettree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual int      Init(const TString &fname);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   // Added methods
   void DeactivateBranches();
   void ActivateBranches(int count, ...); // list n branch names
   void ActivateBranches(const std::vector<TString> &brV); // list of branch names
   void ActivateBranches_forRecHitsEnergyCalc();
   void ActivateBranches_forFitSkim();
   void ActivateBranches_genBasicSet();
   void ActivateBranches_jetID(int leadingJet=1);

   // Added tests
   int passTightJetID(int leadingJet=1) const;
   int passPhotonJetRequirements(int theSet=2) const;

   int passCuts(int leadingJet=1, int theSet=2) const {
     return (passTightJetID(leadingJet) &&
	     passPhotonJetRequirements(theSet)) ? 1:0;
   }

   friend
     std::ostream& operator<<(std::ostream &out, pf_gammajettree &obj) {
     if (out==std::cout) obj.Show();
     else out << "cannot print pf_gammajettree\n";
     return out;
   }

   void PrintSelectedFields(int selection=0);
   // numeric jet identifiers
   Double_t getSumEcalE(int tag, int includeOthers=1) const;
   Double_t getSumHcalE_trackDiffEcal(int leadingJet=1) const;
   std::map<Int_t,Double_t> getHcalEMap(int leadingJet=1,
					double thrContrib=1e-4) const;

   // logical jet identifiers
   Double_t getSumEcalE_id(TJetId_t tag, int includeOthers=1) const {
     double val=-999.;
     switch(tag) {
     case _notJet: val=getSumEcalE(1,includeOthers); break;
     case _leadJet: val=getSumEcalE(0,includeOthers); break;
     case _subleadJet: val=getSumEcalE(2,includeOthers); break;
     default: std::cout << "getSumEcalE_id: not ready for tag=" << tag << "\n";
     }
     return val;
   }

   Double_t getSumHcalE_trackDiffEcal_id(TJetId_t jetId=_leadJet) const {
     const char *fncname="getSumHcalE_trackDiffEcal_id";
     double val=-999;
     switch(jetId) {
     case _notJet: std::cout << fncname << " bad jetId\n"; break;
     case _leadJet: val= getSumHcalE_trackDiffEcal(1); break;
     case _subleadJet: val= getSumHcalE_trackDiffEcal(0); break;
     default: std::cout << fncname << ": not ready for jetId=" << jetId << "\n";
     }
     return val;
   }

   std::map<Int_t,Double_t> getHcalEMap_id(TJetId_t jetId= _leadJet,
					double thrContrib=1e-4) const {
     const char *fncname="getHcalEMap_id";
     std::map<Int_t,Double_t> map;
     switch(jetId) {
     case _notJet: std::cout << fncname << " bad jetId\n"; break;
     case _leadJet: map=getHcalEMap(1,thrContrib); break;
     case _subleadJet: map=getHcalEMap(0,thrContrib); break;
     default: std::cout << fncname << ": not ready for jetId=" << jetId << "\n";
     }
     return map;
   }

};

#endif

// -------------------------------------------------------------
// Implementations
// ---------------------------------------------

#ifdef pf_gammajettree_cxx
pf_gammajettree::pf_gammajettree(TString fname) :
  fChain(new TChain("pf_gammajettree"))
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (fname.Length()==0) {
    fname="/media/ssd/ntuple-data-20141028/PhoJet_tree_Summer12DR53X_GPt170to300.root";
  }
  if (!this->Init(fname)) {
    std::cout << "Initialization failed in constructor" << std::endl;
  }
}

pf_gammajettree::~pf_gammajettree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pf_gammajettree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pf_gammajettree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

int pf_gammajettree::Init(const TString &fname)
{
  if (fname.Length()==0) {
    std::cout << "dijet_PFNtuple::Init: non-empty fname is expected\n";
    return 0;
  }
  fChain->Add(fname);

   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   photonTrig_fired = 0;
   photonTrig_prescale = 0;
   jetTrig_fired = 0;
   jetTrig_prescale = 0;
   tagPho_pfiso_mycharged03 = 0;
   ppfjet_had_E = 0;
   ppfjet_had_px = 0;
   ppfjet_had_py = 0;
   ppfjet_had_pz = 0;
   ppfjet_had_EcalE = 0;
   ppfjet_had_rawHcalE = 0;
   ppfjet_had_emf = 0;
   ppfjet_had_id = 0;
   ppfjet_had_candtrackind = 0;
   ppfjet_had_E_mctruth = 0;
   ppfjet_had_mcpdgId = 0;
   ppfjet_had_ntwrs = 0;
   ppfjet_twr_ieta = 0;
   ppfjet_twr_iphi = 0;
   ppfjet_twr_depth = 0;
   ppfjet_twr_subdet = 0;
   ppfjet_twr_hade = 0;
   ppfjet_twr_frac = 0;
   ppfjet_twr_candtrackind = 0;
   ppfjet_twr_hadind = 0;
   ppfjet_twr_elmttype = 0;
   ppfjet_twr_dR = 0;
   ppfjet_twr_clusterind = 0;
   ppfjet_cluster_eta = 0;
   ppfjet_cluster_phi = 0;
   ppfjet_cluster_dR = 0;
   ppfjet_candtrack_px = 0;
   ppfjet_candtrack_py = 0;
   ppfjet_candtrack_pz = 0;
   ppfjet_candtrack_EcalE = 0;
   pfjet2_had_E = 0;
   pfjet2_had_px = 0;
   pfjet2_had_py = 0;
   pfjet2_had_pz = 0;
   pfjet2_had_EcalE = 0;
   pfjet2_had_rawHcalE = 0;
   pfjet2_had_emf = 0;
   pfjet2_had_id = 0;
   pfjet2_had_candtrackind = 0;
   pfjet2_had_E_mctruth = 0;
   pfjet2_had_mcpdgId = 0;
   pfjet2_had_ntwrs = 0;
   pfjet2_twr_ieta = 0;
   pfjet2_twr_iphi = 0;
   pfjet2_twr_depth = 0;
   pfjet2_twr_subdet = 0;
   pfjet2_twr_hade = 0;
   pfjet2_twr_frac = 0;
   pfjet2_twr_candtrackind = 0;
   pfjet2_twr_hadind = 0;
   pfjet2_twr_elmttype = 0;
   pfjet2_twr_dR = 0;
   pfjet2_twr_clusterind = 0;
   pfjet2_cluster_eta = 0;
   pfjet2_cluster_phi = 0;
   pfjet2_cluster_dR = 0;
   pfjet2_candtrack_px = 0;
   pfjet2_candtrack_py = 0;
   pfjet2_candtrack_pz = 0;
   pfjet2_candtrack_EcalE = 0;

   // clear variables
   RunNumber=LumiBlock=EventNumber=0;
   EventWeight=rho2012=tagPho_pt=pho_2nd_pt=0;
   tagPho_energy=tagPho_eta=tagPho_phi=tagPho_sieie=0;
   tagPho_HoE=tagPho_r9=0;
   tagPho_EcalIsoDR04=tagPho_HcalIsoDR04=tagPho_HcalIsoDR0412=0;
   tagPho_TrkIsoHollowDR04=tagPho_pfiso_myphoton03=tagPho_pfiso_myneutral03=0;
   tagPho_pixelSeed=tagPho_ConvSafeEleVeto=0;
   tagPho_idTight=tagPho_idLoose=0;
   tagPho_genPt=tagPho_genEnergy=0;
   tagPho_genEta=tagPho_genPhi=tagPho_genDeltaR=0;
   nPhotons=nGenJets=nPFJets=0;
   ppfjet_pt=ppfjet_p=ppfjet_E=ppfjet_eta=ppfjet_phi=ppfjet_scale=0;
   ppfjet_E_NPVcorr=ppfjet_area=0;
   ppfjet_NeutralHadronFrac=ppfjet_NeutralEMFrac=0;
   ppfjet_nConstituents=0;
   ppfjet_ChargedHadronFrac=ppfjet_ChargedMultiplicity=ppfjet_ChargedEMFrac=0;
   ppfjet_genpt=ppfjet_genp=ppfjet_genE=ppfjet_gendr=0;
   ppfjet_unkown_E=ppfjet_electron_E=ppfjet_muon_E=ppfjet_photon_E=0;
   ppfjet_unkown_px=ppfjet_electron_px=ppfjet_muon_px=ppfjet_photon_px=0;
   ppfjet_unkown_py=ppfjet_electron_py=ppfjet_muon_py=ppfjet_photon_py=0;
   ppfjet_unkown_pz=ppfjet_electron_pz=ppfjet_muon_pz=ppfjet_photon_pz=0;
   ppfjet_unkown_EcalE=ppfjet_electron_EcalE=0;
   ppfjet_muon_EcalE=ppfjet_photon_EcalE=0;
   ppfjet_unkown_n=ppfjet_electron_n=ppfjet_muon_n=ppfjet_photon_n=0;
   ppfjet_had_n=ppfjet_ntwrs=ppfjet_cluster_n=ppfjet_ncandtracks=0;
   pfjet2_pt=pfjet2_p=pfjet2_E=pfjet2_eta=pfjet2_phi=pfjet2_scale=0;
   pfjet2_E_NPVcorr=pfjet2_area=0;
   pfjet2_NeutralHadronFrac=pfjet2_NeutralEMFrac=0;
   pfjet2_nConstituents=0;
   pfjet2_ChargedHadronFrac=pfjet2_ChargedMultiplicity=pfjet2_ChargedEMFrac=0;
   pfjet2_genpt=pfjet2_genp=pfjet2_genE=pfjet2_gendr=0;
   pfjet2_unkown_E=pfjet2_electron_E=pfjet2_muon_E=pfjet2_photon_E=0;
   pfjet2_unkown_px=pfjet2_electron_px=pfjet2_muon_px=pfjet2_photon_px=0;
   pfjet2_unkown_py=pfjet2_electron_py=pfjet2_muon_py=pfjet2_photon_py=0;
   pfjet2_unkown_pz=pfjet2_electron_pz=pfjet2_muon_pz=pfjet2_photon_pz=0;
   pfjet2_unkown_EcalE=pfjet2_electron_EcalE=0;
   pfjet2_muon_EcalE=pfjet2_photon_EcalE=0;
   pfjet2_unkown_n=pfjet2_electron_n=pfjet2_muon_n=pfjet2_photon_n=0;
   pfjet2_had_n=pfjet2_ntwrs=pfjet2_cluster_n=pfjet2_ncandtracks=0;
   pf_thirdjet_et=pf_thirdjet_pt=pf_thirdjet_p=0;
   pf_thirdjet_px=pf_thirdjet_py=pf_thirdjet_E=0;
   pf_thirdjet_eta=pf_thirdjet_phi=pf_thirdjet_scale=0;
   met_value=met_phi=met_sumEt=0;
   metType1_value=metType1_phi=metType1_sumEt=0;
   pf_NPV=0;

   // Set branch addresses and branch pointers
   //if (!tree) return;
   //fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("photonTrig_fired", &photonTrig_fired, &b_photonTrig_fired);
   fChain->SetBranchAddress("photonTrig_prescale", &photonTrig_prescale, &b_photonTrig_prescale);
   fChain->SetBranchAddress("jetTrig_fired", &jetTrig_fired, &b_jetTrig_fired);
   fChain->SetBranchAddress("jetTrig_prescale", &jetTrig_prescale, &b_jetTrig_prescale);
   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("EventWeight", &EventWeight, &b_EventWeight);
   fChain->SetBranchAddress("rho2012", &rho2012, &b_rho2012);
   fChain->SetBranchAddress("tagPho_pt", &tagPho_pt, &b_tagPho_pt);
   fChain->SetBranchAddress("pho_2nd_pt", &pho_2nd_pt, &b_pho_2nd_pt);
   fChain->SetBranchAddress("tagPho_energy", &tagPho_energy, &b_tagPho_energy);
   fChain->SetBranchAddress("tagPho_eta", &tagPho_eta, &b_tagPho_eta);
   fChain->SetBranchAddress("tagPho_phi", &tagPho_phi, &b_tagPho_phi);
   fChain->SetBranchAddress("tagPho_sieie", &tagPho_sieie, &b_tagPho_sieie);
   fChain->SetBranchAddress("tagPho_HoE", &tagPho_HoE, &b_tagPho_HoE);
   fChain->SetBranchAddress("tagPho_r9", &tagPho_r9, &b_tagPho_r9);
   fChain->SetBranchAddress("tagPho_EcalIsoDR04", &tagPho_EcalIsoDR04, &b_tagPho_EcalIsoDR04);
   fChain->SetBranchAddress("tagPho_HcalIsoDR04", &tagPho_HcalIsoDR04, &b_tagPho_HcalIsoDR04);
   fChain->SetBranchAddress("tagPho_HcalIsoDR0412", &tagPho_HcalIsoDR0412, &b_tagPho_HcalIsoDR0412);
   fChain->SetBranchAddress("tagPho_TrkIsoHollowDR04", &tagPho_TrkIsoHollowDR04, &b_tagPho_TrkIsoHollowDR04);
   fChain->SetBranchAddress("tagPho_pfiso_myphoton03", &tagPho_pfiso_myphoton03, &b_tagPho_pfiso_myphoton03);
   fChain->SetBranchAddress("tagPho_pfiso_myneutral03", &tagPho_pfiso_myneutral03, &b_tagPho_pfiso_myneutral03);
   fChain->SetBranchAddress("tagPho_pfiso_mycharged03", &tagPho_pfiso_mycharged03, &b_tagPho_pfiso_mycharged03);
   fChain->SetBranchAddress("tagPho_pixelSeed", &tagPho_pixelSeed, &b_tagPho_pixelSeed);
   fChain->SetBranchAddress("tagPho_ConvSafeEleVeto", &tagPho_ConvSafeEleVeto, &b_tagPho_ConvSafeEleVeto);
   fChain->SetBranchAddress("tagPho_idTight", &tagPho_idTight, &b_tagPho_idTight);
   fChain->SetBranchAddress("tagPho_idLoose", &tagPho_idLoose, &b_tagPho_idLoose);
   fChain->SetBranchAddress("tagPho_genPt", &tagPho_genPt, &b_tagPho_genPt);
   fChain->SetBranchAddress("tagPho_genEnergy", &tagPho_genEnergy, &b_tagPho_genEnergy);
   fChain->SetBranchAddress("tagPho_genEta", &tagPho_genEta, &b_tagPho_genEta);
   fChain->SetBranchAddress("tagPho_genPhi", &tagPho_genPhi, &b_tagPho_genPhi);
   fChain->SetBranchAddress("tagPho_genDeltaR", &tagPho_genDeltaR, &b_tagPho_genDeltaR);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("nGenJets", &nGenJets, &b_nGenJets);
   fChain->SetBranchAddress("nPFJets", &nPFJets, &b_nPFJets);
   fChain->SetBranchAddress("ppfjet_pt", &ppfjet_pt, &b_ppfjet_pt);
   fChain->SetBranchAddress("ppfjet_p", &ppfjet_p, &b_ppfjet_p);
   fChain->SetBranchAddress("ppfjet_E", &ppfjet_E, &b_ppfjet_E);
   fChain->SetBranchAddress("ppfjet_E_NPVcorr", &ppfjet_E_NPVcorr, &b_ppfjet_E_NPVcorr);
   fChain->SetBranchAddress("ppfjet_area", &ppfjet_area, &b_ppfjet_area);
   fChain->SetBranchAddress("ppfjet_eta", &ppfjet_eta, &b_ppfjet_eta);
   fChain->SetBranchAddress("ppfjet_phi", &ppfjet_phi, &b_ppfjet_phi);
   fChain->SetBranchAddress("ppfjet_scale", &ppfjet_scale, &b_ppfjet_scale);
   fChain->SetBranchAddress("ppfjet_NeutralHadronFrac", &ppfjet_NeutralHadronFrac, &b_ppfjet_NeutralHadronFrac);
   fChain->SetBranchAddress("ppfjet_NeutralEMFrac", &ppfjet_NeutralEMFrac, &b_ppfjet_NeutralEMFrac);
   fChain->SetBranchAddress("ppfjet_nConstituents", &ppfjet_nConstituents, &b_ppfjet_nConstituents);
   fChain->SetBranchAddress("ppfjet_ChargedHadronFrac", &ppfjet_ChargedHadronFrac, &b_ppfjet_ChargedHadronFrac);
   fChain->SetBranchAddress("ppfjet_ChargedMultiplicity", &ppfjet_ChargedMultiplicity, &b_ppfjet_ChargedMultiplicity);
   fChain->SetBranchAddress("ppfjet_ChargedEMFrac", &ppfjet_ChargedEMFrac, &b_ppfjet_ChargedEMFrac);
   fChain->SetBranchAddress("ppfjet_genpt", &ppfjet_genpt, &b_ppfjet_genpt);
   fChain->SetBranchAddress("ppfjet_genp", &ppfjet_genp, &b_ppfjet_genp);
   fChain->SetBranchAddress("ppfjet_genE", &ppfjet_genE, &b_ppfjet_genE);
   fChain->SetBranchAddress("ppfjet_gendr", &ppfjet_gendr, &b_ppfjet_gendr);
   fChain->SetBranchAddress("ppfjet_unkown_E", &ppfjet_unkown_E, &b_ppfjet_unkown_E);
   fChain->SetBranchAddress("ppfjet_electron_E", &ppfjet_electron_E, &b_ppfjet_electron_E);
   fChain->SetBranchAddress("ppfjet_muon_E", &ppfjet_muon_E, &b_ppfjet_muon_E);
   fChain->SetBranchAddress("ppfjet_photon_E", &ppfjet_photon_E, &b_ppfjet_photon_E);
   fChain->SetBranchAddress("ppfjet_unkown_px", &ppfjet_unkown_px, &b_ppfjet_unkown_px);
   fChain->SetBranchAddress("ppfjet_electron_px", &ppfjet_electron_px, &b_ppfjet_electron_px);
   fChain->SetBranchAddress("ppfjet_muon_px", &ppfjet_muon_px, &b_ppfjet_muon_px);
   fChain->SetBranchAddress("ppfjet_photon_px", &ppfjet_photon_px, &b_ppfjet_photon_px);
   fChain->SetBranchAddress("ppfjet_unkown_py", &ppfjet_unkown_py, &b_ppfjet_unkown_py);
   fChain->SetBranchAddress("ppfjet_electron_py", &ppfjet_electron_py, &b_ppfjet_electron_py);
   fChain->SetBranchAddress("ppfjet_muon_py", &ppfjet_muon_py, &b_ppfjet_muon_py);
   fChain->SetBranchAddress("ppfjet_photon_py", &ppfjet_photon_py, &b_ppfjet_photon_py);
   fChain->SetBranchAddress("ppfjet_unkown_pz", &ppfjet_unkown_pz, &b_ppfjet_unkown_pz);
   fChain->SetBranchAddress("ppfjet_electron_pz", &ppfjet_electron_pz, &b_ppfjet_electron_pz);
   fChain->SetBranchAddress("ppfjet_muon_pz", &ppfjet_muon_pz, &b_ppfjet_muon_pz);
   fChain->SetBranchAddress("ppfjet_photon_pz", &ppfjet_photon_pz, &b_ppfjet_photon_pz);
   fChain->SetBranchAddress("ppfjet_unkown_EcalE", &ppfjet_unkown_EcalE, &b_ppfjet_unkown_EcalE);
   fChain->SetBranchAddress("ppfjet_electron_EcalE", &ppfjet_electron_EcalE, &b_ppfjet_electron_EcalE);
   fChain->SetBranchAddress("ppfjet_muon_EcalE", &ppfjet_muon_EcalE, &b_ppfjet_muon_EcalE);
   fChain->SetBranchAddress("ppfjet_photon_EcalE", &ppfjet_photon_EcalE, &b_ppfjet_photon_EcalE);
   fChain->SetBranchAddress("ppfjet_unkown_n", &ppfjet_unkown_n, &b_ppfjet_unkown_n);
   fChain->SetBranchAddress("ppfjet_electron_n", &ppfjet_electron_n, &b_ppfjet_electron_n);
   fChain->SetBranchAddress("ppfjet_muon_n", &ppfjet_muon_n, &b_ppfjet_muon_n);
   fChain->SetBranchAddress("ppfjet_photon_n", &ppfjet_photon_n, &b_ppfjet_photon_n);
   fChain->SetBranchAddress("ppfjet_had_n", &ppfjet_had_n, &b_ppfjet_had_n);
   fChain->SetBranchAddress("ppfjet_had_E", &ppfjet_had_E, &b_ppfjet_had_E);
   fChain->SetBranchAddress("ppfjet_had_px", &ppfjet_had_px, &b_ppfjet_had_px);
   fChain->SetBranchAddress("ppfjet_had_py", &ppfjet_had_py, &b_ppfjet_had_py);
   fChain->SetBranchAddress("ppfjet_had_pz", &ppfjet_had_pz, &b_ppfjet_had_pz);
   fChain->SetBranchAddress("ppfjet_had_EcalE", &ppfjet_had_EcalE, &b_ppfjet_had_EcalE);
   fChain->SetBranchAddress("ppfjet_had_rawHcalE", &ppfjet_had_rawHcalE, &b_ppfjet_had_rawHcalE);
   fChain->SetBranchAddress("ppfjet_had_emf", &ppfjet_had_emf, &b_ppfjet_had_emf);
   fChain->SetBranchAddress("ppfjet_had_id", &ppfjet_had_id, &b_ppfjet_had_id);
   fChain->SetBranchAddress("ppfjet_had_candtrackind", &ppfjet_had_candtrackind, &b_ppfjet_had_candtrackind);
   fChain->SetBranchAddress("ppfjet_had_E_mctruth", &ppfjet_had_E_mctruth, &b_ppfjet_had_E_mctruth);
   fChain->SetBranchAddress("ppfjet_had_mcpdgId", &ppfjet_had_mcpdgId, &b_ppfjet_had_mcpdgId);
   fChain->SetBranchAddress("ppfjet_had_ntwrs", &ppfjet_had_ntwrs, &b_ppfjet_had_ntwrs);
   fChain->SetBranchAddress("ppfjet_ntwrs", &ppfjet_ntwrs, &b_ppfjet_ntwrs);
   fChain->SetBranchAddress("ppfjet_twr_ieta", &ppfjet_twr_ieta, &b_ppfjet_twr_ieta);
   fChain->SetBranchAddress("ppfjet_twr_iphi", &ppfjet_twr_iphi, &b_ppfjet_twr_iphi);
   fChain->SetBranchAddress("ppfjet_twr_depth", &ppfjet_twr_depth, &b_ppfjet_twr_depth);
   fChain->SetBranchAddress("ppfjet_twr_subdet", &ppfjet_twr_subdet, &b_ppfjet_twr_subdet);
   fChain->SetBranchAddress("ppfjet_twr_hade", &ppfjet_twr_hade, &b_ppfjet_twr_hade);
   fChain->SetBranchAddress("ppfjet_twr_frac", &ppfjet_twr_frac, &b_ppfjet_twr_frac);
   fChain->SetBranchAddress("ppfjet_twr_candtrackind", &ppfjet_twr_candtrackind, &b_ppfjet_twr_candtrackind);
   fChain->SetBranchAddress("ppfjet_twr_hadind", &ppfjet_twr_hadind, &b_ppfjet_twr_hadind);
   fChain->SetBranchAddress("ppfjet_twr_elmttype", &ppfjet_twr_elmttype, &b_ppfjet_twr_elmttype);
   fChain->SetBranchAddress("ppfjet_twr_dR", &ppfjet_twr_dR, &b_ppfjet_twr_dR);
   fChain->SetBranchAddress("ppfjet_twr_clusterind", &ppfjet_twr_clusterind, &b_ppfjet_twr_clusterind);
   fChain->SetBranchAddress("ppfjet_cluster_n", &ppfjet_cluster_n, &b_ppfjet_cluster_n);
   fChain->SetBranchAddress("ppfjet_cluster_eta", &ppfjet_cluster_eta, &b_ppfjet_cluster_eta);
   fChain->SetBranchAddress("ppfjet_cluster_phi", &ppfjet_cluster_phi, &b_ppfjet_cluster_phi);
   fChain->SetBranchAddress("ppfjet_cluster_dR", &ppfjet_cluster_dR, &b_ppfjet_cluster_dR);
   fChain->SetBranchAddress("ppfjet_ncandtracks", &ppfjet_ncandtracks, &b_ppfjet_ncandtracks);
   fChain->SetBranchAddress("ppfjet_candtrack_px", &ppfjet_candtrack_px, &b_ppfjet_candtrack_px);
   fChain->SetBranchAddress("ppfjet_candtrack_py", &ppfjet_candtrack_py, &b_ppfjet_candtrack_py);
   fChain->SetBranchAddress("ppfjet_candtrack_pz", &ppfjet_candtrack_pz, &b_ppfjet_candtrack_pz);
   fChain->SetBranchAddress("ppfjet_candtrack_EcalE", &ppfjet_candtrack_EcalE, &b_ppfjet_candtrack_EcalE);
   fChain->SetBranchAddress("pfjet2_pt", &pfjet2_pt, &b_pfjet2_pt);
   fChain->SetBranchAddress("pfjet2_p", &pfjet2_p, &b_pfjet2_p);
   fChain->SetBranchAddress("pfjet2_E", &pfjet2_E, &b_pfjet2_E);
   fChain->SetBranchAddress("pfjet2_E_NPVcorr", &pfjet2_E_NPVcorr, &b_pfjet2_E_NPVcorr);
   fChain->SetBranchAddress("pfjet2_area", &pfjet2_area, &b_pfjet2_area);
   fChain->SetBranchAddress("pfjet2_eta", &pfjet2_eta, &b_pfjet2_eta);
   fChain->SetBranchAddress("pfjet2_phi", &pfjet2_phi, &b_pfjet2_phi);
   fChain->SetBranchAddress("pfjet2_scale", &pfjet2_scale, &b_pfjet2_scale);
   fChain->SetBranchAddress("pfjet2_NeutralHadronFrac", &pfjet2_NeutralHadronFrac, &b_pfjet2_NeutralHadronFrac);
   fChain->SetBranchAddress("pfjet2_NeutralEMFrac", &pfjet2_NeutralEMFrac, &b_pfjet2_NeutralEMFrac);
   fChain->SetBranchAddress("pfjet2_nConstituents", &pfjet2_nConstituents, &b_pfjet2_nConstituents);
   fChain->SetBranchAddress("pfjet2_ChargedHadronFrac", &pfjet2_ChargedHadronFrac, &b_pfjet2_ChargedHadronFrac);
   fChain->SetBranchAddress("pfjet2_ChargedMultiplicity", &pfjet2_ChargedMultiplicity, &b_pfjet2_ChargedMultiplicity);
   fChain->SetBranchAddress("pfjet2_ChargedEMFrac", &pfjet2_ChargedEMFrac, &b_pfjet2_ChargedEMFrac);
   fChain->SetBranchAddress("pfjet2_genpt", &pfjet2_genpt, &b_pfjet2_genpt);
   fChain->SetBranchAddress("pfjet2_genp", &pfjet2_genp, &b_pfjet2_genp);
   fChain->SetBranchAddress("pfjet2_genE", &pfjet2_genE, &b_pfjet2_genE);
   fChain->SetBranchAddress("pfjet2_gendr", &pfjet2_gendr, &b_pfjet2_gendr);
   fChain->SetBranchAddress("pfjet2_unkown_E", &pfjet2_unkown_E, &b_pfjet2_unkown_E);
   fChain->SetBranchAddress("pfjet2_electron_E", &pfjet2_electron_E, &b_pfjet2_electron_E);
   fChain->SetBranchAddress("pfjet2_muon_E", &pfjet2_muon_E, &b_pfjet2_muon_E);
   fChain->SetBranchAddress("pfjet2_photon_E", &pfjet2_photon_E, &b_pfjet2_photon_E);
   fChain->SetBranchAddress("pfjet2_unkown_px", &pfjet2_unkown_px, &b_pfjet2_unkown_px);
   fChain->SetBranchAddress("pfjet2_electron_px", &pfjet2_electron_px, &b_pfjet2_electron_px);
   fChain->SetBranchAddress("pfjet2_muon_px", &pfjet2_muon_px, &b_pfjet2_muon_px);
   fChain->SetBranchAddress("pfjet2_photon_px", &pfjet2_photon_px, &b_pfjet2_photon_px);
   fChain->SetBranchAddress("pfjet2_unkown_py", &pfjet2_unkown_py, &b_pfjet2_unkown_py);
   fChain->SetBranchAddress("pfjet2_electron_py", &pfjet2_electron_py, &b_pfjet2_electron_py);
   fChain->SetBranchAddress("pfjet2_muon_py", &pfjet2_muon_py, &b_pfjet2_muon_py);
   fChain->SetBranchAddress("pfjet2_photon_py", &pfjet2_photon_py, &b_pfjet2_photon_py);
   fChain->SetBranchAddress("pfjet2_unkown_pz", &pfjet2_unkown_pz, &b_pfjet2_unkown_pz);
   fChain->SetBranchAddress("pfjet2_electron_pz", &pfjet2_electron_pz, &b_pfjet2_electron_pz);
   fChain->SetBranchAddress("pfjet2_muon_pz", &pfjet2_muon_pz, &b_pfjet2_muon_pz);
   fChain->SetBranchAddress("pfjet2_photon_pz", &pfjet2_photon_pz, &b_pfjet2_photon_pz);
   fChain->SetBranchAddress("pfjet2_unkown_EcalE", &pfjet2_unkown_EcalE, &b_pfjet2_unkown_EcalE);
   fChain->SetBranchAddress("pfjet2_electron_EcalE", &pfjet2_electron_EcalE, &b_pfjet2_electron_EcalE);
   fChain->SetBranchAddress("pfjet2_muon_EcalE", &pfjet2_muon_EcalE, &b_pfjet2_muon_EcalE);
   fChain->SetBranchAddress("pfjet2_photon_EcalE", &pfjet2_photon_EcalE, &b_pfjet2_photon_EcalE);
   fChain->SetBranchAddress("pfjet2_unkown_n", &pfjet2_unkown_n, &b_pfjet2_unkown_n);
   fChain->SetBranchAddress("pfjet2_electron_n", &pfjet2_electron_n, &b_pfjet2_electron_n);
   fChain->SetBranchAddress("pfjet2_muon_n", &pfjet2_muon_n, &b_pfjet2_muon_n);
   fChain->SetBranchAddress("pfjet2_photon_n", &pfjet2_photon_n, &b_pfjet2_photon_n);
   fChain->SetBranchAddress("pfjet2_had_n", &pfjet2_had_n, &b_pfjet2_had_n);
   fChain->SetBranchAddress("pfjet2_had_E", &pfjet2_had_E, &b_pfjet2_had_E);
   fChain->SetBranchAddress("pfjet2_had_px", &pfjet2_had_px, &b_pfjet2_had_px);
   fChain->SetBranchAddress("pfjet2_had_py", &pfjet2_had_py, &b_pfjet2_had_py);
   fChain->SetBranchAddress("pfjet2_had_pz", &pfjet2_had_pz, &b_pfjet2_had_pz);
   fChain->SetBranchAddress("pfjet2_had_EcalE", &pfjet2_had_EcalE, &b_pfjet2_had_EcalE);
   fChain->SetBranchAddress("pfjet2_had_rawHcalE", &pfjet2_had_rawHcalE, &b_pfjet2_had_rawHcalE);
   fChain->SetBranchAddress("pfjet2_had_emf", &pfjet2_had_emf, &b_pfjet2_had_emf);
   fChain->SetBranchAddress("pfjet2_had_id", &pfjet2_had_id, &b_pfjet2_had_id);
   fChain->SetBranchAddress("pfjet2_had_candtrackind", &pfjet2_had_candtrackind, &b_pfjet2_had_candtrackind);
   fChain->SetBranchAddress("pfjet2_had_E_mctruth", &pfjet2_had_E_mctruth, &b_pfjet2_had_E_mctruth);
   fChain->SetBranchAddress("pfjet2_had_mcpdgId", &pfjet2_had_mcpdgId, &b_pfjet2_had_mcpdgId);
   fChain->SetBranchAddress("pfjet2_had_ntwrs", &pfjet2_had_ntwrs, &b_pfjet2_had_ntwrs);
   fChain->SetBranchAddress("pfjet2_ntwrs", &pfjet2_ntwrs, &b_pfjet2_ntwrs);
   fChain->SetBranchAddress("pfjet2_twr_ieta", &pfjet2_twr_ieta, &b_pfjet2_twr_ieta);
   fChain->SetBranchAddress("pfjet2_twr_iphi", &pfjet2_twr_iphi, &b_pfjet2_twr_iphi);
   fChain->SetBranchAddress("pfjet2_twr_depth", &pfjet2_twr_depth, &b_pfjet2_twr_depth);
   fChain->SetBranchAddress("pfjet2_twr_subdet", &pfjet2_twr_subdet, &b_pfjet2_twr_subdet);
   fChain->SetBranchAddress("pfjet2_twr_hade", &pfjet2_twr_hade, &b_pfjet2_twr_hade);
   fChain->SetBranchAddress("pfjet2_twr_frac", &pfjet2_twr_frac, &b_pfjet2_twr_frac);
   fChain->SetBranchAddress("pfjet2_twr_candtrackind", &pfjet2_twr_candtrackind, &b_pfjet2_twr_candtrackind);
   fChain->SetBranchAddress("pfjet2_twr_hadind", &pfjet2_twr_hadind, &b_pfjet2_twr_hadind);
   fChain->SetBranchAddress("pfjet2_twr_elmttype", &pfjet2_twr_elmttype, &b_pfjet2_twr_elmttype);
   fChain->SetBranchAddress("pfjet2_twr_dR", &pfjet2_twr_dR, &b_pfjet2_twr_dR);
   fChain->SetBranchAddress("pfjet2_twr_clusterind", &pfjet2_twr_clusterind, &b_pfjet2_twr_clusterind);
   fChain->SetBranchAddress("pfjet2_cluster_n", &pfjet2_cluster_n, &b_pfjet2_cluster_n);
   fChain->SetBranchAddress("pfjet2_cluster_eta", &pfjet2_cluster_eta, &b_pfjet2_cluster_eta);
   fChain->SetBranchAddress("pfjet2_cluster_phi", &pfjet2_cluster_phi, &b_pfjet2_cluster_phi);
   fChain->SetBranchAddress("pfjet2_cluster_dR", &pfjet2_cluster_dR, &b_pfjet2_cluster_dR);
   fChain->SetBranchAddress("pfjet2_ncandtracks", &pfjet2_ncandtracks, &b_pfjet2_ncandtracks);
   fChain->SetBranchAddress("pfjet2_candtrack_px", &pfjet2_candtrack_px, &b_pfjet2_candtrack_px);
   fChain->SetBranchAddress("pfjet2_candtrack_py", &pfjet2_candtrack_py, &b_pfjet2_candtrack_py);
   fChain->SetBranchAddress("pfjet2_candtrack_pz", &pfjet2_candtrack_pz, &b_pfjet2_candtrack_pz);
   fChain->SetBranchAddress("pfjet2_candtrack_EcalE", &pfjet2_candtrack_EcalE, &b_pfjet2_candtrack_EcalE);
   fChain->SetBranchAddress("pf_thirdjet_et", &pf_thirdjet_et, &b_pf_thirdjet_et);
   fChain->SetBranchAddress("pf_thirdjet_pt", &pf_thirdjet_pt, &b_pf_thirdjet_pt);
   fChain->SetBranchAddress("pf_thirdjet_p", &pf_thirdjet_p, &b_pf_thirdjet_p);
   fChain->SetBranchAddress("pf_thirdjet_px", &pf_thirdjet_px, &b_pf_thirdjet_px);
   fChain->SetBranchAddress("pf_thirdjet_py", &pf_thirdjet_py, &b_pf_thirdjet_py);
   fChain->SetBranchAddress("pf_thirdjet_E", &pf_thirdjet_E, &b_pf_thirdjet_E);
   fChain->SetBranchAddress("pf_thirdjet_eta", &pf_thirdjet_eta, &b_pf_thirdjet_eta);
   fChain->SetBranchAddress("pf_thirdjet_phi", &pf_thirdjet_phi, &b_pf_thirdjet_phi);
   fChain->SetBranchAddress("pf_thirdjet_scale", &pf_thirdjet_scale, &b_pf_thirdjet_scale);
   fChain->SetBranchAddress("met_value", &met_value, &b_met_value);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_sumEt", &met_sumEt, &b_met_sumEt);
   fChain->SetBranchAddress("metType1_value", &metType1_value, &b_metType1_value);
   fChain->SetBranchAddress("metType1_phi", &metType1_phi, &b_metType1_phi);
   fChain->SetBranchAddress("metType1_sumEt", &metType1_sumEt, &b_metType1_sumEt);
   fChain->SetBranchAddress("pf_NPV", &pf_NPV, &b_pf_NPV);
   Notify();
   return 1;
}

Bool_t pf_gammajettree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void pf_gammajettree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pf_gammajettree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  if (0) std::cout << "entry=" << entry << "\n"; // satisfy the compiler
   return 1;
}

// a few useful methods
void pf_gammajettree::DeactivateBranches() {
  fChain->SetBranchStatus("*",0);
}

void pf_gammajettree::ActivateBranches(int count, ...) {
  va_list vl;
  va_start(vl,count);
  std::cout << "ActivateBranches(" << count << "): ";
  for (int i=0; i<count; ++i) {
    typedef const char* constCharPtr;
    TString brName= TString(va_arg(vl,constCharPtr));
    fChain->SetBranchStatus(brName,1);
    std::cout << " <" << brName << ">";
  }
  std::cout << "\n";
  va_end(vl);
}

void pf_gammajettree::ActivateBranches(const std::vector<TString> &brV) {
  unsigned int count=brV.size();
  std::cout << "ActivateBranches(" << count << "): ";
  for (unsigned int i=0; i<count; ++i) {
    fChain->SetBranchStatus(brV[i],1);
    std::cout << " <" << brV[i] << ">";
  }
  std::cout << "\n";
}
#endif // #ifdef pf_gammajettree_cxx
