//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct  8 12:17:06 2014 by ROOT version 5.34/21
// from TTree pf_dijettree/tree for dijet balancing using PFJets
// found on file: PhoJet_tree.root
//////////////////////////////////////////////////////////

#ifndef PhoJet_h
#define PhoJet_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class PhoJet {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         tagPho_et;
   Float_t         tagPho_energy;
   Float_t         pho_2nd_pt;
   Float_t         tagPho_eta;
   Float_t         tagPho_phi;
   Float_t         tagPho_sieie;
   Float_t         tagPho_HoE;
   Float_t         tagPho_r9;
   Int_t           tagPho_pixelSeed;
   Float_t         tagPho_TrkIsoHollowDR04;
   Float_t         tagPho_EcalIsoDR04;
   Float_t         tagPho_HcalIsoDR04;
   Float_t         tagPho_HcalIsoDR0412;
   Float_t         ppfjet_pt;
   Float_t         ppfjet_p;
   Float_t         ppfjet_E;
   Float_t         ppfjet_eta;
   Float_t         ppfjet_phi;
   Float_t         ppfjet_scale;
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
   Float_t         pf_2ndjet_pt;
   Int_t           pf_Run;
   Int_t           pf_Lumi;
   Int_t           pf_Event;
   Float_t         pf_weight;

   // List of branches
   TBranch        *b_tagPho_et;   //!
   TBranch        *b_tagPho_energy;   //!
   TBranch        *b_pho_2nd_pt;   //!
   TBranch        *b_tagPho_eta;   //!
   TBranch        *b_tagPho_phi;   //!
   TBranch        *b_tagPho_sieie;   //!
   TBranch        *b_tagPho_HoE;   //!
   TBranch        *b_tagPho_r9;   //!
   TBranch        *b_tagPho_pixelSeed;   //!
   TBranch        *b_tagPho_TrkIsoHollowDR04;   //!
   TBranch        *b_tagPho_EcalIsoDR04;   //!
   TBranch        *b_tagPho_HcalIsoDR04;   //!
   TBranch        *b_tagPho_HcalIsoDR0412;   //!
   TBranch        *b_ppfjet_pt;   //!
   TBranch        *b_ppfjet_p;   //!
   TBranch        *b_ppfjet_E;   //!
   TBranch        *b_ppfjet_eta;   //!
   TBranch        *b_ppfjet_phi;   //!
   TBranch        *b_ppfjet_scale;   //!
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
   TBranch        *b_pf_2ndjet_pt;   //!
   TBranch        *b_pf_Run;   //!
   TBranch        *b_pf_Lumi;   //!
   TBranch        *b_pf_Event;   //!
   TBranch        *b_pf_weight;   //!

   PhoJet(TTree *tree=0);
   virtual ~PhoJet();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PhoJet_cxx
PhoJet::PhoJet(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("PhoJet_tree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("PhoJet_tree.root");
      }
      f->GetObject("pf_dijettree",tree);

   }
   Init(tree);
}

PhoJet::~PhoJet()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PhoJet::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PhoJet::LoadTree(Long64_t entry)
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

void PhoJet::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ppfjet_had_E = 0;
   ppfjet_had_px = 0;
   ppfjet_had_py = 0;
   ppfjet_had_pz = 0;
   ppfjet_had_EcalE = 0;
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
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("tagPho_et", &tagPho_et, &b_tagPho_et);
   fChain->SetBranchAddress("tagPho_energy", &tagPho_energy, &b_tagPho_energy);
   fChain->SetBranchAddress("pho_2nd_pt", &pho_2nd_pt, &b_pho_2nd_pt);
   fChain->SetBranchAddress("tagPho_eta", &tagPho_eta, &b_tagPho_eta);
   fChain->SetBranchAddress("tagPho_phi", &tagPho_phi, &b_tagPho_phi);
   fChain->SetBranchAddress("tagPho_sieie", &tagPho_sieie, &b_tagPho_sieie);
   fChain->SetBranchAddress("tagPho_HoE", &tagPho_HoE, &b_tagPho_HoE);
   fChain->SetBranchAddress("tagPho_r9", &tagPho_r9, &b_tagPho_r9);
   fChain->SetBranchAddress("tagPho_pixelSeed", &tagPho_pixelSeed, &b_tagPho_pixelSeed);
   fChain->SetBranchAddress("tagPho_TrkIsoHollowDR04", &tagPho_TrkIsoHollowDR04, &b_tagPho_TrkIsoHollowDR04);
   fChain->SetBranchAddress("tagPho_EcalIsoDR04", &tagPho_EcalIsoDR04, &b_tagPho_EcalIsoDR04);
   fChain->SetBranchAddress("tagPho_HcalIsoDR04", &tagPho_HcalIsoDR04, &b_tagPho_HcalIsoDR04);
   fChain->SetBranchAddress("tagPho_HcalIsoDR0412", &tagPho_HcalIsoDR0412, &b_tagPho_HcalIsoDR0412);
   fChain->SetBranchAddress("ppfjet_pt", &ppfjet_pt, &b_ppfjet_pt);
   fChain->SetBranchAddress("ppfjet_p", &ppfjet_p, &b_ppfjet_p);
   fChain->SetBranchAddress("ppfjet_E", &ppfjet_E, &b_ppfjet_E);
   fChain->SetBranchAddress("ppfjet_eta", &ppfjet_eta, &b_ppfjet_eta);
   fChain->SetBranchAddress("ppfjet_phi", &ppfjet_phi, &b_ppfjet_phi);
   fChain->SetBranchAddress("ppfjet_scale", &ppfjet_scale, &b_ppfjet_scale);
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
   fChain->SetBranchAddress("pf_2ndjet_pt", &pf_2ndjet_pt, &b_pf_2ndjet_pt);
   fChain->SetBranchAddress("pf_Run", &pf_Run, &b_pf_Run);
   fChain->SetBranchAddress("pf_Lumi", &pf_Lumi, &b_pf_Lumi);
   fChain->SetBranchAddress("pf_Event", &pf_Event, &b_pf_Event);
   fChain->SetBranchAddress("pf_weight", &pf_weight, &b_pf_weight);
   Notify();
}

Bool_t PhoJet::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PhoJet::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PhoJet::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PhoJet_cxx
