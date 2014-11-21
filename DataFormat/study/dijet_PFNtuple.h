//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct  7 13:09:05 2014 by ROOT version 5.34/09
// from TTree pf_dijettree/tree for dijet balancing using PFJets
// found on file: /media/spektras/dijetData/tree_1_1_FuQ.root
//////////////////////////////////////////////////////////

#ifndef dijet_PFNtuple_h
#define dijet_PFNtuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Other classes
#include <iostream>
#include <stdarg.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

class dijet_PFNtuple {
public :
  //TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  TChain         *fChain; // modified the type
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         tpfjet_pt;
   Float_t         tpfjet_p;
   Float_t         tpfjet_E;
   Float_t         tpfjet_eta;
   Float_t         tpfjet_phi;
   Float_t         tpfjet_EMfrac;
   Float_t         tpfjet_hadEcalEfrac;
   Float_t         tpfjet_scale;
   Float_t         tpfjet_genpt;
   Float_t         tpfjet_genp;
   Float_t         tpfjet_genE;
   Float_t         tpfjet_gendr;
   Float_t         tpfjet_unkown_E;
   Float_t         tpfjet_electron_E;
   Float_t         tpfjet_muon_E;
   Float_t         tpfjet_photon_E;
   Float_t         tpfjet_unkown_px;
   Float_t         tpfjet_electron_px;
   Float_t         tpfjet_muon_px;
   Float_t         tpfjet_photon_px;
   Float_t         tpfjet_unkown_py;
   Float_t         tpfjet_electron_py;
   Float_t         tpfjet_muon_py;
   Float_t         tpfjet_photon_py;
   Float_t         tpfjet_unkown_pz;
   Float_t         tpfjet_electron_pz;
   Float_t         tpfjet_muon_pz;
   Float_t         tpfjet_photon_pz;
   Float_t         tpfjet_unkown_EcalE;
   Float_t         tpfjet_electron_EcalE;
   Float_t         tpfjet_muon_EcalE;
   Float_t         tpfjet_photon_EcalE;
   Int_t           tpfjet_unkown_n;
   Int_t           tpfjet_electron_n;
   Int_t           tpfjet_muon_n;
   Int_t           tpfjet_photon_n;
   Int_t           tpfjet_had_n;
   vector<float>   *tpfjet_had_E;
   vector<float>   *tpfjet_had_px;
   vector<float>   *tpfjet_had_py;
   vector<float>   *tpfjet_had_pz;
   vector<float>   *tpfjet_had_EcalE;
   vector<float>   *tpfjet_had_rawHcalE;
   vector<float>   *tpfjet_had_emf;
   vector<int>     *tpfjet_had_id;
   vector<int>     *tpfjet_had_candtrackind;
   vector<float>   *tpfjet_had_E_mctruth;
   vector<int>     *tpfjet_had_mcpdgId;
   vector<int>     *tpfjet_had_ntwrs;
   Int_t           tpfjet_ntwrs;
   vector<int>     *tpfjet_twr_ieta;
   vector<int>     *tpfjet_twr_iphi;
   vector<int>     *tpfjet_twr_depth;
   vector<int>     *tpfjet_twr_subdet;
   vector<float>   *tpfjet_twr_hade;
   vector<float>   *tpfjet_twr_frac;
   vector<int>     *tpfjet_twr_candtrackind;
   vector<int>     *tpfjet_twr_hadind;
   vector<int>     *tpfjet_twr_elmttype;
   vector<float>   *tpfjet_twr_dR;
   vector<int>     *tpfjet_twr_clusterind;
   Int_t           tpfjet_cluster_n;
   vector<float>   *tpfjet_cluster_eta;
   vector<float>   *tpfjet_cluster_phi;
   vector<float>   *tpfjet_cluster_dR;
   Int_t           tpfjet_ncandtracks;
   vector<float>   *tpfjet_candtrack_px;
   vector<float>   *tpfjet_candtrack_py;
   vector<float>   *tpfjet_candtrack_pz;
   vector<float>   *tpfjet_candtrack_EcalE;
   Float_t         ppfjet_pt;
   Float_t         ppfjet_p;
   Float_t         ppfjet_E;
   Float_t         ppfjet_eta;
   Float_t         ppfjet_phi;
   Float_t         ppfjet_EMfrac;
   Float_t         ppfjet_hadEcalEfrac;
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
   Float_t         pf_dijet_deta;
   Float_t         pf_dijet_dphi;
   Float_t         pf_dijet_balance;
   Float_t         pf_thirdjet_px;
   Float_t         pf_thirdjet_py;
   Int_t           pf_Run;
   Int_t           pf_Lumi;
   Int_t           pf_Event;
   Float_t         pf_weight;

   // List of branches
   TBranch        *b_tpfjet_pt;   //!
   TBranch        *b_tpfjet_p;   //!
   TBranch        *b_tpfjet_E;   //!
   TBranch        *b_tpfjet_eta;   //!
   TBranch        *b_tpfjet_phi;   //!
   TBranch        *b_tpfjet_EMfrac;   //!
   TBranch        *b_tpfjet_hadEcalEfrac;   //!
   TBranch        *b_tpfjet_scale;   //!
   TBranch        *b_tpfjet_genpt;   //!
   TBranch        *b_tpfjet_genp;   //!
   TBranch        *b_tpfjet_genE;   //!
   TBranch        *b_tpfjet_gendr;   //!
   TBranch        *b_tpfjet_unkown_E;   //!
   TBranch        *b_tpfjet_electron_E;   //!
   TBranch        *b_tpfjet_muon_E;   //!
   TBranch        *b_tpfjet_photon_E;   //!
   TBranch        *b_tpfjet_unkown_px;   //!
   TBranch        *b_tpfjet_electron_px;   //!
   TBranch        *b_tpfjet_muon_px;   //!
   TBranch        *b_tpfjet_photon_px;   //!
   TBranch        *b_tpfjet_unkown_py;   //!
   TBranch        *b_tpfjet_electron_py;   //!
   TBranch        *b_tpfjet_muon_py;   //!
   TBranch        *b_tpfjet_photon_py;   //!
   TBranch        *b_tpfjet_unkown_pz;   //!
   TBranch        *b_tpfjet_electron_pz;   //!
   TBranch        *b_tpfjet_muon_pz;   //!
   TBranch        *b_tpfjet_photon_pz;   //!
   TBranch        *b_tpfjet_unkown_EcalE;   //!
   TBranch        *b_tpfjet_electron_EcalE;   //!
   TBranch        *b_tpfjet_muon_EcalE;   //!
   TBranch        *b_tpfjet_photon_EcalE;   //!
   TBranch        *b_tpfjet_unkown_n;   //!
   TBranch        *b_tpfjet_electron_n;   //!
   TBranch        *b_tpfjet_muon_n;   //!
   TBranch        *b_tpfjet_photon_n;   //!
   TBranch        *b_tpfjet_had_n;   //!
   TBranch        *b_tpfjet_had_E;   //!
   TBranch        *b_tpfjet_had_px;   //!
   TBranch        *b_tpfjet_had_py;   //!
   TBranch        *b_tpfjet_had_pz;   //!
   TBranch        *b_tpfjet_had_EcalE;   //!
   TBranch        *b_tpfjet_had_rawHcalE;   //!
   TBranch        *b_tpfjet_had_emf;   //!
   TBranch        *b_tpfjet_had_id;   //!
   TBranch        *b_tpfjet_had_candtrackind;   //!
   TBranch        *b_tpfjet_had_E_mctruth;   //!
   TBranch        *b_tpfjet_had_mcpdgId;   //!
   TBranch        *b_tpfjet_had_ntwrs;   //!
   TBranch        *b_tpfjet_ntwrs;   //!
   TBranch        *b_tpfjet_twr_ieta;   //!
   TBranch        *b_tpfjet_twr_iphi;   //!
   TBranch        *b_tpfjet_twr_depth;   //!
   TBranch        *b_tpfjet_twr_subdet;   //!
   TBranch        *b_tpfjet_twr_hade;   //!
   TBranch        *b_tpfjet_twr_frac;   //!
   TBranch        *b_tpfjet_twr_candtrackind;   //!
   TBranch        *b_tpfjet_twr_hadind;   //!
   TBranch        *b_tpfjet_twr_elmttype;   //!
   TBranch        *b_tpfjet_twr_dR;   //!
   TBranch        *b_tpfjet_twr_clusterind;   //!
   TBranch        *b_tpfjet_cluster_n;   //!
   TBranch        *b_tpfjet_cluster_eta;   //!
   TBranch        *b_tpfjet_cluster_phi;   //!
   TBranch        *b_tpfjet_cluster_dR;   //!
   TBranch        *b_tpfjet_ncandtracks;   //!
   TBranch        *b_tpfjet_candtrack_px;   //!
   TBranch        *b_tpfjet_candtrack_py;   //!
   TBranch        *b_tpfjet_candtrack_pz;   //!
   TBranch        *b_tpfjet_candtrack_EcalE;   //!
   TBranch        *b_ppfjet_pt;   //!
   TBranch        *b_ppfjet_p;   //!
   TBranch        *b_ppfjet_E;   //!
   TBranch        *b_ppfjet_eta;   //!
   TBranch        *b_ppfjet_phi;   //!
   TBranch        *b_ppfjet_EMfrac;   //!
   TBranch        *b_ppfjet_hadEcalEfrac;   //!
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
   TBranch        *b_pf_dijet_deta;   //!
   TBranch        *b_pf_dijet_dphi;   //!
   TBranch        *b_pf_dijet_balance;   //!
   TBranch        *b_pf_thirdjet_px;   //!
   TBranch        *b_pf_thirdjet_py;   //!
   TBranch        *b_pf_Run;   //!
   TBranch        *b_pf_Lumi;   //!
   TBranch        *b_pf_Event;   //!
   TBranch        *b_pf_weight;   //!

   dijet_PFNtuple(TString fname="");
   virtual ~dijet_PFNtuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual int      Init(const TString &fname);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void DeactivateBranches();
   void ActivateBranches(int count, ...); // list n branch names
   void ActivateBranches(const std::vector<TString> &brV); // list of branch names

   friend
     std::ostream& operator<<(std::ostream &out, dijet_PFNtuple &obj) {
     if (out==std::cout) obj.Show();
     else out << "cannot print dijet_PFNtuple\n";
     return out;
   }

   void PrintSelectedFields(int selection=0);
   // numeric jet identifiers
   Double_t getSumEcalE(int tagJet, int includeOthers=1) const;
   Double_t getSumHcalE_trackDiffEcal(int tagJet) const;
   std::map<Int_t,Double_t> getHcalEMap(int tagJet,
					double thrContrib=1e-4) const;
};

#endif

#ifdef dijet_PFNtuple_cxx
dijet_PFNtuple::dijet_PFNtuple(TString fname) :
  fChain(new TChain("pf_dijettree")) 
{
  if (fname.Length()==0) {
    fname="/media/spektras/dijetData/tree_*.root";
  }
  if (!this->Init(fname)) {
    std::cout << "Initialization failed in constructor" << std::endl;
  }
}

dijet_PFNtuple::~dijet_PFNtuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t dijet_PFNtuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   tpfjet_twr_ieta->clear();
   tpfjet_twr_iphi->clear();
   ppfjet_twr_ieta->clear();
   ppfjet_twr_iphi->clear();
   return fChain->GetEntry(entry);
}
Long64_t dijet_PFNtuple::LoadTree(Long64_t entry)
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

int dijet_PFNtuple::Init(const TString &fname)
{
  if (fname.Length()==0) {
    std::cout << "dijet_PFNtuple::Init: non-empty fname is expected\n";
    return 0;
  }
  fChain->Add(fname);
  std::cout << "chain added the name " << fname << std::endl;

   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   tpfjet_had_E = 0;
   tpfjet_had_px = 0;
   tpfjet_had_py = 0;
   tpfjet_had_pz = 0;
   tpfjet_had_EcalE = 0;
   tpfjet_had_rawHcalE = 0;
   tpfjet_had_emf = 0;
   tpfjet_had_id = 0;
   tpfjet_had_candtrackind = 0;
   tpfjet_had_E_mctruth = 0;
   tpfjet_had_mcpdgId = 0;
   tpfjet_had_ntwrs = 0;
   tpfjet_twr_ieta = 0;
   tpfjet_twr_iphi = 0;
   tpfjet_twr_depth = 0;
   tpfjet_twr_subdet = 0;
   tpfjet_twr_hade = 0;
   tpfjet_twr_frac = 0;
   tpfjet_twr_candtrackind = 0;
   tpfjet_twr_hadind = 0;
   tpfjet_twr_elmttype = 0;
   tpfjet_twr_dR = 0;
   tpfjet_twr_clusterind = 0;
   tpfjet_cluster_eta = 0;
   tpfjet_cluster_phi = 0;
   tpfjet_cluster_dR = 0;
   tpfjet_candtrack_px = 0;
   tpfjet_candtrack_py = 0;
   tpfjet_candtrack_pz = 0;
   tpfjet_candtrack_EcalE = 0;
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
   // Set branch addresses and branch pointers
   //if (!tree) return;
   //fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("tpfjet_pt", &tpfjet_pt, &b_tpfjet_pt);
   fChain->SetBranchAddress("tpfjet_p", &tpfjet_p, &b_tpfjet_p);
   fChain->SetBranchAddress("tpfjet_E", &tpfjet_E, &b_tpfjet_E);
   fChain->SetBranchAddress("tpfjet_eta", &tpfjet_eta, &b_tpfjet_eta);
   fChain->SetBranchAddress("tpfjet_phi", &tpfjet_phi, &b_tpfjet_phi);
   fChain->SetBranchAddress("tpfjet_EMfrac", &tpfjet_EMfrac, &b_tpfjet_EMfrac);
   fChain->SetBranchAddress("tpfjet_hadEcalEfrac", &tpfjet_hadEcalEfrac, &b_tpfjet_hadEcalEfrac);
   fChain->SetBranchAddress("tpfjet_scale", &tpfjet_scale, &b_tpfjet_scale);
   fChain->SetBranchAddress("tpfjet_genpt", &tpfjet_genpt, &b_tpfjet_genpt);
   fChain->SetBranchAddress("tpfjet_genp", &tpfjet_genp, &b_tpfjet_genp);
   fChain->SetBranchAddress("tpfjet_genE", &tpfjet_genE, &b_tpfjet_genE);
   fChain->SetBranchAddress("tpfjet_gendr", &tpfjet_gendr, &b_tpfjet_gendr);
   fChain->SetBranchAddress("tpfjet_unkown_E", &tpfjet_unkown_E, &b_tpfjet_unkown_E);
   fChain->SetBranchAddress("tpfjet_electron_E", &tpfjet_electron_E, &b_tpfjet_electron_E);
   fChain->SetBranchAddress("tpfjet_muon_E", &tpfjet_muon_E, &b_tpfjet_muon_E);
   fChain->SetBranchAddress("tpfjet_photon_E", &tpfjet_photon_E, &b_tpfjet_photon_E);
   fChain->SetBranchAddress("tpfjet_unkown_px", &tpfjet_unkown_px, &b_tpfjet_unkown_px);
   fChain->SetBranchAddress("tpfjet_electron_px", &tpfjet_electron_px, &b_tpfjet_electron_px);
   fChain->SetBranchAddress("tpfjet_muon_px", &tpfjet_muon_px, &b_tpfjet_muon_px);
   fChain->SetBranchAddress("tpfjet_photon_px", &tpfjet_photon_px, &b_tpfjet_photon_px);
   fChain->SetBranchAddress("tpfjet_unkown_py", &tpfjet_unkown_py, &b_tpfjet_unkown_py);
   fChain->SetBranchAddress("tpfjet_electron_py", &tpfjet_electron_py, &b_tpfjet_electron_py);
   fChain->SetBranchAddress("tpfjet_muon_py", &tpfjet_muon_py, &b_tpfjet_muon_py);
   fChain->SetBranchAddress("tpfjet_photon_py", &tpfjet_photon_py, &b_tpfjet_photon_py);
   fChain->SetBranchAddress("tpfjet_unkown_pz", &tpfjet_unkown_pz, &b_tpfjet_unkown_pz);
   fChain->SetBranchAddress("tpfjet_electron_pz", &tpfjet_electron_pz, &b_tpfjet_electron_pz);
   fChain->SetBranchAddress("tpfjet_muon_pz", &tpfjet_muon_pz, &b_tpfjet_muon_pz);
   fChain->SetBranchAddress("tpfjet_photon_pz", &tpfjet_photon_pz, &b_tpfjet_photon_pz);
   fChain->SetBranchAddress("tpfjet_unkown_EcalE", &tpfjet_unkown_EcalE, &b_tpfjet_unkown_EcalE);
   fChain->SetBranchAddress("tpfjet_electron_EcalE", &tpfjet_electron_EcalE, &b_tpfjet_electron_EcalE);
   fChain->SetBranchAddress("tpfjet_muon_EcalE", &tpfjet_muon_EcalE, &b_tpfjet_muon_EcalE);
   fChain->SetBranchAddress("tpfjet_photon_EcalE", &tpfjet_photon_EcalE, &b_tpfjet_photon_EcalE);
   fChain->SetBranchAddress("tpfjet_unkown_n", &tpfjet_unkown_n, &b_tpfjet_unkown_n);
   fChain->SetBranchAddress("tpfjet_electron_n", &tpfjet_electron_n, &b_tpfjet_electron_n);
   fChain->SetBranchAddress("tpfjet_muon_n", &tpfjet_muon_n, &b_tpfjet_muon_n);
   fChain->SetBranchAddress("tpfjet_photon_n", &tpfjet_photon_n, &b_tpfjet_photon_n);
   fChain->SetBranchAddress("tpfjet_had_n", &tpfjet_had_n, &b_tpfjet_had_n);
   fChain->SetBranchAddress("tpfjet_had_E", &tpfjet_had_E, &b_tpfjet_had_E);
   fChain->SetBranchAddress("tpfjet_had_px", &tpfjet_had_px, &b_tpfjet_had_px);
   fChain->SetBranchAddress("tpfjet_had_py", &tpfjet_had_py, &b_tpfjet_had_py);
   fChain->SetBranchAddress("tpfjet_had_pz", &tpfjet_had_pz, &b_tpfjet_had_pz);
   fChain->SetBranchAddress("tpfjet_had_EcalE", &tpfjet_had_EcalE, &b_tpfjet_had_EcalE);
   fChain->SetBranchAddress("tpfjet_had_rawHcalE", &tpfjet_had_rawHcalE, &b_tpfjet_had_rawHcalE);
   fChain->SetBranchAddress("tpfjet_had_emf", &tpfjet_had_emf, &b_tpfjet_had_emf);
   fChain->SetBranchAddress("tpfjet_had_id", &tpfjet_had_id, &b_tpfjet_had_id);
   fChain->SetBranchAddress("tpfjet_had_candtrackind", &tpfjet_had_candtrackind, &b_tpfjet_had_candtrackind);
   fChain->SetBranchAddress("tpfjet_had_E_mctruth", &tpfjet_had_E_mctruth, &b_tpfjet_had_E_mctruth);
   fChain->SetBranchAddress("tpfjet_had_mcpdgId", &tpfjet_had_mcpdgId, &b_tpfjet_had_mcpdgId);
   fChain->SetBranchAddress("tpfjet_had_ntwrs", &tpfjet_had_ntwrs, &b_tpfjet_had_ntwrs);
   fChain->SetBranchAddress("tpfjet_ntwrs", &tpfjet_ntwrs, &b_tpfjet_ntwrs);
   fChain->SetBranchAddress("tpfjet_twr_ieta", &tpfjet_twr_ieta, &b_tpfjet_twr_ieta);
   fChain->SetBranchAddress("tpfjet_twr_iphi", &tpfjet_twr_iphi, &b_tpfjet_twr_iphi);
   fChain->SetBranchAddress("tpfjet_twr_depth", &tpfjet_twr_depth, &b_tpfjet_twr_depth);
   fChain->SetBranchAddress("tpfjet_twr_subdet", &tpfjet_twr_subdet, &b_tpfjet_twr_subdet);
   fChain->SetBranchAddress("tpfjet_twr_hade", &tpfjet_twr_hade, &b_tpfjet_twr_hade);
   fChain->SetBranchAddress("tpfjet_twr_frac", &tpfjet_twr_frac, &b_tpfjet_twr_frac);
   fChain->SetBranchAddress("tpfjet_twr_candtrackind", &tpfjet_twr_candtrackind, &b_tpfjet_twr_candtrackind);
   fChain->SetBranchAddress("tpfjet_twr_hadind", &tpfjet_twr_hadind, &b_tpfjet_twr_hadind);
   fChain->SetBranchAddress("tpfjet_twr_elmttype", &tpfjet_twr_elmttype, &b_tpfjet_twr_elmttype);
   fChain->SetBranchAddress("tpfjet_twr_dR", &tpfjet_twr_dR, &b_tpfjet_twr_dR);
   fChain->SetBranchAddress("tpfjet_twr_clusterind", &tpfjet_twr_clusterind, &b_tpfjet_twr_clusterind);
   fChain->SetBranchAddress("tpfjet_cluster_n", &tpfjet_cluster_n, &b_tpfjet_cluster_n);
   fChain->SetBranchAddress("tpfjet_cluster_eta", &tpfjet_cluster_eta, &b_tpfjet_cluster_eta);
   fChain->SetBranchAddress("tpfjet_cluster_phi", &tpfjet_cluster_phi, &b_tpfjet_cluster_phi);
   fChain->SetBranchAddress("tpfjet_cluster_dR", &tpfjet_cluster_dR, &b_tpfjet_cluster_dR);
   fChain->SetBranchAddress("tpfjet_ncandtracks", &tpfjet_ncandtracks, &b_tpfjet_ncandtracks);
   fChain->SetBranchAddress("tpfjet_candtrack_px", &tpfjet_candtrack_px, &b_tpfjet_candtrack_px);
   fChain->SetBranchAddress("tpfjet_candtrack_py", &tpfjet_candtrack_py, &b_tpfjet_candtrack_py);
   fChain->SetBranchAddress("tpfjet_candtrack_pz", &tpfjet_candtrack_pz, &b_tpfjet_candtrack_pz);
   fChain->SetBranchAddress("tpfjet_candtrack_EcalE", &tpfjet_candtrack_EcalE, &b_tpfjet_candtrack_EcalE);
   fChain->SetBranchAddress("ppfjet_pt", &ppfjet_pt, &b_ppfjet_pt);
   fChain->SetBranchAddress("ppfjet_p", &ppfjet_p, &b_ppfjet_p);
   fChain->SetBranchAddress("ppfjet_E", &ppfjet_E, &b_ppfjet_E);
   fChain->SetBranchAddress("ppfjet_eta", &ppfjet_eta, &b_ppfjet_eta);
   fChain->SetBranchAddress("ppfjet_phi", &ppfjet_phi, &b_ppfjet_phi);
   fChain->SetBranchAddress("ppfjet_EMfrac", &ppfjet_EMfrac, &b_ppfjet_EMfrac);
   fChain->SetBranchAddress("ppfjet_hadEcalEfrac", &ppfjet_hadEcalEfrac, &b_ppfjet_hadEcalEfrac);
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
   fChain->SetBranchAddress("pf_dijet_deta", &pf_dijet_deta, &b_pf_dijet_deta);
   fChain->SetBranchAddress("pf_dijet_dphi", &pf_dijet_dphi, &b_pf_dijet_dphi);
   fChain->SetBranchAddress("pf_dijet_balance", &pf_dijet_balance, &b_pf_dijet_balance);
   fChain->SetBranchAddress("pf_thirdjet_px", &pf_thirdjet_px, &b_pf_thirdjet_px);
   fChain->SetBranchAddress("pf_thirdjet_py", &pf_thirdjet_py, &b_pf_thirdjet_py);
   fChain->SetBranchAddress("pf_Run", &pf_Run, &b_pf_Run);
   fChain->SetBranchAddress("pf_Lumi", &pf_Lumi, &b_pf_Lumi);
   fChain->SetBranchAddress("pf_Event", &pf_Event, &b_pf_Event);
   fChain->SetBranchAddress("pf_weight", &pf_weight, &b_pf_weight);
   Notify();
   return 1;
}

Bool_t dijet_PFNtuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void dijet_PFNtuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t dijet_PFNtuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  if (entry==0) std::cout << "\n"; // harmless; to prevent compiler complaints
   return 1;
}

void dijet_PFNtuple::DeactivateBranches() {
  fChain->SetBranchStatus("*",0);
}

void dijet_PFNtuple::ActivateBranches(int count, ...) {
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

void dijet_PFNtuple::ActivateBranches(const std::vector<TString> &brV) {
  unsigned int count=brV.size();
  std::cout << "ActivateBranches(" << count << "): ";
  for (unsigned int i=0; i<count; ++i) {
    fChain->SetBranchStatus(brV[i],1);
    std::cout << " <" << brV[i] << ">";
  }
  std::cout << "\n";
}

#endif // #ifdef dijet_PFNtuple_cxx
