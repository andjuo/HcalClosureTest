//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 13 15:38:02 2015 by ROOT version 5.34/09
// from TTree misc_tree/tree for misc.info
// found on file: PhoJet_tree_CHS_from_AOD.root
//////////////////////////////////////////////////////////

#ifndef misc_gammajettree_h
#define misc_gammajettree_h

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

class misc_gammajettree {
public :
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Bool_t          ignoreHLT;
   Bool_t          doCaloJets;
   Bool_t          doPFJets;
   Bool_t          doGenJets;
   Bool_t          workOnAOD;
   vector<string>  *photonTriggerNames;
   vector<string>  *jetTriggerNames;
   ULong64_t       nProcessed;

   // List of branches
   TBranch        *b_ignoreHLT;   //!
   TBranch        *b_doCaloJets;   //!
   TBranch        *b_doPFJets;   //!
   TBranch        *b_doGenJets;   //!
   TBranch        *b_workOnAOD;   //!
   TBranch        *b_photonTriggerNames;   //!
   TBranch        *b_jetTriggerNames;   //!
   TBranch        *b_nProcessed;   //!

   misc_gammajettree(TString fname="");
   virtual ~misc_gammajettree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual int      Init(const TString &fname);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   // Added methods
   friend
     std::ostream& operator<<(std::ostream &out, misc_gammajettree &obj) {
     if (out==std::cout) obj.Show();
     else out << "cannot print pf_gammajettree\n";
     return out;
   }

};

#endif

#ifdef misc_gammajettree_cxx
misc_gammajettree::misc_gammajettree(TString fname) :
  fChain(new TChain("miscItems/misc_tree"))
{
  if (fname.Length()==0) return;
  if (!this->Init(fname)) {
    std::cout << "Initialization failed in constructor" << std::endl;
  }
}

misc_gammajettree::~misc_gammajettree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t misc_gammajettree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t misc_gammajettree::LoadTree(Long64_t entry)
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

int misc_gammajettree::Init(const TString &fname)
{
  if (fname.Length()==0) {
    std::cout << "misc_gammajettree::Init: non-empty fname is expected\n";
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
   photonTriggerNames = 0;
   jetTriggerNames = 0;

   // Clear variables
   ignoreHLT=doCaloJets=doPFJets=doGenJets=workOnAOD=false;
   nProcessed=0;

   // Set branch addresses and branch pointers
   fCurrent = -1;
   //fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ignoreHLT", &ignoreHLT, &b_ignoreHLT);
   fChain->SetBranchAddress("doCaloJets", &doCaloJets, &b_doCaloJets);
   fChain->SetBranchAddress("doPFJets", &doPFJets, &b_doPFJets);
   fChain->SetBranchAddress("doGenJets", &doGenJets, &b_doGenJets);
   fChain->SetBranchAddress("workOnAOD", &workOnAOD, &b_workOnAOD);
   fChain->SetBranchAddress("photonTriggerNames", &photonTriggerNames, &b_photonTriggerNames);
   fChain->SetBranchAddress("jetTriggerNames", &jetTriggerNames, &b_jetTriggerNames);
   fChain->SetBranchAddress("nProcessed", &nProcessed, &b_nProcessed);
   Notify();
   return 1;
}

Bool_t misc_gammajettree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void misc_gammajettree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t misc_gammajettree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  if (0) std::cout << "entry=" << entry << "\n"; // dummy statement
   return 1;
}
#endif // #ifdef misc_gammajettree_cxx
