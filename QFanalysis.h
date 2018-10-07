//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct  5 16:43:24 2015 by ROOT version 5.34/30
// from TTree h2000/TileCal-Ntuple
// found on file: small61.root
//////////////////////////////////////////////////////////

#ifndef QFanalysis_h
#define QFanalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class QFanalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Short_t         sample[4][64][48][7];
   Float_t         ene[4][64][48];
   Float_t         time[4][64][48];
   Float_t         chi2[4][64][48];
   Float_t         eMF[4][64][48][7];
   Float_t         tMF[4][64][48][7];
   Float_t         chi2MF[4][64][48];
   Float_t         pedMF[4][64][48];

   // List of branches
   TBranch        *b_sample;   //!
   TBranch        *b_ene;   //!
   TBranch        *b_time;   //!
   TBranch        *b_chi2;   //!
   TBranch        *b_eMF;   //!
   TBranch        *b_tMF;   //!
   TBranch        *b_chi2MF;   //!
   TBranch        *b_pedMF;   //!

   QFanalysis(TTree *tree=0);
   virtual ~QFanalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef QFanalysis_cxx
QFanalysis::QFanalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("small_all.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("small_all.root");
      }
      f->GetObject("h2000",tree);

   }
   Init(tree);
}

QFanalysis::~QFanalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t QFanalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t QFanalysis::LoadTree(Long64_t entry)
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

void QFanalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("sample", sample, &b_sample);
   fChain->SetBranchAddress("ene", ene, &b_ene);
   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("chi2", chi2, &b_chi2);
   fChain->SetBranchAddress("eMF", eMF, &b_eMF);
   fChain->SetBranchAddress("tMF", tMF, &b_tMF);
   fChain->SetBranchAddress("chi2MF", chi2MF, &b_chi2MF);
   fChain->SetBranchAddress("pedMF", pedMF, &b_pedMF);
   Notify();
}

Bool_t QFanalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void QFanalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t QFanalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef QFanalysis_cxx
