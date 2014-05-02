//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 30 20:07:26 2014 by ROOT version 5.34/11
// from TTree tree_passedEvents/tree_passedEvents
// found on file: PosAn_BTF_92_20140430-020137_beam.root
//////////////////////////////////////////////////////////

#ifndef hodoScan_h
#define hodoScan_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxhodox = 1;
const Int_t kMaxhodoy = 1;
const Int_t kMaxbgo = 1;
const Int_t kMaxcef3 = 1;
const Int_t kMaxbgo_corr = 1;
const Int_t kMaxcef3_corr = 1;
const Int_t kMaxhodox_corr = 1;
const Int_t kMaxhodoy_corr = 1;

class hodoScan {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         evtNumber;
   UInt_t          adcData;
   UInt_t          adcBoard;
   UInt_t          adcChannel;
   Int_t           nHodoFibersX;
   Int_t           nHodoFibersY;
   Int_t           nHodoFibersCorrX;
   Int_t           nHodoFibersCorrY;
   Int_t           hodox_chan;
   Int_t           hodoy_chan;
   Int_t           cef3_chan;
   Int_t           bgo_chan;
   Float_t         hodox[8];   //[hodox_chan]
   Float_t         hodoy[8];   //[hodoy_chan]
   Float_t         bgo[8];   //[bgo_chan]
   Float_t         cef3[4];   //[cef3_chan]
   Float_t         bgo_corr[8];   //[bgo_chan]
   Float_t         cef3_corr[4];   //[cef3_chan]
   Float_t         hodox_corr[8];   //[hodox_chan]
   Float_t         hodoy_corr[8];   //[hodoy_chan]

   // List of branches
   TBranch        *b_evtNumber;   //!
   TBranch        *b_adcData;   //!
   TBranch        *b_adcBoard;   //!
   TBranch        *b_adcChannel;   //!
   TBranch        *b_nHodoFibersX;   //!
   TBranch        *b_nHodoFibersY;   //!
   TBranch        *b_nHodoFibersCorrX;   //!
   TBranch        *b_nHodoFibersCorrY;   //!
   TBranch        *b_hodox_chan;   //!
   TBranch        *b_hodoy_chan;   //!
   TBranch        *b_cef3_chan;   //!
   TBranch        *b_bgo_chan;   //!
   TBranch        *b_hodox;   //!
   TBranch        *b_hodoy;   //!
   TBranch        *b_bgo;   //!
   TBranch        *b_cef3;   //!
   TBranch        *b_bgo_corr;   //!
   TBranch        *b_cef3_corr;   //!
   TBranch        *b_hodox_corr;   //!
   TBranch        *b_hodoy_corr;   //!

   hodoScan(TTree *tree=0);
   virtual ~hodoScan();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef hodoScan_cxx
hodoScan::hodoScan(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("PosAn_BTF_92_20140430-020137_beam.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("PosAn_BTF_92_20140430-020137_beam.root");
      }
      f->GetObject("tree_passedEvents",tree);

   }
   Init(tree);
}

hodoScan::~hodoScan()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t hodoScan::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t hodoScan::LoadTree(Long64_t entry)
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

void hodoScan::Init(TTree *tree)
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

   fChain->SetBranchAddress("evtNumber", &evtNumber, &b_evtNumber);
   fChain->SetBranchAddress("adcData", &adcData, &b_adcData);
   fChain->SetBranchAddress("adcBoard", &adcBoard, &b_adcBoard);
   fChain->SetBranchAddress("adcChannel", &adcChannel, &b_adcChannel);
   fChain->SetBranchAddress("nHodoFibersX", &nHodoFibersX, &b_nHodoFibersX);
   fChain->SetBranchAddress("nHodoFibersY", &nHodoFibersY, &b_nHodoFibersY);
   fChain->SetBranchAddress("nHodoFibersCorrX", &nHodoFibersCorrX, &b_nHodoFibersCorrX);
   fChain->SetBranchAddress("nHodoFibersCorrY", &nHodoFibersCorrY, &b_nHodoFibersCorrY);
   fChain->SetBranchAddress("hodox_chan", &hodox_chan, &b_hodox_chan);
   fChain->SetBranchAddress("hodoy_chan", &hodoy_chan, &b_hodoy_chan);
   fChain->SetBranchAddress("cef3_chan", &cef3_chan, &b_cef3_chan);
   fChain->SetBranchAddress("bgo_chan", &bgo_chan, &b_bgo_chan);
   fChain->SetBranchAddress("hodox", hodox, &b_hodox);
   fChain->SetBranchAddress("hodoy", hodoy, &b_hodoy);
   fChain->SetBranchAddress("bgo", bgo, &b_bgo);
   fChain->SetBranchAddress("cef3", cef3, &b_cef3);
   fChain->SetBranchAddress("bgo_corr", bgo_corr, &b_bgo_corr);
   fChain->SetBranchAddress("cef3_corr", cef3_corr, &b_cef3_corr);
   fChain->SetBranchAddress("hodox_corr", hodox_corr, &b_hodox_corr);
   fChain->SetBranchAddress("hodoy_corr", hodoy_corr, &b_hodoy_corr);
   Notify();
}

Bool_t hodoScan::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void hodoScan::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t hodoScan::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef hodoScan_cxx
