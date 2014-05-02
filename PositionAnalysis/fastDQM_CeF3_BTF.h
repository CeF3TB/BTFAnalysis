//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 20 23:56:33 2014 by ROOT version 5.34/12
// from TTree eventRawData/VME event RAW Data
// found on file: run0000021_20140220.root
//////////////////////////////////////////////////////////

#ifndef fastDQM_CeF3_BTF_h
#define fastDQM_CeF3_BTF_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>

#define NUM_ADC_CHANNELS 40

// ADC BOARD 0 -> V265
// ADC BOARD 1 -> V792

#define BGO_ADC_BOARD 0
#define BGO_CHANNELS 8
#define BGO_ADC_START_CHANNEL 0
#define BGO_ADC_END_CHANNEL 7

#define CEF3_ADC_BOARD 1
#define CEF3_CHANNELS 4
#define CEF3_ADC_START_CHANNEL 2
#define CEF3_ADC_END_CHANNEL 5

#define HODOX_ADC_BOARD 1
#define HODOX_CHANNELS 8
#define HODOX_ADC_START_CHANNEL 16
#define HODOX_ADC_END_CHANNEL 23

#define HODOY_ADC_BOARD 1
#define HODOY_CHANNELS 8
#define HODOY_ADC_START_CHANNEL 24
#define HODOY_ADC_END_CHANNEL 31

#define CENTERX_TAGGER_ADC_BOARD 1
#define CENTERX_TAGGER_CHANNELS 2 
#define CENTERX_TAGGER_ADC_START_CHANNEL 7
#define CENTERX_TAGGER_ADC_END_CHANNEL 8

#define CENTERY_TAGGER_ADC_BOARD 1
#define CENTERY_TAGGER_CHANNELS 1 //One channel is missing
#define CENTERY_TAGGER_ADC_START_CHANNEL 6
#define CENTERY_TAGGER_ADC_END_CHANNEL 6

#define SCINT_FRONT_ADC_BOARD 1
#define SCINT_FRONT_ADC_CHANNEL 0

#define SCINT_BACK_ADC_BOARD 1
#define SCINT_BACK_ADC_CHANNEL 9

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class fastDQM_CeF3_BTF {
public :
  std::map<TString,TObject*> outObjects;
  std::string          outFile;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          evtNumber;
   UInt_t          evtTimeDist;
   UInt_t          evtTimeStart;
   UInt_t          evtTime;
   UInt_t          boardTriggerBit;
   UInt_t          triggerWord;
   UInt_t          nAdcChannels;
   UInt_t          adcBoard[NUM_ADC_CHANNELS];   //[nAdcChannels]
   UInt_t          adcChannel[NUM_ADC_CHANNELS];   //[nAdcChannels]
   UInt_t          adcData[NUM_ADC_CHANNELS];   //[nAdcChannels]
   UInt_t          nTdcChannels;
   UInt_t          tdcBoard[1];   //[nTdcChannels]
   UInt_t          tdcChannel[1];   //[nTdcChannels]
   UInt_t          tdcData[1];   //[nTdcChannels]

   // List of branches
   TBranch        *b_evtNumber;   //!
   TBranch        *b_evtTimeDist;   //!
   TBranch        *b_evtTimeStart;   //!
   TBranch        *b_evtTime;   //!
   TBranch        *b_boardTriggerBit;   //!
   TBranch        *b_triggerWord;   //!
   TBranch        *b_nAdcChannels;   //!
   TBranch        *b_adcBoard;   //!
   TBranch        *b_adcChannel;   //!
   TBranch        *b_adcData;   //!
   TBranch        *b_nTdcChannels;   //!
   TBranch        *b_tdcBoard;   //!
   TBranch        *b_tdcChannel;   //!
   TBranch        *b_tdcData;   //!

   fastDQM_CeF3_BTF(TTree *tree=0);
   virtual ~fastDQM_CeF3_BTF();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef fastDQM_CeF3_BTF_cxx
fastDQM_CeF3_BTF::fastDQM_CeF3_BTF(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("run0000021_20140220.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("run0000021_20140220.root");
      }
      f->GetObject("eventRawData",tree);

   }
   Init(tree);
}

fastDQM_CeF3_BTF::~fastDQM_CeF3_BTF()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fastDQM_CeF3_BTF::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fastDQM_CeF3_BTF::LoadTree(Long64_t entry)
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

void fastDQM_CeF3_BTF::Init(TTree *tree)
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
   fChain->SetBranchAddress("evtTimeDist", &evtTimeDist, &b_evtTimeDist);
   fChain->SetBranchAddress("evtTimeStart", &evtTimeStart, &b_evtTimeStart);
   fChain->SetBranchAddress("evtTime", &evtTime, &b_evtTime);
   fChain->SetBranchAddress("boardTriggerBit", &boardTriggerBit, &b_boardTriggerBit);
   fChain->SetBranchAddress("triggerWord", &triggerWord, &b_triggerWord);
   fChain->SetBranchAddress("nAdcChannels", &nAdcChannels, &b_nAdcChannels);
   fChain->SetBranchAddress("adcBoard", adcBoard, &b_adcBoard);
   fChain->SetBranchAddress("adcChannel", adcChannel, &b_adcChannel);
   fChain->SetBranchAddress("adcData", adcData, &b_adcData);
   fChain->SetBranchAddress("nTdcChannels", &nTdcChannels, &b_nTdcChannels);
   fChain->SetBranchAddress("tdcBoard", &tdcBoard, &b_tdcBoard);
   fChain->SetBranchAddress("tdcChannel", &tdcChannel, &b_tdcChannel);
   fChain->SetBranchAddress("tdcData", &tdcData, &b_tdcData);
   Notify();
}

Bool_t fastDQM_CeF3_BTF::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void fastDQM_CeF3_BTF::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t fastDQM_CeF3_BTF::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef fastDQM_CeF3_BTF_cxx
