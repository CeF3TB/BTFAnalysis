#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"

#include "TMVA/Reader.h"

#include "fastDQM_CeF3_BTF.h"
#include "interface/RunHelper.h"
#include "interface/PositionTools.h"
#include "interface/DrawTools.h"
#include "TGAxis.h"
#include "TLegend.h"
#include "TPaletteAxis.h"

int main( int argc, char* argv[] ) {


  std::string runName = "precalib_BGO_pedestal_noSource";
  if( argc>1 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
  }

  std::string tag = "V00";
  if( argc>2 ) {
    std::string tag_str(argv[2]);
    tag = tag_str;
  }

  std::string cuts="";
  if( argc>3 ) {
    std::string cuts_str(argv[3]);
    cuts = cuts_str;
  }
  std::cout<<"cuts="<<cuts<<std::endl;

  TString cuts_tstr(cuts);

  TString runName_tstr(runName);
  bool isOnlyRunNumber = !(runName_tstr.BeginsWith("BTF_"));


  TChain* tree = new TChain("recoTree");
  if( isOnlyRunNumber ) {
    std::cout << "-> We believe you are passing the program only the run number!" << std::endl;
    std::cout << "-> So for instance you are passing '246' for run 'BTF_246_20140501-212512_beam'" << std::endl;
    std::cout << "(if this is not the case this means TROUBLE)" << std::endl;
    std::cout << "-> Will look for runs matching run number: " << runName << std::endl;
    tree->Add(Form("analysisTrees_%s/Reco_BTF_%s_*beam.root/recoTree", tag.c_str(), runName.c_str()) );
    if( tree->GetEntries()==0 ) {
      std::cout << "WARNING! Didn't find any events matching run: " << runName << std::endl;
      std::cout << "Exiting" << std::endl;
      exit(1913);
    }
  } else {
    std::string fileName = "analysisTrees_"+tag+"/Reco_" + runName + ".root";
    TFile* file = TFile::Open(fileName.c_str());
    if( file==0 ) {
      std::cout << "ERROR! Din't find file " << fileName << std::endl;
      std::cout << "Exiting." << std::endl;
      exit(11);
    }
    tree = (TChain*)file->Get("recoTree");
  }



   UInt_t          runNumber;
   UInt_t          evtNumber;
   UInt_t          adcData;
   UInt_t          adcBoard;
   UInt_t          adcChannel;
   Int_t           nHodoFibersX;
   Int_t           nHodoFibersY;
   Int_t           nHodoClustersX;
   Float_t         pos_hodoClustX[CEF3_CHANNELS];   //[nHodoClustersX]
   Int_t           nFibres_hodoClustX[CEF3_CHANNELS];   //[nHodoClustersX]
   Int_t           nHodoClustersY;
   Float_t         pos_hodoClustY[CEF3_CHANNELS];   //[nHodoClustersY]
   Int_t           nFibres_hodoClustY[CEF3_CHANNELS];   //[nHodoClustersY]
   Int_t           hodox_chan;
   Int_t           hodoy_chan;
   Int_t           cef3_chan;
   Int_t           bgo_chan;
   Float_t         hodox[BGO_CHANNELS];   //[hodox_chan]
   Float_t         hodoy[BGO_CHANNELS];   //[hodoy_chan]
   Float_t         bgo[BGO_CHANNELS];   //[bgo_chan]
   Float_t         cef3[CEF3_CHANNELS];   //[cef3_chan]
   Float_t         cef3_pedSubtracted[CEF3_CHANNELS];   //[cef3_chan]
   Float_t         bgo_corr[BGO_CHANNELS];   //[bgo_chan]
   Float_t         bgo_pedSubtracted[BGO_CHANNELS];   //[bgo_chan]
   Float_t         cef3_corr[CEF3_CHANNELS];   //[cef3_chan]
   Float_t         hodox_corr[BGO_CHANNELS];   //[hodox_chan]
   Float_t         hodoy_corr[BGO_CHANNELS];   //[hodoy_chan]
   Float_t         scintFront;
   Bool_t          isSingleEle_scintFront;
   Float_t         xBeam;
   Float_t         yBeam;
   Bool_t          cef3_ok;
   Bool_t          cef3_corr_ok;
   Bool_t          bgo_ok;
   Bool_t          bgo_corr_ok;

   tree->SetBranchAddress("runNumber", &runNumber);
   tree->SetBranchAddress("evtNumber", &evtNumber);
   tree->SetBranchAddress("adcData", &adcData);
   tree->SetBranchAddress("adcBoard", &adcBoard);
   tree->SetBranchAddress("adcChannel", &adcChannel);
   tree->SetBranchAddress("nHodoFibersX", &nHodoFibersX);
   tree->SetBranchAddress("nHodoFibersY", &nHodoFibersY);
   tree->SetBranchAddress("nHodoClustersX", &nHodoClustersX);
   tree->SetBranchAddress("pos_hodoClustX", pos_hodoClustX);
   tree->SetBranchAddress("nFibres_hodoClustX", nFibres_hodoClustX);
   tree->SetBranchAddress("nHodoClustersY", &nHodoClustersY);
   tree->SetBranchAddress("pos_hodoClustY", pos_hodoClustY);
   tree->SetBranchAddress("nFibres_hodoClustY", nFibres_hodoClustY);
   tree->SetBranchAddress("hodox_chan", &hodox_chan);
   tree->SetBranchAddress("hodoy_chan", &hodoy_chan);
   tree->SetBranchAddress("cef3_chan", &cef3_chan);
   tree->SetBranchAddress("bgo_chan", &bgo_chan);
   tree->SetBranchAddress("hodox", hodox);
   tree->SetBranchAddress("hodoy", hodoy);
   tree->SetBranchAddress("bgo", bgo);
   tree->SetBranchAddress("cef3", cef3);
   tree->SetBranchAddress("cef3_pedSubtracted", cef3_pedSubtracted);
   tree->SetBranchAddress("bgo_corr", bgo_corr);
   tree->SetBranchAddress("bgo_pedSubtracted", bgo_pedSubtracted);
   tree->SetBranchAddress("cef3_corr", cef3_corr);
   tree->SetBranchAddress("hodox_corr", hodox_corr);
   tree->SetBranchAddress("hodoy_corr", hodoy_corr);
   tree->SetBranchAddress("scintFront", &scintFront);
   tree->SetBranchAddress("isSingleEle_scintFront", &isSingleEle_scintFront);
   tree->SetBranchAddress("xBeam", &xBeam);
   tree->SetBranchAddress("yBeam", &yBeam);
   tree->SetBranchAddress("cef3_ok", &cef3_ok);
   tree->SetBranchAddress("cef3_corr_ok", &cef3_corr_ok);
   tree->SetBranchAddress("bgo_ok", &bgo_ok);
   tree->SetBranchAddress("bgo_corr_ok", &bgo_corr_ok);


   int nentries = tree->GetEntries();

  if( isOnlyRunNumber ) {
    // modify runname in such a way that it's useful for getBeamPosition and outfile:
    runName = "BTF_" + runName + "_beam";
  }

  std::string outputdir = "2DAnTrees_"+tag;

  system( Form("mkdir -p %s", outputdir.c_str()) );
  std::string outfileName = outputdir + "/2DAn_" + runName + ".root";
  TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );

  TStyle* style = DrawTools::setStyle();
  //  tree->Print();

  TH1F* cef3SpectrumTotal=new TH1F("cef3SpectrumTotal","cef3SpectrumTotal",300,0,11999/1000.);
  TH1F* cef3SpectrumSingleEle=new TH1F("cef3SpectrumSingleEle","cef3SpectrumSingleEle",300,0,11999/1000.);
  TH1F* cef3SpectrumSingleEleHodo=new TH1F("cef3SpectrumSingleEleHodo","cef3SpectrumSingleEleHodo",300,0,11999/1000.);

  cef3SpectrumTotal->SetLineWidth(2);
  cef3SpectrumSingleEle->SetLineWidth(2);
  cef3SpectrumSingleEleHodo->SetLineWidth(2);

  //  TH2F* cef3VsScintFront= new TH2F("cef3VsScintFront","cef3VsScintFront",40,0.,4000./1000.,40,0.,12000.);
  TH2F* cef3VsScintFront= new TH2F("cef3VsScintFront","cef3VsScintFront",100,0.,4000./1000.,100,0.,12000./1000.);
  TH2F* cef3TopVsBottom= new TH2F("cefTopVsBottom","cefTopVsBottom",100,0.,8000./1000.,100,0.,8000./1000.);


  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {
    tree->GetEntry(iEntry);
    if( iEntry % 10000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;
    if(cef3_corr_ok){
      cef3SpectrumTotal->Fill((cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3])/1000.);
      if(isSingleEle_scintFront)      cef3SpectrumSingleEle->Fill((cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3])/1000.);
      if(isSingleEle_scintFront*(nHodoClustersX==1 && nHodoClustersY==1))      cef3SpectrumSingleEleHodo->Fill((cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3])/1000.);
    }

      if(cuts=="hodo"){
	if(!(nHodoClustersX>=1 && nHodoClustersY>=1))continue;
      }else if (cuts=="signalFibers"){
	if(!(cef3_corr[0]>10 || cef3_corr[1]>10 || cef3_corr[2]>10 || cef3_corr[3]>0))continue;
      }else if (cuts=="bgo"){
	if(!((bgo_corr[0]+bgo_corr[1]+bgo_corr[2]+bgo_corr[3]+bgo_corr[4]+bgo_corr[5]+bgo_corr[6]+bgo_corr[7])<1000))continue;
      }

    if(cuts!=""){if(scintFront<150.) continue;}
    if(scintFront<0.5 ||(cef3[0]+cef3[1]+cef3[2]+cef3[3])<1.5*1000. ) continue;
    cef3VsScintFront->Fill(scintFront/1000.,(cef3[0]+cef3[1]+cef3[2]+cef3[3])/1000.);
    cef3TopVsBottom->Fill((cef3_corr[2]+cef3_corr[3])/1000.,(cef3_corr[0]+cef3_corr[1])/1000.);
  }

  style->SetPadTopMargin(0.05);
  style->SetPadBottomMargin(0.13);//0.13);                                                                                                                   
  style->SetPadLeftMargin(0.17);//0.16);                                                                                                                     
  style->SetPadRightMargin(0.13);//0.02); 

  TCanvas* c1= new TCanvas("c1"," ", 600,600);
  c1->cd();



  cef3VsScintFront->SetYTitle("CeF_{3} [ADC Channel / 1000]");
  //  cef3VsScintFront->SetTitleOffset(1.7,"Y");
  
  ((TGaxis*)cef3VsScintFront->GetXaxis())->SetMaxDigits(3);

  //  cef3VsScintFront->SetAxisRange(10,cef3VsScintFront->GetMaximum(),"Z");
  //  cef3VsScintFront->SetXTitle("Front Scintillator [ADC Channel / 1000]");
  TPaveText* labelXaxis= new TPaveText(0.1,0.05,0.965,0., "brNDC"); 
  labelXaxis->SetFillColor(kWhite);
  labelXaxis->SetTextSize(cef3VsScintFront->GetXaxis()->GetTitleSize());
  labelXaxis->SetTextAlign(31); // align right
  labelXaxis->SetTextFont(cef3VsScintFront->GetXaxis()->GetTitleFont());
  labelXaxis->AddText("Front Scintillator [ADC Channel / 1000]");
  //  cef3VsScintFront->SetAxisRange(9.99,10000.,"Z");

//  cef3VsScintFront->SetAxisRange(0.5,3.9,"X");
//  cef3VsScintFront->SetAxisRange(1.5*1000,10.*1000,"Y");
  cef3VsScintFront->Draw("colz");  

  labelXaxis->Draw("same");
  ((TGaxis*)cef3VsScintFront->GetYaxis())->Draw("same");
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)cef3VsScintFront->GetListOfFunctions()->FindObject("palette");
//  palette->SetX1NDC(0.88);
//  palette->SetX2NDC(0.92);
//  palette->SetY1NDC(0.2);
//  palette->SetY2NDC(0.95);

  gPad->Modified();
  gPad->Update();


  TPaveText* labelTop = DrawTools::getLabelTop_2D();
  labelTop->Draw("same");
  //   c1->SetLogz();
  c1->SaveAs("cef3VsScintFront_"+cuts_tstr+".C");  
  c1->SaveAs("cef3VsScintFront_"+cuts_tstr+".png");
  c1->SaveAs("cef3VsScintFront_"+cuts_tstr+".eps");
  c1->SaveAs("cef3VsScintFront_"+cuts_tstr+".pdf");

  c1->Clear();
  c1->cd();
  cef3TopVsBottom->SetYTitle("Top Fibres [ADC Channel / 1000]");
  //  cef3TopVsBottom->SetTitleOffset(1.7,"Y");
  //  cef3TopVsBottom->SetXTitle("CeF_{3} Bottom Fibres [ADC Channel / 1000]");
  //  cef3TopVsBottom->SetAxisRange(9.99,10000.,"Z");
  cef3TopVsBottom->Draw("colz");  
  TPaveText* labelXaxisTop= new TPaveText(0.1,0.05,0.965,0., "brNDC"); 
  labelXaxisTop->SetFillColor(kWhite);
  labelXaxisTop->SetTextSize(cef3VsScintFront->GetXaxis()->GetTitleSize());
  labelXaxisTop->SetTextAlign(31); // align right
  labelXaxisTop->SetTextFont(cef3VsScintFront->GetXaxis()->GetTitleFont());
  labelXaxisTop->AddText("Bottom Fibres [ADC Channel / 1000]");
  labelXaxisTop->Draw("same");

  labelTop->Draw("same");

  TPaveText* additional_palette=DrawTools::getCef3LabelLeft();
  additional_palette->Draw("same");

  gPad->Update();
  //  TPaletteAxis *palette2 = (TPaletteAxis*)cef3TopVsBottom->GetListOfFunctions()->FindObject("palette");
//  palette2->SetX1NDC(0.88);
//  palette2->SetX2NDC(0.92);
//  palette2->SetY1NDC(0.13);
//  palette2->SetY2NDC(0.95);

  gPad->Modified();
  gPad->Update();


  //  c1->SetLogz();
  c1->SaveAs("cef3TopVsBottom_"+cuts_tstr+".C");  
  c1->SaveAs("cef3TopVsBottom_"+cuts_tstr+".png");
  c1->SaveAs("cef3TopVsBottom_"+cuts_tstr+".eps");
  c1->SaveAs("cef3TopVsBottom_"+cuts_tstr+".pdf");

  c1->Clear();
  c1->cd();
  cef3SpectrumTotal->Draw();
  cef3SpectrumTotal->SetXTitle("ADC Channel / 1000");
  cef3SpectrumTotal->SetYTitle("Entries / (40 ADC Channels)");
  //  TPaveText *l1 = new  TPaveText(0.1,0.1,0.3,0.3,"brNDC");

  bool drawAuthors=false;
  TText *l1 = new TText(12.8,1.0,"P. Meridiani et al., to be publ. in Proc. IEEE NSS 2014");
  //  l1->AddText(
  //  l1->SetFillColor(kWhite);
  l1->SetTextSize(0.035);
  //  l1->SetTextAlign(11); // align right
  l1->SetTextFont(42);
  
  l1->SetTextAngle(90);
  //  l1->SetTextAlign(22);
  if(drawAuthors)  l1->Draw("same");
  TPaveText* labelTop2 = DrawTools::getLabelTop_2D();
  labelTop2->Draw("same");
  TLegend* legend = new TLegend(0.50, 0.80, 0.85, 0.85,"W-CeF_{3} Single Tower");
  legend->SetTextSize(0.036);
  legend->SetFillColor(kWhite);
  legend->SetLineColor(kWhite);
  legend->SetFillStyle(0);
  legend->Draw("same");
  c1->SaveAs("cef3SpectrumTotal.C");
  c1->SaveAs("cef3SpectrumTotal.png");
  c1->SaveAs("cef3SpectrumTotal.eps");
  c1->SaveAs("cef3SpectrumTotal.pdf");
  c1->SetLogy();
  c1->SaveAs("cef3SpectrumTotal_log.C");
  c1->SaveAs("cef3SpectrumTotal_log.png");
  c1->SaveAs("cef3SpectrumTotal_log.eps");
  c1->SaveAs("cef3SpectrumTotal_log.pdf");
  c1->SetLogy(0);



  c1->Clear();
  c1->cd();
  cef3SpectrumSingleEle->Draw();
  cef3SpectrumSingleEle->SetXTitle("ADC Channel / 1000");
  cef3SpectrumSingleEle->SetYTitle("Entries / (40 ADC Channels)");
  legend->Draw("same");
  labelTop2->Draw("same");
  c1->SaveAs("cef3SpectrumSingleEle.png");
  c1->SaveAs("cef3SpectrumSingleEle.eps");
  c1->SaveAs("cef3SpectrumSingleEle.pdf");
  c1->SetLogy();
  c1->SaveAs("cef3SpectrumSingleEle_log.png");
  c1->SaveAs("cef3SpectrumSingleEle_log.eps");
  c1->SaveAs("cef3SpectrumSingleEle_log.pdf");

  c1->Clear();
  c1->SetLogy(0);
  c1->cd();
  cef3SpectrumSingleEleHodo->Draw();
  legend->Draw("same");
  cef3SpectrumSingleEleHodo->SetXTitle("ADC Channel / 1000");
  cef3SpectrumSingleEleHodo->SetYTitle("Entries / (40 ADC Channels)");
  labelTop2->Draw("same");
  c1->SaveAs("cef3SpectrumSingleEleHodo.png");
  c1->SaveAs("cef3SpectrumSingleEleHodo.eps");
  c1->SaveAs("cef3SpectrumSingleEleHodo.pdf");
  c1->SetLogy();
  c1->SaveAs("cef3SpectrumSingleEleHodo_log.png");
  c1->SaveAs("cef3SpectrumSingleEleHodo_log.eps");
  c1->SaveAs("cef3SpectrumSingleEleHodo_log.pdf");
  c1->SetLogy(0);

  c1->Clear();


  TStyle* style2 = DrawTools::setStyle();
  TCanvas* c2= new TCanvas("c2"," ", 600,600);
  c2->cd();

  TLegend* legend2 = new TLegend(0.50, 0.60, 0.85, 0.85,"W-CeF_{3} Single Tower");
  legend2->SetTextSize(0.034);
  legend2->SetFillColor(kWhite);
  legend2->SetLineColor(0);
  legend2->SetFillStyle(0);

  TLegend* legend3 = new TLegend(0.58, 0.72, 0.85, 0.92,"W-CeF_{3} Single Tower");
  legend3->SetTextSize(0.029);
  legend3->SetFillColor(kWhite);
  legend3->SetLineColor(0);
  legend3->SetFillStyle(0);



  cef3SpectrumTotal->Draw();
  cef3SpectrumSingleEle->SetLineColor(kRed);
  cef3SpectrumSingleEle->Draw("same");
  cef3SpectrumSingleEleHodo->SetLineColor(kBlue);
  legend2->AddEntry(cef3SpectrumTotal,"All events","l");
  legend2->AddEntry(cef3SpectrumSingleEle,"Single e^{-} selection","l");
  legend2->AddEntry(cef3SpectrumSingleEleHodo,"Central 8x8 mm^{2}","l");
  if(drawAuthors)  l1->Draw("same");
  legend2->Draw("same");
  labelTop2->Draw("same");
  cef3SpectrumSingleEleHodo->Draw("same");
  gPad->RedrawAxis();
//  c2->SaveAs("cef3SpectrumSuperimposed.C");
//  c2->SaveAs("cef3SpectrumSuperimposed.png");
//  c2->SaveAs("cef3SpectrumSuperimposed.eps");
//  c2->SaveAs("cef3SpectrumSuperimposed.pdf");
  c2->SetLogy();
  cef3SpectrumTotal->GetYaxis()->SetRangeUser(1,5.*cef3SpectrumTotal->GetMaximum());
  legend2->Clear();
  legend3->AddEntry(cef3SpectrumTotal,"All events","l");
  legend3->AddEntry(cef3SpectrumSingleEle,"Single e^{-} selection","l");
  legend3->AddEntry(cef3SpectrumSingleEleHodo,"Central 8x8 mm^{2}","l");
  legend3->Draw("same");
  cef3SpectrumTotal->Draw("same");
  if(drawAuthors)  l1->Draw("same");
  c2->SaveAs("cef3SpectrumSuperimposed_log.C");
  c2->SaveAs("cef3SpectrumSuperimposed_log.png");
  c2->SaveAs("cef3SpectrumSuperimposed_log.eps");
  c2->SaveAs("cef3SpectrumSuperimposed_log.pdf");
  c2->SetLogy(0);


}
