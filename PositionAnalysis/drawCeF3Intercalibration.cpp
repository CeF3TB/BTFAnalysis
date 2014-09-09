#include <iostream>
#include <string>
#include <stdlib.h>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TEllipse.h"
#include "TString.h"

#include "interface/DrawTools.h"

#include "TApplication.h"


int main ( int argc, char* argv[] ) {

  TApplication* a = new TApplication("a", 0, 0);

  TStyle* style = DrawTools::setStyle();
  style->cd();

  std::string inputDir = "./analysisTrees";
  std::string runName = "BTF_259_20140502-012847_beam";
  std::string tag = "default";



  if( argc == 3 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
    std::string tag_str(argv[2]);
    tag = tag_str;
  }else if (argc==4){
    std::string inputDir_str(argv[2]);
    inputDir =inputDir_str;
  } else{
    std::cout<<"Usage:"<<std::endl;
    std::cout<<"./calibrateCef3 BTF_XX tag "<<std::endl;
    exit(12345);
  }

  if(tag!="default") inputDir = inputDir + "_"+tag;


  std::string outputdir = "CeF3IntercalPlots/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );



  TFile* file = TFile::Open(Form("%s/Reco_%s.root", inputDir.c_str(),runName.c_str()));
  std::cout<<"opening file:"<<file->GetName();

  TTree* tree = (TTree*)file->Get("recoTree");

  TCanvas* canny = new TCanvas("canny", "CeF3 Intercalibration", 1200,800);
  canny->Divide(2,1);
  canny->cd(1);

  int nBins = 75;
  float xMin = 200.;
  float xMax = 2500.;

  TH1D* hu0 = new TH1D("hu0","",nBins,xMin,xMax);
  hu0->SetLineColor(kBlue);
  TH1D* hu1 = new TH1D("hu1","",nBins,xMin,xMax);
 hu1->SetLineColor(kAzure+8);
  TH1D* hu2 = new TH1D("hu2","",nBins,xMin,xMax);
  hu2->SetLineColor(kGreen+1);
  TH1D* hu3 = new TH1D("hu3","",nBins,xMin,xMax);
  hu3->SetLineColor(kMagenta-4);


 
  TH1D* h0 = new TH1D("h0","",nBins,xMin,xMax);
 h0->SetLineColor(kBlue);
  TH1D* h1 = new TH1D("h1","",nBins,xMin,xMax);
  h1->SetLineColor(kAzure+8);
  TH1D* h2 = new TH1D("h2","",nBins,xMin,xMax);  
  h2->SetLineColor(kGreen+1);
  TH1D* h3 = new TH1D("h3","",nBins,xMin,xMax);
  h3->SetLineColor(kMagenta-4);


  tree->Draw("cef3[0]>>hu0");
  tree->Draw("cef3[1]>>hu1","","same");
  tree->Draw("cef3[2]>>hu2","","same");
  tree->Draw("cef3[3]>>hu3","","same");

  canny->cd(2);

  tree->Draw("cef3_corr[0]>>h0");
  tree->Draw("cef3_corr[1]>>h1","","same");
  tree->Draw("cef3_corr[2]>>h2","","same");
  tree->Draw("cef3_corr[3]>>h3","","same");

  TLegend* legend2 = new TLegend( 0.6, 0.9, 0.9, 0.6 );
  legend2->SetFillStyle(0);
  legend2->SetLineColor(0);
  legend2->SetLineWidth(0);
  legend2->SetTextSize(0.035);
  legend2->AddEntry( h0, "CeF3_corr[0]" );
  legend2->AddEntry( h1, "CeF3_corr[1]" );
  legend2->AddEntry( h2, "CeF3_corr[2]" );
  legend2->AddEntry( h3, "CeF3_corr[3]" );
  legend2->Draw("same");


//////////////////////////////////////////////////////////////////////////////
//only central events 

  TCanvas* canny2 = new TCanvas("canny2","",1200,800);
  canny2->cd();


  TH1D* FibX=new TH1D("FibX","",6,0,6);
  FibX->SetLineColor(kAzure+8);
  TH1D* FibY=new TH1D("FibY","",6,0,6);
  FibY->SetLineColor(kBlue);

  tree->Draw("nFibres_hodoClustX>>FibX","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1)");
  tree->Draw("nFibres_hodoClustY>>FibY","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1)","same");

  FibX->GetXaxis()->SetTitle("nFibres");

  TLegend* legend = new TLegend( 0.6, 0.9, 0.9, 0.6 );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
  legend->SetLineWidth(0);
  legend->SetTextSize(0.035);
  legend->AddEntry("FibX", "nFibres_hodoClustX" );
  legend->AddEntry("FibY" , "nFibres_hodoClustY" );
  legend->Draw("same");


  /////////////////////////////////////////////////
  //POSITION of the cluster and cut to it to get only central events

  TCanvas* canny3 = new TCanvas("canny3","",600,800);
  canny3->cd();


  TH1D* PosX=new TH1D("PosX","",18,-5,5);
  PosX->SetLineColor(kRed);
  TH1D* PosY=new TH1D("PosY","",18,-5,5);
  PosY->SetLineColor(kGreen+1); 
  //Clusters with only 2 Fibres:
  TH1D* PosX_2=new TH1D("PosX_2","",18,-5,5);
  PosX_2->SetLineColor(kMagenta+2);
  TH1D* PosY_2=new TH1D("PosY_2","",18,-5,5);
  PosY_2->SetLineColor(kTeal+3);


  tree->Draw("pos_hodoClustX>>PosX","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1)");
  tree->Draw("pos_hodoClustY>>PosY","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1)","same");  

  tree->Draw("pos_hodoClustX>>PosX_2","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1 && nFibres_hodoClustX< 3. && nFibres_hodoClustY< 3. )","same");
  tree->Draw("pos_hodoClustY>>PosY_2","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1  && nFibres_hodoClustX< 3. && nFibres_hodoClustY< 3.)","same");

  PosX->GetXaxis()->SetTitle("Position");

 PosX->Scale(1./PosX->Integral());
 PosY->Scale(1./PosY->Integral());
  PosX_2->Scale(1./PosX_2->Integral());
  PosY_2->Scale(1./PosY_2->Integral());
  
  PosX->Rebin();
  PosY->Rebin();
  PosX_2->Rebin();
  PosY_2->Rebin();
  
  TLegend* legend3 = new TLegend( 0.7, 0.9, 0.9, 0.7 );
  legend3->SetFillColor(0);
  legend3->SetFillStyle(0);
  legend3->SetLineColor(0);
  legend3->SetLineWidth(0);
  legend3->SetTextSize(0.035);
  legend3->AddEntry("PosX", "Cl_PosX" );
  legend3->AddEntry("PosY", "Cl_PosY" );
  legend3->AddEntry("PosX_2", "Cl_PosX" );
  legend3->AddEntry((TObject*)0, "max 2Fibres", "");
  legend3->AddEntry("PosY_2", "Cl_PosY" );
  legend3->AddEntry((TObject*)0, "max 2Fibres", "");
  legend3->Draw("same");


  ////////////////Saving the plots//////////

  canny3->SaveAs(Form("%s/Cluster_Position_%s_%s.pdf", outputdir.c_str(), runName.c_str(), tag.c_str() ));

  canny2->SaveAs(Form("%s/nFibres_%s_%s.pdf", outputdir.c_str(), runName.c_str(), tag.c_str() ));

  canny->SaveAs(Form("%s/Intercalibration_%s_%s.pdf", outputdir.c_str() , runName.c_str(), tag.c_str() ));
		canny->Clear();

		delete canny;
		delete canny2;
		delete canny3;
  
    return 0;
}





TStyle* setStyle() {

  // set the TStyle
  TStyle* style_ = new TStyle("DrawBaseStyle", "");
  style_->SetCanvasColor(0);
  style_->SetPadColor(0);
  style_->SetFrameFillColor(0);
  style_->SetStatColor(0);
  style_->SetOptStat(0);
  style_->SetTitleFillColor(0);
  style_->SetCanvasBorderMode(0);
  style_->SetPadBorderMode(0);
  style_->SetFrameBorderMode(0);
  style_->SetPadBottomMargin(0.12);
  style_->SetPadLeftMargin(0.12);
  style_->cd();

  // For the canvas:
  style_->SetCanvasBorderMode(0);
  style_->SetCanvasColor(kWhite);
  style_->SetCanvasDefH(600); //Height of canvas
  style_->SetCanvasDefW(600); //Width of canvas
  style_->SetCanvasDefX(0); //POsition on screen
  style_->SetCanvasDefY(0);

  // For the Pad:
  style_->SetPadBorderMode(0);
  style_->SetPadColor(kWhite);
  style_->SetPadGridX(false);
  style_->SetPadGridY(false);
  style_->SetGridColor(0);
  style_->SetGridStyle(3);
  style_->SetGridWidth(1);

  // For the frame:
  style_->SetFrameBorderMode(0);
  style_->SetFrameBorderSize(1);
  style_->SetFrameFillColor(0);
  style_->SetFrameFillStyle(0);
  style_->SetFrameLineColor(1);
  style_->SetFrameLineStyle(1);
  style_->SetFrameLineWidth(1);


  // Margins:
  style_->SetPadTopMargin(0.05);
  style_->SetPadBottomMargin(0.15);//0.13);
  style_->SetPadLeftMargin(0.15);//0.16);
  style_->SetPadRightMargin(0.05);//0.02);

  // For the Global title:

  style_->SetOptTitle(0);
  style_->SetTitleFont(42);
  style_->SetTitleColor(1);
  style_->SetTitleTextColor(1);
  style_->SetTitleFillColor(10);
  style_->SetTitleFontSize(0.05);

  // For the axis titles:

  style_->SetTitleColor(1, "XYZ");
  style_->SetTitleFont(42, "XYZ");
  style_->SetTitleSize(0.05, "XYZ");
  style_->SetTitleXOffset(1.15);//0.9);
  style_->SetTitleYOffset(1.4); // => 1.15 if exponents

  // For the axis labels:

  style_->SetLabelColor(1, "XYZ");
  style_->SetLabelFont(42, "XYZ");
  style_->SetLabelOffset(0.007, "XYZ");
  style_->SetLabelSize(0.045, "XYZ");

  // For the axis:

  style_->SetAxisColor(1, "XYZ");
  style_->SetStripDecimals(kTRUE);
  style_->SetTickLength(0.03, "XYZ");
  style_->SetNdivisions(510, "XYZ");
  style_->SetPadTickX(1); // To get tick marks on the opposite side of the frame
  style_->SetPadTickY(1);

  style_->cd();

  return style_;

}
