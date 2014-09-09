#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>

#include "TH1.h"
#include "TString.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h" 
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TLine.h"
#include "TStyle.h"

#include "interface/DrawTools.h"

#include "fastDQM_CeF3_BTF.h"
#include "interface/HodoCluster.h"
#include "interface/RunHelper.h"
#include "interface/CalibrationUtility.h"
#include "interface/EnergyCalibration.h"

#include "TApplication.h"

//works for the 3 high stats runs 92,246,259 (otherwise change lable in leg)
//works for the 5 central runs at different energies 259,286,293,308,314

int main(int argc, char* argv[]){
TApplication* a = new TApplication("a", 0, 0);
 TStyle* style=  DrawTools::setStyle();
 style->cd();

  std::string outputdir = "CalibComparisonPlots/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  std::string runName0 = "BTF_259_20140502-012847_beam";
  std::string runName1 = "BTF_259_20140502-012847_beam";
  std::string runName2 = "BTF_259_20140502-012847_beam";
  std::string runName3 = "BTF_259_20140502-012847_beam";
  std::string runName4 = "BTF_259_20140502-012847_beam";


  //std::string tag = "default";
  std::string tag = "V03";

  std::string inputDir = "./CeF3Calibration";


  if( argc == 5 ) {
    std::string runName_str0(argv[1]);
    runName0 = runName_str0;
    std::string runName_str1(argv[2]);
    runName1 = runName_str1;
    std::string runName_str2(argv[3]);
    runName2 = runName_str2;
    std::string tag_str(argv[4]);
    tag = tag_str;
  } else if(argc == 7){
       std::string runName_str0(argv[1]);
    runName0 = runName_str0;
    std::string runName_str1(argv[2]);
    runName1 = runName_str1;
    std::string runName_str2(argv[3]);
    runName2 = runName_str2;
    std::string runName_str3(argv[4]);
    runName3 = runName_str3;
    std::string runName_str4(argv[5]);
    runName4 = runName_str4;
    std::string tag_str(argv[6]);
    tag = tag_str;
  } else{
    std::cout<<"Usage:"<<std::endl;
    std::cout<<"./calibrateCef3 BTF_XX BTF_XX BTF_XX tag "<<std::endl;
    exit(12345);
  }




TCanvas* canny = new TCanvas("canny", "",200,200);
 canny->cd();

  ifstream in;

  float x;
  int nlines;
  TNtuple *ntuple = new TNtuple("ntuple","data from .txt file", "x");

  double corr0[4];
  double corr_uncert0[4];
  double corr1[4];
  double corr_uncert1[4];
  double corr2[4];
  double corr_uncert2[4];
  double corr3[4];
  double corr_uncert3[4];
  double corr4[4];
  double corr_uncert4[4];

 
    int const ntot = 10;
    TString openname[ntot] = {
      Form("%s/constants_%s_%s.txt",inputDir.c_str(),runName0.c_str(), tag.c_str() ), Form("%s/constants_uncert_%s_%s.txt", inputDir.c_str(),runName0.c_str(), tag.c_str()),   Form("%s/constants_%s_%s.txt", inputDir.c_str(),runName1.c_str(), tag.c_str()), Form("%s/constants_uncert_%s_%s.txt", inputDir.c_str(),runName1.c_str(), tag.c_str()),   Form("%s/constants_%s_%s.txt",inputDir.c_str(), runName2.c_str(), tag.c_str() ),Form("%s/constants_uncert_%s_%s.txt", inputDir.c_str(), runName2.c_str(), tag.c_str() ),  Form("%s/constants_%s_%s.txt", inputDir.c_str(), runName3.c_str(), tag.c_str() ), Form("%s/constants_uncert_%s_%s.txt", inputDir.c_str(), runName3.c_str(), tag.c_str() ),     Form("%s/constants_%s_%s.txt", inputDir.c_str(), runName4.c_str(), tag.c_str() ),Form("%s/constants_uncert_%s_%s.txt", inputDir.c_str(), runName4.c_str(), tag.c_str() ) };


  for(int i=0; i<ntot+1; ++i){

    in.open(openname[i-1]);


    while (nlines<4) {
      in>>x;
      if(!in.good()) break;
      ntuple->Fill(x);
      nlines++;
    }

in.close();
 nlines=0;
}

    ntuple->SetBranchAddress("x",&x);

    for(int j=0; j<4; ++j){
      ntuple->GetEntry(j);
      corr0[j]= x;
      ntuple->GetEntry(j+4);
      corr_uncert0[j]=x;

      ntuple->GetEntry(j+8);
      corr1[j]= x;
      ntuple->GetEntry(j+12);
      corr_uncert1[j]=x;

      ntuple->GetEntry(j+16);
      corr2[j]= x;
      ntuple->GetEntry(j+20);
      corr_uncert2[j]=x;

      if(argc==7){

      ntuple->GetEntry(j+24);
      corr3[j]= x;
      ntuple->GetEntry(j+28);
      corr_uncert3[j]=x;

      ntuple->GetEntry(j+32);
      corr4[j]= x;
      ntuple->GetEntry(j+36);
      corr_uncert4[j]=x;
}
    }



  double ch0[4] = {-0.1,0.9,1.9,2.9};
  double ch1[4]={-0.05,0.95,1.95, 2.95};
  double ch2[4] = {0,1,2,3};
  double ch3[4] = {0.05,1.05,2.05,3.05};
  double ch4[4] ={0.1,1.1,2.1,3.1};
  double cherr[4]={0,0,0,0};


  TGraphErrors* graf0 = new TGraphErrors(4,ch0,corr0,cherr,corr_uncert0);
  graf0->SetMarkerStyle(9);
  graf0->SetMarkerSize(0.5);
  graf0->SetMarkerColor(kGreen+1);
  TGraphErrors* graf1 = new TGraphErrors(4,ch1,corr1,cherr,corr_uncert1);
  graf1->SetMarkerStyle(9);
  graf1->SetMarkerSize(0.5);
  graf1->SetMarkerColor(kBlue);
  TGraphErrors* graf2 = new TGraphErrors(4,ch2,corr2,cherr,corr_uncert2);
  graf2->SetMarkerStyle(9);
  graf2->SetMarkerSize(0.5);
  graf2->SetMarkerColor(kMagenta);



TMultiGraph *multi=new TMultiGraph();
 multi->Add(graf0); multi->Add(graf1); multi->Add(graf2);


 TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);

 if(argc==5){
  leg->AddEntry(graf0,"Run 92","P");
  leg->AddEntry(graf1,"Run 246","P");
  leg->AddEntry(graf2,"Run 259","P");


 multi->SetTitle(";Channel Nr.;Correction Factor");
  multi->Draw("AP");
  multi->GetYaxis()->SetRangeUser(0.975,1.025);
  multi->GetXaxis()->SetNdivisions(4);




}


  if(argc==7){ 
  graf1->SetMarkerColor(kCyan+1);
  graf2->SetMarkerColor(kBlue);
 
  TGraphErrors* graf3 = new TGraphErrors(4,ch3,corr3,cherr,corr_uncert3);
  graf3->SetMarkerStyle(9);
  graf3->SetMarkerSize(0.5);
  graf3->SetMarkerColor(kViolet-1);
  TGraphErrors* graf4 = new TGraphErrors(4,ch4,corr4,cherr,corr_uncert4);
  graf4->SetMarkerStyle(9);
  graf4->SetMarkerSize(0.5);
  graf4->SetMarkerColor(kPink+7);


multi->Add(graf3); 
multi->Add(graf4);
multi->SetTitle(";Channel Nr.;Correction Factor");


 multi->Draw("AP");

  leg->AddEntry(graf0,"98.3 MeV","P");
  leg->AddEntry(graf1,"147.4 MeV","P");
  leg->AddEntry(graf2,"196.5 MeV","P");
  leg->AddEntry(graf3,"294.8 MeV","P");
  leg->AddEntry(graf4,"491.4 MeV","P");

  //multi->GetYaxis()->SetRangeUser(0.9,1.1);



  multi->GetXaxis()->SetNdivisions(4);


  }



 multi->SetTitle(";Channel Nr.;Correction Factor");




multi->Draw("AP");
canny->Update();



  leg->SetFillColor(0);
  leg->Draw("same");

  TLine* lin = new TLine(-0.25,1.,3.25,1.);
  lin->SetLineColor(kRed);
  lin->Draw();



  if(argc==5){canny->SaveAs( Form( "%s/comparison_%s.pdf", outputdir.c_str(), tag.c_str() ));
  }  else if(argc==7){  canny->SaveAs( Form( "%s/comparison_energy_%s.pdf", outputdir.c_str(), tag.c_str() ));}

 delete canny;

return 0;
}
