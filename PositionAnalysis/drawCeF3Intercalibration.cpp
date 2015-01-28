#include <iostream>
#include <string>
#include <stdlib.h>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TEllipse.h"
#include "TString.h"
#include "TGaxis.h" 
#include "TAxis.h" 

#include "interface/DrawTools.h"

#include "TApplication.h"

//This serves the purpose to see if the intercalibration was sort of successful,
//including some cuts on the nubmer of fibres per cluster and position of said.
//Furthermore it produces a nice plot of the mulitiplicity of the cluster position.
//The last (outcommented part) is just some random asymmetry tests.

int main ( int argc, char* argv[] ) {

  TApplication* a = new TApplication("a", 0, 0);

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

  DrawTools::setStyle();


  TFile* file = TFile::Open(Form("%s/Reco_%s.root", inputDir.c_str(),runName.c_str()));
  std::cout<<"opening file:"<<file->GetName();

  TTree* tree = (TTree*)file->Get("recoTree");

  TCanvas* canny = new TCanvas("canny", "CeF3 Intercalibration", 1200,800);
 
 canny->Divide(2,1);
  canny->cd(1);

  gPad->SetLogy();

  int nBins = 75;
  float xMin = 200.;
  float xMax = 3500.;

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

  /*
  tree->Draw("cef3[0]>>hu0");
  tree->Draw("cef3[1]>>hu1","","same");
  tree->Draw("cef3[2]>>hu2","","same");
  tree->Draw("cef3[3]>>hu3","","same");
  */

  tree->Draw("cef3_corr[0]>>hu0","isSingleEle_scintFront==1","");
  tree->Draw("cef3_corr[1]>>hu1","isSingleEle_scintFront==1","same");
  tree->Draw("cef3_corr[2]>>hu2","isSingleEle_scintFront==1","same");
  tree->Draw("cef3_corr[3]>>hu3","isSingleEle_scintFront==1","same");

  TLegend* legend2 = new TLegend( 0.5, 0.9, 0.9, 0.6 );
  legend2->SetFillStyle(0);
  legend2->SetLineColor(0);
  legend2->SetLineWidth(0);
  legend2->SetTextSize(0.035);
  legend2->AddEntry((TObject*)0,"scintFront only","");
  legend2->AddEntry( hu0, Form("mean0= %1.f\n #pm %1.f\n",hu0->GetMean(),hu0->GetMeanError()) );
  legend2->AddEntry( hu1,Form("mean1= %1.f\n #pm %1.f\n",hu1->GetMean(),hu1->GetMeanError()) );
  legend2->AddEntry( hu2, Form("mean2= %1.f\n #pm %1.f\n",hu2->GetMean(),hu2->GetMeanError()) );
  legend2->AddEntry( hu3, Form("mean3= %1.f\n #pm %1.f\n",hu3->GetMean() ,hu3->GetMeanError() ));
  legend2->Draw("same");

  canny->cd(2);


  gPad->SetLogy();

  tree->Draw("cef3_corr[0]>>h0","isSingleEle_scintFront==1&& nHodoClustersX==1 && nHodoClustersY==1");
  tree->Draw("cef3_corr[1]>>h1","isSingleEle_scintFront==1&& nHodoClustersX==1 && nHodoClustersY==1","same");
  tree->Draw("cef3_corr[2]>>h2","isSingleEle_scintFront==1&& nHodoClustersX==1 && nHodoClustersY==1","same");
  tree->Draw("cef3_corr[3]>>h3","isSingleEle_scintFront==1&& nHodoClustersX==1 && nHodoClustersY==1","same");

  TLegend* legend2f = new TLegend( 0.5, 0.9, 0.9, 0.6 );
  legend2f->SetFillStyle(0);
  legend2f->SetLineColor(0);
  legend2f->SetLineWidth(0);
  legend2f->SetTextSize(0.035);
  legend2f->AddEntry((TObject*)0,"scintFront&1Cluster","");

  legend2f->AddEntry( h0, Form("mean0= %1.f\n #pm %1.f\n",h0->GetMean(),h0->GetMeanError()) );
  legend2f->AddEntry( h1,Form("mean1= %1.f\n #pm %1.f\n",h1->GetMean(),h1->GetMeanError()) );
  legend2f->AddEntry( h2, Form("mean2= %1.f\n #pm %1.f\n",h2->GetMean(),h2->GetMeanError()) );
  legend2f->AddEntry( h3, Form("mean3= %1.f\n #pm %1.f\n",h3->GetMean() ,h3->GetMeanError() ));

  legend2f->Draw("same");



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

    std::vector<std::string> runs; 
  std::vector<float> beamEnergy;
 
  runs.push_back("BTF_314_20140503-024715_beam");
  beamEnergy.push_back(98.3);
  runs.push_back("BTF_308_20140503-002534_beam");
  beamEnergy.push_back(147.4);
  runs.push_back("BTF_293_20140502-180258_beam");
  beamEnergy.push_back(196.5);
  runs.push_back("BTF_286_20140502-153528_beam");
  beamEnergy.push_back(294.8);
  runs.push_back("BTF_259_20140502-012847_beam");
  beamEnergy.push_back(491.4);
  
  TGraphErrors* gr_mean_x = new TGraphErrors(0);
  TGraphErrors* gr_mean_y = new TGraphErrors(0);

  TGraphErrors* gr_rms_x = new TGraphErrors(0);
  TGraphErrors* gr_rms_y = new TGraphErrors(0);

  TCanvas* canny3 = new TCanvas("canny3","",600,600);
  canny3->cd();
  DrawTools::setStyle();
  TGaxis::SetMaxDigits(3);


  for( unsigned i=0; i<runs.size(); ++i ) {

    TFile* file2 = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(), runs[i].c_str()));
    TTree* tree2 = (TTree*)file2->Get("recoTree");

    float energy = beamEnergy[i];
    float energyErr = 5.;


  TH2D* h2_axes2 = new TH2D( "axes", "", 18, -4.5,4.5, 10, 0., 17000 );
  h2_axes2->SetXTitle("Position [mm]");
  h2_axes2->SetYTitle("Entries");
  h2_axes2->Draw("");

  TH1D* PosX=new TH1D("PosX","",20,-5,5);
  PosX->SetLineColor(46);
  TH1D* PosY=new TH1D("PosY","",20,-5,5);
  PosY->SetLineColor(38); 
  //Clusters with only 2 Fibres:
  //  TH1D* PosX_2=new TH1D("PosX_2","",20,-5,5); //  PosX_2->SetLineColor(kMagenta+2);
  //  TH1D* PosY_2=new TH1D("PosY_2","",18,-5,5); //  PosY_2->SetLineColor(kTeal+3);

  tree2->Draw("pos_hodoClustX>>PosX","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1)","same");
  // tree2->Draw("pos_hodoClustY>>PosY","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1 && pos_hodoClustY<2.5 && -2.5<pos_hodoClustY)","same");
    tree2->Draw("pos_hodoClustY>>PosY","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1 && pos_hodoClustY<3. && -3.<pos_hodoClustY)","same");  
  // tree->Draw("pos_hodoClustX>>PosX_2","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1 && nFibres_hodoClustX< 3. && nFibres_hodoClustY< 3. )","same");   // tree->Draw("pos_hodoClustY>>PosY_2","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1  && nFibres_hodoClustX< 3. && nFibres_hodoClustY< 3.)","same");
  PosX->Rebin();   PosY->Rebin();
  PosX->SetLineWidth(2.);  PosY->SetLineWidth(2.);

  TPaveText* label_top = new TPaveText();
  label_top = DrawTools::getLabelTop(Form("%.0f MeV Electron Beam", energy));
  label_top->Draw("same");

  TLegend* legend3 = new TLegend( 0.6, 0.9, 0.9, 0.75 );
  legend3->SetFillColor(0);
  legend3->SetFillStyle(0);
  legend3->SetLineColor(0);
  legend3->SetLineWidth(0);
  legend3->SetTextSize(0.035);  //  legend3->AddEntry("PosX","X" );
  legend3->AddEntry("PosX", Form("#mu_{X} = %.2f #pm %.2f",PosX->GetMean(),PosX->GetMeanError()), "L");  // legend3->AddEntry("PosY", "Y");
  legend3->AddEntry("PosY", Form("#mu_{Y} = %.2f #pm %.2f",PosY->GetMean(),PosY->GetMeanError()), "L");
  legend3->Draw("same");

  canny3->SaveAs(Form("%s/Cluster_Position_%.0f_%s.pdf", outputdir.c_str(), energy, tag.c_str() ));

  gr_mean_x->SetPoint( i, energy, PosX->GetMean() );
  gr_mean_x->SetPointError( i, energyErr, PosX->GetMeanError() );
  gr_mean_y->SetPoint( i, energy, PosY->GetMean() );
  gr_mean_y->SetPointError( i, energyErr, PosY->GetMeanError() );

  gr_rms_x->SetPoint( i, energy, PosX->GetRMS() );
  gr_rms_x->SetPointError( i, energyErr, PosX->GetRMSError() );
  gr_rms_y->SetPoint( i, energy, PosY->GetRMS() );
  gr_rms_y->SetPointError( i, energyErr, PosY->GetRMSError() );


  std::cout << "mean x = " << PosX->GetMean() << std::endl;

	canny3->Clear();
}
  canny3->cd();




  TH2D* h2_axes = new TH2D( "axes", "", 10, 0., 550, 10, -1.,1. );
  h2_axes->SetXTitle("Electron Beam Energy [MeV]");
  h2_axes->SetYTitle("Mean Position [mm]");
  h2_axes->Draw("");

  gr_mean_x->SetMarkerStyle(20);
  gr_mean_x->SetMarkerSize(1.6);
  gr_mean_x->SetMarkerColor(46);
  gr_mean_x->Draw("p same");

  gr_mean_y->SetMarkerStyle(20);
  gr_mean_y->SetMarkerSize(1.6);
  gr_mean_y->SetMarkerColor(38);
  gr_mean_y->Draw("p same");

  TLegend* leg0 = new TLegend(0.75, 0.65, 0.85, 0.85);
  leg0->SetTextSize(0.038);
  leg0->AddEntry(gr_mean_x,"X","p");

  leg0->AddEntry(gr_mean_y,"Y","p");

  leg0->SetFillColor(0);
  leg0->Draw("same");

  TPaveText* label_top2 = new TPaveText();
  label_top2 = DrawTools::getLabelTop("Single Electron Beam");
  label_top2->Draw("same");


 canny3->SaveAs(Form("%s/Cluster_Position.pdf", outputdir.c_str() ));


 /////////// RMS PLOT //////////////////////////////////
 canny3->Clear();
  TH2D* h2_axes5 = new TH2D( "axes", "", 10, 0., 550, 10, 0.,3. );
  h2_axes5->SetXTitle("Electron Beam Energy [MeV]");
  h2_axes5->SetYTitle("RMS [mm]");
  h2_axes5->Draw("");

  gr_rms_x->SetMarkerStyle(20);
  gr_rms_x->SetMarkerSize(1.6);
  gr_rms_x->SetMarkerColor(46);
  gr_rms_x->Draw("p same");

  gr_rms_y->SetMarkerStyle(20);
  gr_rms_y->SetMarkerSize(1.6);
  gr_rms_y->SetMarkerColor(38);
  gr_rms_y->Draw("p same");

  TLegend* leg00 = new TLegend(0.75, 0.7, 0.85, 0.85);
  leg00->SetTextSize(0.038);
  leg00->AddEntry(gr_rms_x,"X ","p");
  leg00->AddEntry(gr_rms_y,"Y ","p");
  leg00->SetFillColor(0);
  leg00->Draw("same");

  label_top2->Draw("same");


  canny3->SaveAs(Form("%s/Cluster_PositionRMS.pdf", outputdir.c_str() ));


 	canny3->Clear();





  ////////////////Saving the plots//////////
	/*
 

  //  canny2->SaveAs(Form("%s/nFibres_%s_%s.pdf", outputdir.c_str(), runName.c_str(), tag.c_str() ));

  //canny->SaveAs(Form("%s/Intercalibration_%s_%s.pdf", outputdir.c_str() , runName.c_str(), tag.c_str() ));
		canny->Clear();


		canny3->Clear();
		canny->cd();

		canny->Divide(3,1);

		/*
        xMin= -1.1;
        xMax = 1.1;
	int nBinss = 30;		

canny->cd(1);
	gPad->SetLogy();
	TH1D* histo = new TH1D("histo","",nBinss,xMin,xMax);
	TH1D* histo2 = new TH1D("histo2","",nBinss,xMin,xMax);
tree->Draw("(cef3_corr[0]+cef3_corr[3]-cef3_corr[1]-cef3_corr[2])/(cef3_corr[0]+cef3_corr[3]+cef3_corr[1]+cef3_corr[2])>>histo","cef3_corr>5.","");
tree->Draw("(cef3_corr[0]-cef3_corr[3]+cef3_corr[1]-cef3_corr[2])/(cef3_corr[0]+cef3_corr[3]+cef3_corr[1]+cef3_corr[2])>>histo2","cef3_corr>5.","same");

  histo->SetLineColor(kBlue);
  histo2->SetLineColor(kMagenta);

  TLegend* legend5 = new TLegend( 0.4, 0.45, 0.6, 0.15 );
  legend5->SetFillStyle(0);
  legend5->SetLineColor(0);
  legend5->SetLineWidth(0);
  legend5->SetTextSize(0.035);
  legend5->AddEntry( histo, "L/R asymmetry" );
  legend5->AddEntry((TObject*)0, Form("meanLR = %.5f\n #pm %.5f\n",histo->GetMean(),histo->GetMeanError()), "");
  legend5->AddEntry( histo2, "U/D asymmetry");
 legend5->AddEntry((TObject*)0, Form("meanUD = %.5f\n #pm %.5f\n",histo2->GetMean(),histo2->GetMeanError()), "");
  legend5->Draw("same");


canny->cd(2);

 	gPad->SetLogy();
	TH1D* histos = new TH1D("histos","",nBinss,xMin,xMax);
	TH1D* histo2s = new TH1D("histo2s","",nBinss,xMin,xMax);
  histos->SetLineColor(kBlue);
  histo2s->SetLineColor(kMagenta);

tree->Draw("(cef3_corr[0]+cef3_corr[3]-cef3_corr[1]-cef3_corr[2])/(cef3_corr[0]+cef3_corr[3]+cef3_corr[1]+cef3_corr[2])>>histos","cef3_corr>5.&&isSingleEle_scintFront==1","");
tree->Draw("(cef3_corr[0]-cef3_corr[3]+cef3_corr[1]-cef3_corr[2])/(cef3_corr[0]+cef3_corr[3]+cef3_corr[1]+cef3_corr[2])>>histo2s","cef3_corr>5.&&isSingleEle_scintFront==1","same");
  histo->SetLineColor(kBlue);
  histo2->SetLineColor(kMagenta);

 TLegend* legend5s= new TLegend( 0.4, 0.45, 0.6, 0.15 );
  legend5s->SetFillStyle(0);
  legend5s->SetLineColor(0);
  legend5s->SetLineWidth(0);
  legend5s->SetTextSize(0.035);
  legend5s->AddEntry((TObject*)0, "scintFront", "");
  legend5s->AddEntry( histos, "L/R asymmetry" );
  legend5s->AddEntry((TObject*)0, Form("meanLR = %.5f\n #pm %.5f\n",histos->GetMean(),histos->GetMeanError()), "");
  legend5s->AddEntry( histo2s, "U/D asymmetry");
 legend5s->AddEntry((TObject*)0, Form("meanUD = %.5f\n #pm %.5f\n",histo2s->GetMean(),histo2s->GetMeanError()), "");

  legend5s->Draw("same");


canny->cd(3);
	gPad->SetLogy();
	TH1D* histoh = new TH1D("histoh","",nBinss,xMin,xMax);
	TH1D* histo2h = new TH1D("histo2h","",nBinss,xMin,xMax);
  histoh->SetLineColor(kBlue);
  histo2h->SetLineColor(kMagenta);

tree->Draw("(cef3_corr[0]+cef3_corr[3]-cef3_corr[1]-cef3_corr[2])/(cef3_corr[0]+cef3_corr[3]+cef3_corr[1]+cef3_corr[2])>>histoh","cef3_corr>5.&&isSingleEle_scintFront==1&& nHodoClustersX==1 && nHodoClustersY==1","");
tree->Draw("(cef3_corr[0]-cef3_corr[3]+cef3_corr[1]-cef3_corr[2])/(cef3_corr[0]+cef3_corr[3]+cef3_corr[1]+cef3_corr[2])>>histo2h","cef3_corr>5.&&isSingleEle_scintFront==1&& nHodoClustersX==1 && nHodoClustersY==1","same");
  histoh->SetLineColor(kBlue);
  histo2h->SetLineColor(kMagenta);

 TLegend* legend5h = new TLegend( 0.4, 0.45, 0.6, 0.15 );
  legend5h->SetFillStyle(0);
  legend5h->SetLineColor(0);
  legend5h->SetLineWidth(0);
  legend5h->SetTextSize(0.035);
  legend5h->AddEntry((TObject*)0, "scintFront", "");
  legend5h->AddEntry((TObject*)0, "+nHodoClust", "");
  legend5h->AddEntry( histos, "L/R asymmetry" );
  legend5h->AddEntry((TObject*)0, Form("meanLR = %.5f\n #pm %.5f\n",histoh->GetMean(),histoh->GetMeanError()), "");
  legend5h->AddEntry( histo2s, "U/D asymmetry");
 legend5h->AddEntry((TObject*)0, Form("meanUD = %.5f\n #pm %.5f\n",histo2h->GetMean(),histo2h->GetMeanError()), "");

  legend5h->Draw("same");

  canny->SaveAs(Form("%s/asymmetry_%s_%s.pdf", outputdir.c_str(), runName.c_str(), tag.c_str() ));





  canny3->Clear();
  canny->Clear();
  canny->Divide(2,1);
 
  canny->cd(1);

  gPad->SetLogy();

  legend2->Clear();
  legend2f->Clear();

  tree->Draw("cef3_corr[0]+cef3_corr[3]>>hu0","isSingleEle_scintFront==1","");
  tree->Draw("cef3_corr[1]+cef3_corr[2]>>hu1","isSingleEle_scintFront==1","same");


  legend2->SetFillStyle(0);
  legend2->SetLineColor(0);
  legend2->SetLineWidth(0);
  legend2->SetTextSize(0.035);
  legend2->AddEntry((TObject*)0,"scintFront only","");
  legend2->AddEntry( hu0, Form("0+3(L) m= %.0f\n #pm %.0f\n",hu0->GetMean(),hu0->GetMeanError()) );
  legend2->AddEntry( hu1, Form("1+2(R) m= %.0f\n #pm %.0f\n",hu1->GetMean(),hu1->GetMeanError()) );

  legend2->Draw("same");

  canny->cd(2);


  gPad->SetLogy();

  tree->Draw("cef3_corr[0]+cef3_corr[3]>>h0","isSingleEle_scintFront==1&& nHodoClustersX==1 && nHodoClustersY==1");
  tree->Draw("cef3_corr[1]+cef3_corr[2]>>h1","isSingleEle_scintFront==1&& nHodoClustersX==1 && nHodoClustersY==1","same");


  legend2f->SetFillStyle(0);
  legend2f->SetLineColor(0);
  legend2f->SetLineWidth(0);
  legend2f->SetTextSize(0.035);
  legend2f->AddEntry((TObject*)0,"scintFront&1Cluster","");
  legend2f->AddEntry( h0, Form("0+3(L) m= %.0f\n #pm %.0f\n",h0->GetMean(),h0->GetMeanError()) );
  legend2f->AddEntry( h1, Form("1+2(R) m= %.0f\n #pm %.0f\n",h1->GetMean(),h1->GetMeanError()) );


  legend2f->Draw("same");

  canny->SaveAs(Form("%s/asymmetryLRseparate_%s_%s.pdf", outputdir.c_str(), runName.c_str(), tag.c_str() ));






		delete canny;
		delete canny2;
		delete canny3;
		*/
    return 0;
}



/*

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
*/
