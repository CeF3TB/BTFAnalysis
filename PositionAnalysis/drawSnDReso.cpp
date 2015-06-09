#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGaxis.h"

#include "interface/DrawTools.h"
#include "interface/FitTools.h"

#include "TApplication.h"

//This draws the lovely comparision of Data (1x1), MC (1x1) and MC (5x5) resolution plot



struct ResoStruct {

  float resp;
  float resp_error;

  float reso;
  float reso_error;

};

ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, const std::string& name, float energyS , bool withHodo);
ResoStruct getRespAndReso( TF1* f1, float energyS );
float getRatioError( float num, float denom, float numErr, float denomErr );
ResoStruct addPhotoStatistics( ResoStruct rs, bool withHodo );

void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name, int niter, float nSigma ) ;
TF1* fitSingleElectronPeak( const std::string& outputdir, const std::string& name, TTree* tree, int niter, float nSigma , bool withHodo);

TH1D* projectHodo ( TTree* tree, const std::string& name, bool withHodo);
TH1D* projectHodoMC ( TTree* tree, const std::string& name, bool withHodo);



int main( int argc, char* argv[] ) {

  TApplication* a = new TApplication("a",0,0);

  bool withHodo=0;
  if( argc>2 ) {
    withHodo=1;
  }

  DrawTools::setStyle();

  TGaxis::SetMaxDigits(3);
 
  std::string outputdir = "SimNDataResolution/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  std::vector<std::string> runs; 
  std::vector<float> beamEnergy;
 
  std::vector<std::string> simulation; 
  std::vector<float> beamEnergySimulation; 

  std::vector<std::string> simulationIdeal3; 
  std::vector<std::string> simulationIdeal5; 
  std::vector<std::string> simulationIdeal1m; 

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

  simulation.push_back("Simulation98");
  beamEnergySimulation.push_back(98.);
  simulation.push_back("Simulation147");
  beamEnergySimulation.push_back(147.);
  simulation.push_back("Simulation196");
  beamEnergySimulation.push_back(196.);
  simulation.push_back("Simulation295");
  beamEnergySimulation.push_back(295.);
  simulation.push_back("Simulation491");
  beamEnergySimulation.push_back(491.);
 
  simulationIdeal3.push_back("SimulationIdeal_98");
  simulationIdeal3.push_back("SimulationIdeal_147");
  simulationIdeal3.push_back("SimulationIdeal_196");
  simulationIdeal3.push_back("SimulationIdeal_295");
  simulationIdeal3.push_back("SimulationIdeal_491");

  simulationIdeal5.push_back("SimulationIdeal_98");
  simulationIdeal5.push_back("SimulationIdeal_147");
  simulationIdeal5.push_back("SimulationIdeal_196");
  simulationIdeal5.push_back("SimulationIdeal_295");
  simulationIdeal5.push_back("SimulationIdeal_491");

  simulationIdeal1m.push_back("SimulationIdeal_98");
  simulationIdeal1m.push_back("SimulationIdeal_147");
  simulationIdeal1m.push_back("SimulationIdeal_196");
  simulationIdeal1m.push_back("SimulationIdeal_295");
  simulationIdeal1m.push_back("SimulationIdeal_491");


 TGraphErrors* gr_reso_data = new TGraphErrors(0);

 TGraphErrors* gr_reso_simul = new TGraphErrors(0);
 TGraphErrors* gr_reso_simul3 = new TGraphErrors(0);
 TGraphErrors* gr_reso_simul5 = new TGraphErrors(0);
 TGraphErrors* gr_reso_simul1m = new TGraphErrors(0);

 TGraphErrors* gr_effsdata = new TGraphErrors(0);
 TGraphErrors* gr_effssimul = new TGraphErrors(0);
 TGraphErrors* gr_effssimulEnergy = new TGraphErrors(0);
 TGraphErrors* gr_effssimulPerfBeam = new TGraphErrors(0);
 TGraphErrors* gr_effssimulIdealBeam = new TGraphErrors(0);

 TGraphErrors* gr_effsHodo10 = new TGraphErrors(0);
 TGraphErrors* gr_effsBeam45 = new TGraphErrors(0);
 
 TGraphErrors* gr_effsXOffset = new TGraphErrors(0);
  
 // gStyle->SetOptFit(1);


 TFile* fileMean = TFile::Open(Form("analysisTrees_V03/Reco_%s.root",  runs[4].c_str()));
 TTree* treeMean = (TTree*)fileMean->Get("recoTree");
 float xMax = 550.;
 
 //for CeF3 data:
 TF1* energyfuncC = fitSingleElectronPeak( outputdir, runs[4], treeMean,5,1.4, withHodo );
 float adcEnergyC = energyfuncC->GetParameter(1);
 
 
 //for the simulated Cef3:
 TFile* energyfileS = TFile::Open(Form("OriginalSimulationData/newReal/Reco_%s.root", simulation[4].c_str()));
 TTree* energytreeS = (TTree*)energyfileS->Get("recoTree");
 ResoStruct energyrsS = getResponseResolutionMC( outputdir, energytreeS, simulation[4], beamEnergySimulation[4], withHodo );
 float energy491S = energyrsS.resp;
 
 
 for( unsigned i=0; i<runs.size(); ++i ) {

    ////////////THE DATA//////////////
    TFile* file = TFile::Open(Form("analysisTrees_V03/Reco_%s.root", runs[i].c_str()));
    TTree* tree = (TTree*)file->Get("recoTree");

    TF1* thisFunc = fitSingleElectronPeak( outputdir, runs[i], tree,5, 1.4, withHodo );

    float energy = beamEnergy[i];
    float energyErr = 5.;
 
    float mean = thisFunc->GetParameter(1);
    float meanErr = thisFunc->GetParError(1);

    float rms = thisFunc->GetParameter(2);
    float rmsErr = thisFunc->GetParError(2);

    float reso = 100.* rms/mean ;
    float resoErr = 100* getRatioError( rms, mean, rmsErr, meanErr);

    gr_reso_data->SetPoint(i, energy, reso );
    gr_reso_data->SetPointError(i, energyErr, resoErr );


    //////////////Efficiencies////////////////
    TH1F* histo = new TH1F("histo","", 500, 0., 500. );
    TH1F* histo2 = new TH1F("histo2","", 500, 0., 500. );

    tree->Project( "histo" ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "isSingleEle_scintFront" );
    tree->Project("histo2","cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(nHodoClustersX==1 && nHodoClustersY==1&& isSingleEle_scintFront )" );

    float scints = histo->GetEntries();
    float final= histo2->GetEntries();
    
    float eff = final/scints;
    
    gr_effsdata->SetPoint(i, energy, eff*100);
  

    /////THE REAL (1x1) SIMULATION ///////////////////
    TFile* fileS = TFile::Open(Form("OriginalSimulationData/newReal/Reco_%s.root",  simulation[i].c_str())); 
    //TFile* fileS = TFile::Open(Form("OriginalSimulationData/RealWHodo/Reco_%s.root",  simulation[i].c_str()));
    
    TTree* treeS = (TTree*)fileS->Get("recoTree");
    
    ResoStruct rs = getResponseResolutionMC( outputdir, treeS, simulation[i], beamEnergySimulation[i] ,withHodo );
    ResoStruct rs_ps = addPhotoStatistics( rs,withHodo  );
    
    gr_reso_simul->SetPoint( i, beamEnergySimulation[i], rs_ps.reso );
    gr_reso_simul->SetPointError( i,0,  rs_ps.reso_error );


    //For its Efficiencies////////////////
    TH1F* histos = new TH1F( "histos", "", 500, 0., 500. );
    TH1F* histo2s = new TH1F( "histo2s", "", 500, 0., 500. );
  
    treeS->Project( "histos" ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "isSingleEle_scintFront" );

    treeS->Project( "histo2s","cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(nHodoClustersX==1 && nHodoClustersY==1 && isSingleEle_scintFront) " );

    float scintss = histos->GetEntries();
    float finals= histo2s->GetEntries();
    
    gr_effssimul->SetPoint(i,beamEnergySimulation[i], finals/scintss *100);
 

    /////THE REAL ENERGY SMEARED (1x1) SIMULATION ///////////////////
    TFile* fileSE = TFile::Open(Form("OriginalSimulationData/RealEnergy/Reco_%s.root",  simulation[i].c_str())); 
    TTree* treeSE = (TTree*)fileSE->Get("recoTree");

    TH1F* histosE = new TH1F( "histosE", "", 500, 0., 500. );
    TH1F* histo2sE = new TH1F( "histo2sE", "", 500, 0., 500. );
    
    treeSE->Project( "histosE" ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "isSingleEle_scintFront" );
    
    treeSE->Project( "histo2sE","cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(nHodoClustersX==1 && nHodoClustersY==1 && isSingleEle_scintFront) " );
    
    float scintssE = histosE->GetEntries();
    float finalsE= histo2sE->GetEntries();
  
    gr_effssimulEnergy->SetPoint(i,beamEnergySimulation[i], finalsE/scintssE *100);
 

    /////THE Perfect Beam (1x1) SIMULATION ///////////////////
    TFile* filePB = TFile::Open(Form("OriginalSimulationData/newRealWOSmear/Reco_%s.root",  simulation[i].c_str())); 
    TTree* treePB = (TTree*)filePB->Get("recoTree");
    
    
    TH1F* histosPB = new TH1F( "histosPB", "", 500, 0., 500. );
    TH1F* histo2sPB = new TH1F( "histo2sPB", "", 500, 0., 500. );
    
    treePB->Project( "histosPB" ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "isSingleEle_scintFront" );
    
    treePB->Project( "histo2sPB","cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(nHodoClustersX==1 && nHodoClustersY==1 && isSingleEle_scintFront) " );
    
    float scintssPB = histosPB->GetEntries();
    float finalsPB= histo2sPB->GetEntries();
    
    gr_effssimulPerfBeam->SetPoint(i,beamEnergySimulation[i], finalsPB/scintssPB *100);
    
    

    /////THE Ideal  Beam (1x1) SIMULATION  without air///////////////////
    TFile* fileIB = TFile::Open(Form("OriginalSimulationData/newRealWOAirSmear/Reco_%s.root",  simulation[i].c_str())); 
    
    TTree* treeIB = (TTree*)fileIB->Get("recoTree");
        
    TH1F* histosIB = new TH1F( "histosIB", "", 500, 0., 500. );
    TH1F* histo2sIB = new TH1F( "histo2sIB", "", 500, 0., 500. );
    
    treeIB->Project( "histosIB" ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "isSingleEle_scintFront" );

    treeIB->Project( "histo2sIB","cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(nHodoClustersX==1 && nHodoClustersY==1 && isSingleEle_scintFront) " );

    float scintssIB = histosIB->GetEntries();
    float finalsIB= histo2sIB->GetEntries();
    
    gr_effssimulIdealBeam->SetPoint(i,beamEnergySimulation[i], finalsIB/scintssIB *100);
    
    
    
    /*  ////THE (3x3) SIMULATION/////////////////
	TFile* file3 = TFile::Open(Form("OriginalSimulationData/big3/Reco_%s.root",  simulationIdeal3[i].c_str())); TTree* tree3 = (TTree*)file3->Get("recoTree");
	ResoStruct rs3 = getResponseResolutionMC( outputdir, tree3, simulationIdeal3[i], beamEnergySimulation[i] ); ResoStruct rs_ps3 = addPhotoStatistics( rs3 );
	gr_reso_simul3->SetPoint( i, beamEnergySimulation[i], rs_ps3.reso );
	gr_reso_simul3->SetPointError( i,0,  rs_ps3.reso_error );
    */
    
 
    ////THE (5x5) SIMULATION/////////////////
    TFile* file5 = TFile::Open(Form("OriginalSimulationData/newMat5/Reco_%s.root",  simulationIdeal5[i].c_str()));
    //   TFile* file5 = TFile::Open(Form("OriginalSimulationData/big5/Reco_%s.root",  simulationIdeal5[i].c_str()));
    TTree* tree5 = (TTree*)file5->Get("recoTree");
    
    ResoStruct rs5 = getResponseResolutionMC( outputdir, tree5, simulationIdeal5[i], beamEnergySimulation[i] ,withHodo );
    ResoStruct rs_ps5 = addPhotoStatistics( rs5 ,withHodo );
    
    gr_reso_simul5->SetPoint( i, beamEnergySimulation[i], rs_ps5.reso );
    gr_reso_simul5->SetPointError( i,0,  rs_ps5.reso_error );
    
    

    
    /////THE HODO 10x10mm (1x1) SIMULATION ///////////////////
    TFile* fileH= TFile::Open(Form("OriginalSimulationData/newReal10Hodo/Reco_%s.root",  simulation[i].c_str())); 
    TTree* treeH = (TTree*)fileH->Get("recoTree");
    
    TH1F* histoH = new TH1F( "histoH", "", 500, 0., 500. );
    TH1F* histo2H = new TH1F( "histo2H", "", 500, 0., 500. );
    
    treeH->Project( "histoH" ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "isSingleEle_scintFront" );
    
    treeH->Project( "histo2H","cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(nHodoClustersX==1 && nHodoClustersY==1 && isSingleEle_scintFront) " );
    
    float scintsH = histoH->GetEntries();
    float finalsH= histo2H->GetEntries();
    
    gr_effsHodo10->SetPoint(i,beamEnergySimulation[i], finalsH/scintsH *100);
    
    



    /////Beam smeared 4x5 mm (1x1) SIMULATION ///////////////////
    TFile* fileBS= TFile::Open(Form("OriginalSimulationData/newRealBeam45/Reco_%s.root",  simulation[i].c_str())); 
    TTree* treeBS = (TTree*)fileBS->Get("recoTree");

    TH1F* histoBS = new TH1F( "histoBS", "", 500, 0., 500. );
    TH1F* histo2BS = new TH1F( "histo2BS", "", 500, 0., 500. );
    
    treeBS->Project( "histoBS" ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "isSingleEle_scintFront" );
    
    treeBS->Project( "histo2BS","cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(nHodoClustersX==1 && nHodoClustersY==1 && isSingleEle_scintFront) " );
    
    float scintsBS = histoBS->GetEntries();
    float finalsBS = histo2BS->GetEntries();
    
    gr_effsBeam45->SetPoint(i,beamEnergySimulation[i], finalsBS/scintsBS *100);
    
    




    /////Beam Offset in X  (1x1) SIMULATION ///////////////////
    TFile* fileX= TFile::Open(Form("OriginalSimulationData/XOffset/Reco_%s.root",  simulation[i].c_str())); 
    TTree* treeX = (TTree*)fileX->Get("recoTree");

    TH1F* histoX = new TH1F( "histoX", "", 500, 0., 500. );
    TH1F* histo2X = new TH1F( "histo2X", "", 500, 0., 500. );
    
    treeX->Project( "histoX" ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "isSingleEle_scintFront" );
    
    treeX->Project( "histo2X","cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(nHodoClustersX==1 && nHodoClustersY==1 && isSingleEle_scintFront) " );
    
    float scintsX = histoX->GetEntries();
    float finalsX = histo2X->GetEntries();
    
    gr_effsXOffset->SetPoint(i,beamEnergySimulation[i], finalsX/scintsX *100);
    

  }


 TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
 c1->cd();
 

 ///////////////////////RESOLUTION///////////////////////////////
 TH2D* h2_axes2 = new TH2D( "axes", "", 10, 0., xMax, 10, 0., 70. );
 h2_axes2->SetXTitle("Electron Beam Energy [MeV]");
 h2_axes2->SetYTitle("Energy Resolution [%]");
 h2_axes2->Draw("");
 
 // Data (1x1)
 gr_reso_data->SetMarkerStyle(20);
 gr_reso_data->SetMarkerSize(1.6);
 gr_reso_data->SetMarkerColor(kBlue);
 gr_reso_data->Draw("p same");
 
 //TF1 *fun= new TF1("fun", "sqrt([0]*[0]/(x))",50, xMax);
 TF1 *fun= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]/(x*x))",50, xMax);
 fun->SetParameter(1,50.);
 gr_reso_data->Fit(fun,"RN");
 fun->SetLineWidth(1.);
 fun->SetLineColor(46);
 // fun->Draw("L same");
 
 // MC (1x1)
 gr_reso_simul->SetMarkerStyle(24);
 gr_reso_simul->SetMarkerSize(1.5);
 gr_reso_simul->SetMarkerColor(kBlue);
 gr_reso_simul->Draw("p same");
 
 TF1 *fun1= new TF1("fun1",  "sqrt([0]*[0]/x+[1]*[1]/(x*x))",50, xMax);
 fun1->SetParameter(1,50.);
 gr_reso_simul->Fit(fun1,"RN");
 fun1->SetLineWidth(1.);
 fun1->SetLineColor(38);
 // fun1->Draw("L same");

 //MC (5x5)
 gr_reso_simul5->SetMarkerStyle(21);
 gr_reso_simul5->SetMarkerSize(1.2);
 gr_reso_simul5->SetMarkerColor(kBlack);
 gr_reso_simul5->Draw("p same");
 
 TF1 *fun5= new TF1("fun5", "sqrt([0]*[0]/(x))",23, xMax);
 gr_reso_simul5->Fit(fun5,"RN");
 fun5->SetLineWidth(1.2);
 fun5->SetLineColor(kBlack);
 fun5->Draw("L same");
 
 
 TLegend* leg4 = new TLegend(0.5, 0.55, 0.78, 0.9);
 leg4->SetTextSize(0.038);
 leg4->AddEntry(gr_reso_data,"Data (1x1)","p");
 // leg4->AddEntry(fun,Form("S = %.2f\n #pm %.2f\n %s",1./sqrt(1000)*(fun->GetParameter(0)), 1./sqrt(1000)*(fun->GetParError(0)),"%" ),"L");
 //  leg4->AddEntry( (TObject*)0,Form("N =   %.2f\n #pm %.2f\n %s",1./(1000)*(fun->GetParameter(1)), 1./(1000)*(fun->GetParError(1)),"%" ),"");
 // leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun->GetChisquare(), fun->GetNDF() ), "");
 leg4->AddEntry(gr_reso_simul,"MC (1x1)","p");
 //leg4->AddEntry(fun1,Form("S = %.2f\n #pm %.2f\n %s",1./sqrt(1000)*(fun1->GetParameter(0)), 1./sqrt(1000)*(fun1->GetParError(0)),"%" ),"L");
 // leg4->AddEntry((TObject*)0,Form("N =   %.2f\n #pm %.2f\n %s",1./(1000)*(fun1->GetParameter(1)), 1./(1000)*(fun1->GetParError(1)),"%" ),"");
 // leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun1->GetChisquare(), fun1->GetNDF() ), "");
 /*
   leg4->AddEntry(gr_reso_simul3,"MC (3x3)","p");
   //leg4->AddEntry(fun1,Form("S = %.2f\n #pm %.2f\n ",1./sqrt(1000)*(fun3->GetParameter(0)), 1./sqrt(1000)*(fun3->GetParError(0)) ),"L");
   leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun3->GetChisquare(), fun3->GetNDF() ), "");
 */
 leg4->AddEntry(gr_reso_simul5,"MC (5x5)","p");
 leg4->AddEntry(fun5,Form("S = %.2f\n %s / #sqrt{E [GeV]}",1./sqrt(1000)*(fun5->GetParameter(0)), "%" ),"L");
 //  leg4->AddEntry(fun5,Form("S = %.2f\n #pm %.2f\n %s /#sqrt{E [GeV]}",1./sqrt(1000)*(fun5->GetParameter(0)), 1./sqrt(1000)*(fun5->GetParError(0)),"%" ),"L");
 // leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun5->GetChisquare(), fun5->GetNDF() ), "");
 
  //  leg4->AddEntry(fun3,"S = 10 ","L");
 
 leg4->SetFillColor(0);
 leg4->Draw("same");
 
 TPaveText* label_low = new TPaveText(0.165,0.175,0.5,0.21, "brNDC");
 label_low->SetFillColor(kWhite);
 label_low->SetTextSize(0.038);
 label_low->SetTextAlign(11); // align right
 label_low->SetTextFont(62);
 label_low->AddText( "W-CeF_{3} Single Tower");
 label_low->Draw("same");
 
 TPaveText* label_top = new TPaveText();
 label_top = DrawTools::getLabelTop("Single Electron Beam");
 label_top->Draw("same");

 
 if(withHodo==0){  c1->SaveAs( Form( "%s/resolutionWOHodo.pdf", outputdir.c_str() ) );}
 else if(withHodo==1){  c1->SaveAs( Form( "%s/resolutionWHodo.pdf", outputdir.c_str() ) );}
 
 
 ///////////EFFICIENCIES //////////////////////////

 c1->Clear();
 
 TH2D* h2_axes = new TH2D( "axes", "", 10, 0., xMax, 10, 0.,102. );
 h2_axes->SetXTitle("Electron Beam Energy [MeV]");
 h2_axes->SetYTitle("Efficiency (Hodo|scintFront) [%]");
 h2_axes->Draw("");
 
 // Data (1x1)
 gr_effsdata->SetMarkerStyle(20);
 gr_effsdata->SetMarkerSize(1.6);
 gr_effsdata->SetMarkerColor(46);
 gr_effsdata->Draw("p same");
 // MC (1x1)
 gr_effssimul->SetMarkerStyle(24);
 gr_effssimul->SetMarkerSize(1.6);
 gr_effssimul->SetMarkerColor(46);
 gr_effssimul->Draw("p same");
 // MC (1x1) Energy Smeared
 gr_effssimulEnergy->SetMarkerStyle(25);
 gr_effssimulEnergy->SetMarkerSize(1.6);
 gr_effssimulEnergy->SetMarkerColor(38);
 //  gr_effssimulEnergy->Draw("p same");
 // MC (1x1) Perfect Beam
 gr_effssimulPerfBeam->SetMarkerStyle(25);
 gr_effssimulPerfBeam->SetMarkerSize(1.6);
 gr_effssimulPerfBeam->SetMarkerColor(38);
 gr_effssimulPerfBeam->Draw("p same");
 // MC (1x1) Ideal Beam
 gr_effssimulIdealBeam->SetMarkerStyle(21);
 gr_effssimulIdealBeam->SetMarkerSize(1.6);
 gr_effssimulIdealBeam->SetMarkerColor(38);
 gr_effssimulIdealBeam->Draw("p same");
 
 // MC (1x1) 10 mm HodoScope
 gr_effsHodo10->SetMarkerStyle(22);
 gr_effsHodo10->SetMarkerSize(1.6);
 gr_effsHodo10->SetMarkerColor(38);
 // gr_effsHodo10->Draw("p same");
 // MC (1x1) 4x5mm Beam
 gr_effsBeam45->SetMarkerStyle(26);
 gr_effsBeam45->SetMarkerSize(1.6);
 gr_effsBeam45->SetMarkerColor(38);
 // gr_effsBeam45->Draw("p same");
 
 // MC (1x1) 4x5mm Beam
 gr_effsXOffset->SetMarkerStyle(22);
 gr_effsXOffset->SetMarkerSize(1.6);
 gr_effsXOffset->SetMarkerColor(38);
 gr_effsXOffset->Draw("p same");
 
 
 TLegend* leg = new TLegend(0.55, 0.18, 0.78, 0.4);
 leg->SetTextSize(0.038);
 leg->AddEntry(gr_effsdata,"Data","p");
 leg->AddEntry(gr_effssimul,"MC","p");
 // leg->AddEntry(gr_effssimulEnergy,"MC (1x1) smeared E","p");
 leg->AddEntry(gr_effssimulPerfBeam,"MC Pencil Beam","p");
 leg->AddEntry(gr_effssimulIdealBeam,"MC no Air","p");
 //leg->AddEntry(gr_effsHodo10,"MC 10mm Hodo","p");
 // leg->AddEntry(gr_effsBeam45,"MC 4x5mm Beam","p");
 leg->AddEntry(gr_effsXOffset,"MC XOffset Corr","p");
 
 leg->SetFillColor(0);
 leg->Draw("same");
 
 TPaveText* label_top2 = new TPaveText();
 label_top2 = DrawTools::getLabelTop("Single Electron Beam");
 label_top2->Draw("same");
 
 
 c1->SaveAs( Form( "%s/effs.pdf", outputdir.c_str() ) );


 return 0;
}






ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, const std::string& name, float energyS, bool withHodo  ) {

  TH1D* h1 = new TH1D( name.c_str(), "", 500, 0., 500. );

  h1 =  projectHodoMC(tree, name.c_str(), withHodo);

  TF1* f1 = new TF1( Form("f1_%s", name.c_str()), "gaus", 5., 500.);

  FitTools::doSingleFit( h1, f1, outputdir, name,5,1.5 );

  ResoStruct rs = getRespAndReso( f1, energyS );

  return rs;

}



ResoStruct getRespAndReso( TF1* f1, float energyS) {

  // float energy = energyS;
  float energyErrS = 0;

  float meanS = f1->GetParameter(1);
  float meanErrS = f1->GetParError(1);

  float rmsS = f1->GetParameter(2);
  float rmsErrS = f1->GetParError(2);

  // float resoS = 100.* rmsS/energyS; //in percent
  float resoS = 100.* rmsS/meanS; //in percent
  float resoErrS = 100.* getRatioError( rmsS, meanS, meanErrS, rmsErrS );

  ResoStruct rs;
  rs.resp = meanS;
  rs.resp_error = meanErrS;
  rs.reso = resoS;
  rs.reso_error = resoErrS;

  return rs;

}


float getRatioError( float num, float denom, float numErr, float denomErr ) {

  return sqrt( numErr*numErr/(denom*denom) + denomErr*denomErr*num*num/(denom*denom*denom*denom) );

}


ResoStruct addPhotoStatistics( ResoStruct rs, bool withHodo  ) {

  // MC response is already in MeV, 0.49 is the number of p.e. per MeV
  // RM = 90% of energy, so assume 90% of energy (0.9*491) is deposited in central channel
  //float nADC = rs.resp/173.097*3235.76; //with hodo cut
  float nADC;

  if(withHodo==1){
    nADC = rs.resp/173.097*3235.76; //with hodo cut
  }  else{
    nADC = rs.resp/169.445*3224.36;
  }

  float nPhotoElectrons = nADC/27.3;
  
  float poissonError = 100./sqrt( nPhotoElectrons ); // in percent
  
  float resoUnsmeared = rs.reso;
  
  float resoSmeared =  sqrt( rs.reso*rs.reso + poissonError*poissonError );
  
  rs.reso = resoSmeared;

  rs.reso_error = rs.reso_error * resoSmeared / resoUnsmeared; 

  return rs;

} 







TF1* fitSingleElectronPeak( const std::string& outputdir, const std::string& name, TTree* tree, int niter, float nSigma, bool withHodo  ) {


  //  TH1D* h1 = new TH1D(name.c_str(), "", 200, 0., 6000.);
  TH1D* h1 ;
  if(name == "BTF_314_20140503-024715_beam")
    h1 = new TH1D(name.c_str(), "", 200, 0., 1500.);
  else if( name == "259_20140502-012847_beam")
    h1 = new TH1D(name.c_str(), "", 200, 1000., 6000.);
  else 
    h1 = new TH1D(name.c_str(), "", 200, 0., 60000.);


  h1 =  projectHodo(tree, name.c_str(), withHodo );


  // gStyle->SetOptFit(1);

  TF1* f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 300., 5000.);
  f1->SetParameter(0, h1->Integral() );
  f1->SetParameter(1, h1->GetMean() );
  f1->SetParameter(2, h1->GetRMS() );

  f1->SetParError(1, h1->GetMeanError() );
  f1->SetParError(2, h1->GetRMSError() );

  doSingleFit( h1, f1, outputdir, name, niter, nSigma );
  //  FitTools::doSingleFit( h1, f1, outputdir, name, niter, nSigma );

  return f1;

}




TH1D* projectHodoMC ( TTree* tree, const std::string& name , bool withHodo){
  
  TH1D* h1 = new TH1D(name.c_str(), "", 200, 5., 500.);
  
  if(withHodo==1){
    tree->Project( name.c_str() ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 )");}
  else{
    tree->Project( name.c_str() ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]","isSingleEle_scintFront");}
  
  
  return h1;
}



TH1D* projectHodo ( TTree* tree, const std::string& name , bool withHodo){


  TH1D* h1 ;
  if(name == "BTF_314_20140503-024715_beam")
    h1 = new TH1D(name.c_str(), "", 100, 75., 1275.);
  else if( name == "BTF_259_20140502-012847_beam")
    h1 = new TH1D(name.c_str(), "", 100, 1500., 6600.);
  else 
    h1 = new TH1D(name.c_str(), "", 100, 100., 4000.);

  //  TH1D* h1 = new TH1D(name.c_str(), "", 200, 100., 6000.);
  
  if(withHodo==1){
    tree->Project( name.c_str() ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 )");}
  else{
    tree->Project( name.c_str() ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]","isSingleEle_scintFront");}
  

  return h1;
}













void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name, int niter, float nSigma ) {
  
  //gStyle->SetOptFit(1);

  h1->Fit( f1, "RQN" );
  f1->SetLineColor(kRed);
  f1->SetLineWidth(3);
  f1->SetParLimits(1,2.,5000);
  f1->SetParLimits(2,2.,1000.);

  for( int iter=0; iter<niter; iter++ ) {

    float mean  = f1->GetParameter(1);
    float sigma = f1->GetParameter(2);
    float fitMin = mean - nSigma*sigma;
    //float fitMin = mean - nSigma*sigma/5.;
    //    if((mean-nSigma*sigma) <0. && mean >100.) {fitMin=110.;}
    float fitMax = mean + nSigma*sigma;
    f1->SetRange( fitMin, fitMax );
    if( iter==(niter-1) ){
      h1->Fit( f1, "RQ+" );
    }
    else{
      h1->Fit( f1, "RQN" ); }
  }
  h1->SetLineWidth(2);

  int nBins = h1->GetNbinsX();
  float lowerEdge = h1->GetMinimum();
  float upperEdge = h1->GetMaximum();
  TCanvas* c1 = new TCanvas("cX", "", 600, 600);
  c1->cd();
  h1->GetXaxis()->SetTitle("ADC Channel");
  h1->GetYaxis()->SetTitle(Form("Events / (%.0f ADC Channel)", (upperEdge-lowerEdge)/nBins));
  h1->Draw();

  TPaveText* label_fit = new TPaveText(0.52,0.89-0.06*3,0.89,0.89, "brNDC");
  label_fit->SetFillColor(kWhite);
  label_fit->SetTextSize(0.038);
  label_fit->SetTextAlign(10); // align right
  label_fit->SetTextFont(62);
  label_fit->AddText("W-CeF_{3} Single Tower");

  std::string mu_str = Form("#mu = %.0f #pm %.0f ", f1->GetParameter(1), f1->GetParError(1) );
  label_fit->AddText(mu_str.c_str());
  std::string sigma_str = Form("#sigma = %.0f #pm %.0f ", f1->GetParameter(2), f1->GetParError(2) );
  label_fit->AddText(sigma_str.c_str());
  label_fit->Draw("same");

  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.09);
  gPad->SetBottomMargin(0.15);
  //  gStyle->SetTitleYOffset(0.7); // => 1.15 if exponents
  

  
  TPaveText* label_top =  new TPaveText(0.3,0.953,0.945,0.975, "brNDC");
  label_top->SetFillColor(kWhite);
  label_top->SetTextSize(0.038);
  label_top->SetTextAlign(31); // align right
  label_top->SetTextFont(62);

  if(name == "BTF_314_20140503-024715_beam")
    label_top->AddText("98 MeV Electron Beam");
  else if( name == "BTF_259_20140502-012847_beam")
    label_top->AddText("491 MeV Electron Beam");
  else 
    label_top->AddText("Electron Beam");
  

  //  label_top->AddText("Electron Beam");
  label_top->Draw("same");
  

  c1->SaveAs( Form("%s/fit_%s.png", outputdir.c_str(), name.c_str()) );
  c1->SaveAs( Form("%s/fit_%s.pdf", outputdir.c_str(), name.c_str()) );
  c1->SaveAs( Form("%s/fit_%s.eps", outputdir.c_str(), name.c_str()) );
  delete c1;

}

