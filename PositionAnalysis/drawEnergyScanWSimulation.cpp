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

//If you want a plot of the linearty and resolution of data vs simulation
//this is what you want to run. For the resolution plot with the 5x5 simulation
//see drawSnDReso.cpp. 
//Also draws a BGO asymmetry plot, but it commented out atm, as
//BGO 4 for some reason was just a bit overeager.

//Usage: 
// ./drawEnergyScanWSimulation tag = without hodo cut
// ./drawEnergyScanWSimulation tag somestring = with hodo cut

struct ResoStruct {

  float resp;
  float resp_error;

  float reso;
  float reso_error;

};

ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, const std::string& name, float energyS, bool withHodo );
ResoStruct getRespAndReso( TF1* f1, float energyS );
float getRatioError( float num, float denom, float numErr, float denomErr );
ResoStruct addPhotoStatistics( ResoStruct rs , bool withHodo);



TF1* fitSingleElectronPeakBGO( const std::string& outputdir, const std::string& name, TTree* tree, int niter, float nSigma, std::string whattoProject, std::string cut );
void doSingleFitBGO( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name, int niter, float nSigma ) ;
TF1* fitSingleElectronPeakBGOSimul( const std::string& outputdir, const std::string& name, TTree* tree, int niter, float nSigma );
void doSingleFitBGOSimul( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name, int niter, float nSigma );

TF1* fitSingleElectronPeak( const std::string& outputdir, const std::string& name, TTree* tree, int niter, float nSigma, bool withHodo );
TH1D* projectHodo ( TTree* tree, const std::string& name, bool withHodo);
TH1D* projectHodoMC ( TTree* tree, const std::string& name, bool withHodo);

int main( int argc, char* argv[] ) {

  TApplication* a = new TApplication("a",0,0);


  std::string tag="V00";
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }

  bool withHodo=0;
  if( argc>2 ) {
    withHodo=1;
  }
 
  std::cout << "-> Using tag: " << tag << std::endl;

  std::string outputdir = "EnergyScanSimulation/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  DrawTools::setStyle();

  TGaxis::SetMaxDigits(3);
 
  std::vector<std::string> runs; 
  std::vector<float> beamEnergy;
 
  std::vector<std::string> simulation; 
  std::vector<float> beamEnergySimulation; 

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

  std::vector<std::string> simulationIdeal; 
  simulationIdeal.push_back("SimulationIdeal_98");
  simulationIdeal.push_back("SimulationIdeal_147");
  simulationIdeal.push_back("SimulationIdeal_196");
  simulationIdeal.push_back("SimulationIdeal_295");
  simulationIdeal.push_back("SimulationIdeal_491");


  TGraphErrors* gr_resp_vs_energy = new TGraphErrors(0);
  TGraphErrors* gr_reso_vs_energy = new TGraphErrors(0);

  TGraphErrors* gr_resp_vs_energy_simul = new TGraphErrors(0);
  TGraphErrors* gr_reso_vs_energy_simul = new TGraphErrors(0);

  TGraphErrors* gr_respSimul = new TGraphErrors(0);
  TGraphErrors* gr_respIdeal = new TGraphErrors(0);
  TGraphErrors* gr_flipped5 = new TGraphErrors(0);

  TGraphErrors* gr_simideal = new TGraphErrors(0);

  TGraphErrors* gr_sigmaData = new TGraphErrors(0);
  TGraphErrors* gr_sigmaSimul = new TGraphErrors(0);
  TGraphErrors* gr_sigmaSimulPS = new TGraphErrors(0);

  TGraphErrors* gr_resp_vs_energy_dev = new TGraphErrors(0);
  TGraphErrors* gr_resp_vs_energy_simul_dev = new TGraphErrors(0);


  TFile* fileMean = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(), runs[4].c_str()));
  TTree* treeMean = (TTree*)fileMean->Get("recoTree");
 
  float xMax = 550.;
 
  /*   
  float fitSigmas = 0.9;
  int iterations = 4;

  TGraphErrors* gr_bgo_asymm= new TGraphErrors(0);

  TGraphErrors* gr_resp_vs_energy_bgo = new TGraphErrors(0);
  TGraphErrors* gr_resp_vs_energy_bgo_simul = new TGraphErrors(0);
  TGraphErrors* gr_resolutionCorr = new TGraphErrors(0);
  TGraphErrors* gr_resp_vs_energy_bgo_corr = new TGraphErrors(0);
  TGraphErrors* gr_responseCorr = new TGraphErrors(0);
  TGraphErrors* gr_through0 = new TGraphErrors(0);
  

  //491MeV=peak pos of run 259 in adc for BGO!
  TF1* energyfunc = fitSingleElectronPeakBGO( outputdir, runs[4], treeMean,iterations,fitSigmas,"bgo_corr[0]+bgo_corr[1]+bgo_corr[2]+bgo_corr[3]+bgo_corr[4]+bgo_corr[5]+bgo_corr[6]+bgo_corr[7]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 )" );
  float energy491= beamEnergy[4];
  float adcEnergy = energyfunc->GetParameter(1);
  // float adcEnergyErr = energyfunc->GetParError(1);

  /*
  gr_through0->SetPoint(0, 0, 0);
  gr_through0->SetPoint(1,energy491, adcEnergy);
  gr_through0->SetPointError(1, energy491*0.01, adcEnergyErr);
  TF1* f1_line_0 = new TF1("line0", "[0]*x", 0., xMax);
  gr_through0->Fit(f1_line_0, "RN");
  float slope = f1_line_0->GetParameter(0);
  float slopeErr = f1_line_0->GetParError(0);
  */ 

  //and for CeF3:
  TF1* energyfuncC = fitSingleElectronPeak( outputdir, runs[4], treeMean,5,1.4, withHodo );

  float energy491C= beamEnergy[4];
  float adcEnergyC = energyfuncC->GetParameter(1);
  //   float adcEnergyErrC= energyfuncC->GetParError(1);


  //for the simulated Cef3:
  //TFile* energyfileS = TFile::Open(Form("OriginalSimulationData/XOffset/Reco_%s.root", simulation[4].c_str()));
  TFile* energyfileS = TFile::Open(Form("OriginalSimulationData/newReal/Reco_%s.root", simulation[4].c_str()));
  
  TTree* energytreeS = (TTree*)energyfileS->Get("recoTree");
    
  ResoStruct energyrsS = getResponseResolutionMC( outputdir, energytreeS, simulation[4], beamEnergySimulation[4], withHodo );

  float energy491S = energyrsS.resp;

  std::cout << "EnergyADC = " << adcEnergyC << " , EnergyPeakSimulation = " << energy491S << std::endl;

    /*
   ///  LOOOOOOOP //////////////////////////////
for( unsigned i=0; i<runs.size(); ++i ) {
   TFile* file = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(), runs[i].c_str()));
   TTree* tree = (TTree*)file->Get("recoTree");
   float energy = beamEnergy[i];
   
   //////BGO ASYMMMETRY /////////////
   TF1* funcBGO3 = fitSingleElectronPeakBGO(outputdir, runs[i],tree,iterations , fitSigmas, "bgo_corr[3]" ,"isSingleEle_scintFront" );
   float meanBGO3 = funcBGO3->GetParameter(1);
   float meanBGOerr3 = funcBGO3->GetParError(1);
   
   TF1* funcBGO4 = fitSingleElectronPeakBGO(outputdir, runs[i],tree,iterations , fitSigmas, "bgo_corr[4]" ,"isSingleEle_scintFront" );
   float meanBGO4 = funcBGO4->GetParameter(1);
   float meanBGOerr4 = funcBGO4->GetParError(1);
   //   float energyCorrected =  (meanBGO/slope);
   //gr_resp_vs_energy_bgo->SetPoint(i, energyCorrected, meanBGO);
   
   gr_bgo_asymm->SetPoint(i, energy, (meanBGO3-meanBGO4)/(meanBGO3+meanBGO4));
   
   TCanvas* canny = new TCanvas( "canny", "", 600, 600 );
   canny->cd();
   
   TH1D* h1 = new TH1D("h1", "", 200,-1.,1.);
   tree->Project("h1","(bgo_corr[3]-bgo_corr[4])/(bgo_corr[3]+bgo_corr[4])","bgo_corr[3]>0 && bgo_corr[4]>0 && isSingleEle_scintFront==1 && nHodoClustersX==1 && nHodoClustersY==1");
   TH1D* h2 = new TH1D("h2", "", 200,-1.,1.);
   tree->Project("h2","(bgo_corr[3]-bgo_corr[4])/(bgo_corr[3]+bgo_corr[4])","bgo_corr[3]>0 && bgo_corr[4]>0 && isSingleEle_scintFront==1 && nHodoClustersX==0 && nHodoClustersY==0 ");
   
   h1->Rebin(2);   h2->Rebin(2);
   
   if(energy<400.){  h1->Rebin(2);  h2->Rebin(2);  }
   
   float yMax= 3* h1->GetMaximum()/h1->Integral();
   
   TH2D* axes = new TH2D( "axes", "", 10, -1.1, 1.1, 10, 0., yMax);
   axes->SetXTitle("BGO Asymmetry");
   axes->SetYTitle("");
   axes->Draw("");
   
   h1->SetLineColor(46);  h2->SetLineColor(38);
   h1->SetLineWidth(2);  h2->SetLineWidth(2);
   h1->Scale(1./h1->Integral());  h2->Scale(1./h2->Integral());
   h1->Draw("same");  h2->Draw("same");
   
   TLegend* legs = new TLegend(0.30, 0.7, 0.855, 0.9);
   legs->SetTextSize(0.038);
   legs->AddEntry(h1,"ScintFront && HodoClust == 1","L");
   legs->AddEntry( (TObject*)0, Form("#mu = %f",h1->GetMean() ), "");
   legs->AddEntry(h2,"ScintFront && HodoClust == 0","L");
   legs->AddEntry( (TObject*)0, Form("#mu = %f",h2->GetMean() ), "");
   legs->SetFillColor(0);
   legs->Draw("same");
   
   canny->SaveAs( Form( "%s/BGO_asym_%d.pdf", outputdir.c_str(),int(energy) ) );
   canny->Clear();
   
   TH1D* h1u = new TH1D("h1u", "", 200,-1.,1.);
   tree->Project("h1u","(bgo_corr[1]-bgo_corr[6])/(bgo_corr[1]+bgo_corr[7])","bgo_corr[1]>0 && bgo_corr[6]>0 && isSingleEle_scintFront==1 && nHodoClustersX==1 && nHodoClustersY==1");

   TH1D* h2d = new TH1D("h2d", "", 200,-1.,1.);
   tree->Project("h2d","(bgo_corr[1]-bgo_corr[6])/(bgo_corr[1]+bgo_corr[6])","bgo_corr[1]>0 && bgo_corr[6]>0 && isSingleEle_scintFront==1 && nHodoClustersX==0 && nHodoClustersY==0 ");
   
   h1u->Rebin(2);
   h2d->Rebin(2);
      
   axes->Draw("");
   
   h1u->SetLineColor(46);
   h2d->SetLineColor(38);
   
   h1u->SetLineWidth(2);
   h2d->SetLineWidth(2);
   
   h1u->Scale(1./h1u->Integral());
   h2d->Scale(1./h2d->Integral());
   
   h1u->Draw("same");
   h2d->Draw("same");
   
   TLegend* legsf = new TLegend(0.30, 0.7, 0.855, 0.9);
   legsf->SetTextSize(0.038);
   legsf->AddEntry(h1u,"ScintFront && HodoClust == 1","L");
   legsf->AddEntry(h2d,"ScintFront && HodoClust == 0","L");
   legsf->SetFillColor(0);
   legsf->Draw("same");

   canny->SaveAs( Form( "%s/BGO_asym_UpDown_%d.pdf", outputdir.c_str(),int(energy) ) );
  

   //tree->Project( histoName.c_str(),"bgo_corr[0]+bgo_corr[1]+bgo_corr[2]+bgo_corr[3]+bgo_corr[4]+bgo_corr[5]+bgo_corr[6]+bgo_corr[7]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 )");


/*
//For the simulated BGOs
    TFile* fileS = TFile::Open(Form("OriginalSimulationData/XOffset/Reco_%s.root", tag.c_str(), simulation[i].c_str()));
    TTree* treeS = (TTree*)fileS->Get("recoTree");
    TF1* thisFuncS = fitSingleElectronPeakBGOSimul( outputdir, simulation[i], treeS,10, 0.9);
    float energyS = beamEnergySimulation[i];
    float energyErrS = 0.;   float meanS = thisFuncS->GetParameter(1);
    float meanSerr = thisFuncS->GetParError(1);   

    std::cout << "meanS bgo = " << meanS << std::endl;    
    std::cout << "mean Energy ADC = " << adcEnergy << std::endl;    
    std::cout << adcEnergy / meanS << std::endl;   
    meanS = meanS /52.1293 * adcEnergy;

    gr_resp_vs_energy_bgo_simul->SetPoint(i, energyS, meanS);   
    gr_resp_vs_energy_bgo_simul->SetPointError(i, energyErrS, meanSerr);


}

/*
  TF1* f1_line_bgo = new TF1("line", "[0] + [1]*x", 0., xMax );  gr_resp_vs_energy_bgo->Fit(f1_line_bgo, "RN");
  f1_line_bgo->SetLineWidth(1.);  f1_line_bgo->SetLineColor(46);
  float m2 = f1_line_bgo->GetParameter(1);  float m2_err = f1_line_bgo->GetParError(1);
  float q2 = f1_line_bgo->GetParameter(0);  float q2_err = f1_line_bgo->GetParError(0);
*/

    ///THE LOOOOOOOOP///
  for( unsigned i=0; i<runs.size(); ++i ) {

    TFile* file = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(), runs[i].c_str()));
    TTree* tree = (TTree*)file->Get("recoTree");

    TF1* thisFunc = fitSingleElectronPeak( outputdir, runs[i], tree,5, 1.4 , withHodo);

    float energy = beamEnergy[i];
    float energyErr = 5.;
 
    float mean = thisFunc->GetParameter(1);
    float meanErr = thisFunc->GetParError(1);

    float rms = thisFunc->GetParameter(2);
    float rmsErr = thisFunc->GetParError(2);

    float reso = 100.* rms/mean ;
    float resoErr = 100* getRatioError( rms, mean, rmsErr, meanErr);

    gr_resp_vs_energy->SetPoint( i, energy, mean );
    gr_resp_vs_energy->SetPointError( i, energyErr, meanErr );

    gr_reso_vs_energy->SetPoint( i, energy, reso );
    gr_reso_vs_energy->SetPointError( i, energyErr, resoErr );

    gr_sigmaData->SetPoint( i, energy, rms);
    gr_sigmaData->SetPointError( i, energyErr, rmsErr);

    gr_resp_vs_energy_dev->SetPoint(i, energy, mean/adcEnergyC*energy491S /energy);


    ///////////////// (1x1) Shashlik ("real setup") //////////////////////////////
    TFile* fileS = TFile::Open(Form("OriginalSimulationData/newReal/Reco_%s.root", simulation[i].c_str()));
    //  TFile* fileS = TFile::Open(Form("OriginalSimulationData/XOffset/Reco_%s.root", simulation[i].c_str()));

    TTree* treeS = (TTree*)fileS->Get("recoTree");

    ResoStruct rs = getResponseResolutionMC( outputdir, treeS, simulation[i], beamEnergySimulation[i],withHodo );
    ResoStruct rs_ps = addPhotoStatistics( rs , withHodo);

    gr_respSimul->SetPoint( i, beamEnergySimulation[i], rs.resp );
    gr_respSimul->SetPointError( i,0, rs.resp_error);

    gr_resp_vs_energy_simul->SetPoint( i, beamEnergySimulation[i], rs.resp * adcEnergyC/energy491S );
    gr_resp_vs_energy_simul->SetPointError( i,0, 0);

    gr_resp_vs_energy_simul_dev->SetPoint(i, beamEnergySimulation[i], rs.resp /beamEnergySimulation[i] );

    gr_reso_vs_energy_simul->SetPoint( i, beamEnergySimulation[i], rs_ps.reso );
    gr_reso_vs_energy_simul->SetPointError( i,0,  rs_ps.reso_error );

    gr_sigmaSimul->SetPoint(i, beamEnergySimulation[i], rs.reso*rs.resp *adcEnergyC/energy491S /100.);
  
    float nADC = rs.resp /energy491S * adcEnergyC;
    float nPhotoElectrons = nADC/27.3;
    float poissonError = 1./sqrt( nPhotoElectrons ); 
    float resoPS= sqrt( rs.reso*rs.reso*rs.resp*rs.resp *adcEnergyC/energy491S *adcEnergyC/energy491S  /100./100. + poissonError*poissonError  *nADC *nADC );

    gr_sigmaSimulPS->SetPoint(i, beamEnergySimulation[i], resoPS );

    std::cout << "Resolution = " << rs.reso << std::endl;
    std::cout << "Resolution.PE = " << rs_ps.reso << std::endl;

    //gr_simideal->SetPoint( i, beamEnergySimulation[i], rs_ps.reso );

    ///////////////// 5x5 matrix //////////////////////////////
    TFile* fileF= TFile::Open(Form("OriginalSimulationData/newMat5/Reco_%s.root",  simulationIdeal[i].c_str()));
    TTree* treeF = (TTree*)fileF->Get("recoTree");

    ResoStruct rsF = getResponseResolutionMC( outputdir, treeF, simulationIdeal[i], beamEnergySimulation[i] ,withHodo );
    ResoStruct rs_psF = addPhotoStatistics( rsF, withHodo );

    gr_simideal->SetPoint( i, beamEnergySimulation[i], rs_psF.reso );
    gr_simideal->SetPointError( i,0,  rs_psF.reso_error );

    gr_respIdeal->SetPoint( i, beamEnergySimulation[i], rsF.resp);
    gr_respIdeal->SetPointError( i,0, rsF.resp_error );

    ///////////////// 5x5 matrix REVERSED ABS ACT ORDER //////////////////////////////
    TFile* fileA= TFile::Open(Form("OriginalSimulationData/ActFirstIdeal/Reco_%s.root",  simulationIdeal[i].c_str()));
    TTree* treeA = (TTree*)fileA->Get("recoTree");

    ResoStruct rsA = getResponseResolutionMC( outputdir, treeA, simulationIdeal[i], beamEnergySimulation[i] ,withHodo );
    ResoStruct rs_psA = addPhotoStatistics( rsA, withHodo );

    gr_flipped5->SetPoint( i, beamEnergySimulation[i], rsA.resp);
    gr_flipped5->SetPointError( i,0, rsA.resp_error );

    /*
    //////BGO ///////////// TF1* funcBGO = fitSingleElectronPeakBGO(outputdir, runs[i],tree,4, 0.9); float meanBGO = funcBGO->GetParameter(1);   float meanBGOerr = funcBGO->GetParError(1);
   gr_resp_vs_energy_bgo->SetPoint(i,energy,meanBGO);gr_resp_vs_energy_bgo->SetPointError(i,energyErr,meanBGOerr);
 */
  }

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D( "axes", "", 10, 0., xMax, 10, 0., 3500. );
  h2_axes->SetXTitle("Electron Beam Energy [MeV]");
  h2_axes->SetYTitle("CeF3 Response [ADC Counts]");
  h2_axes->Draw("");

  gr_resp_vs_energy->SetMarkerStyle(20);
  gr_resp_vs_energy->SetMarkerSize(1.6);
  gr_resp_vs_energy->SetMarkerColor(46);
  gr_resp_vs_energy->Draw("p same");

  gr_resp_vs_energy_simul->SetMarkerStyle(20);
  gr_resp_vs_energy_simul->SetMarkerSize(1.6);
  gr_resp_vs_energy_simul->SetMarkerColor(38);
  gr_resp_vs_energy_simul->Draw("p same");

  TF1* f1_line = new TF1("line", "[0] + [1]*x", 0., xMax );
  gr_resp_vs_energy->Fit(f1_line, "RN");
  f1_line->SetLineColor(46);
  f1_line->SetLineWidth(1.);
  f1_line->Draw("L same");
  TF1* f1_lines = new TF1("lines", "[0] + [1]*x", 0., xMax );
  gr_resp_vs_energy_simul->Fit(f1_lines,"RN");
  f1_lines->SetLineWidth(1.);
  f1_lines->SetLineColor(38);
  f1_lines->Draw("L same");

  TLegend* leg0 = new TLegend(0.55, 0.2, 0.85, 0.45);
  leg0->SetTextSize(0.038);
  leg0->AddEntry(gr_resp_vs_energy,"Data","p");
  leg0->AddEntry(f1_line,Form("Offset = %.0f\n #pm %.0f\n",f1_line->GetParameter(0), f1_line->GetParError(0) ),"L");
  leg0->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",f1_line->GetChisquare(), f1_line->GetNDF() ), "");
  leg0->AddEntry(gr_resp_vs_energy_simul,"MC","p");
  leg0->AddEntry(f1_lines,Form("Offset = %.0f\n #pm %.0f\n",f1_lines->GetParameter(0), f1_lines->GetParError(0) ),"L");
  leg0->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.0f\n / %d",f1_lines->GetChisquare(), f1_lines->GetNDF() ), "");
  leg0->SetFillColor(0);
  leg0->Draw("same");

  TPaveText* label_top = new TPaveText();
  label_top = DrawTools::getLabelTop("Single Electron Beam");
  label_top->Draw("same");


  if(withHodo==0){  c1->SaveAs( Form( "%s/resp_vs_energyWOHodo.pdf", outputdir.c_str() ) );}
  else if(withHodo==1){  c1->SaveAs( Form( "%s/resp_vs_energyWHodo.pdf", outputdir.c_str() ) );}

  c1->Clear();

  /*
  /////////// REVERSED SETUP (= active layer first) LINEARITYYYYYYYYY////////////////
  TH2D* h2_axes5 = new TH2D( "axes", "", 10, 0., xMax, 10, 0., 250. );
  h2_axes5->SetXTitle("Electron Beam Energy [MeV]");
  h2_axes5->SetYTitle("CeF_{3} Response [MeV]");
  h2_axes5->Draw("");

  gr_respSimul->SetMarkerStyle(20);
  gr_respSimul->SetMarkerSize(1.6);
  gr_respSimul->SetMarkerColor(46);
  gr_respSimul->Draw("p same");

  gr_respIdeal->SetMarkerStyle(20);
  gr_respIdeal->SetMarkerSize(1.6);
  gr_respIdeal->SetMarkerColor(38);
  gr_respIdeal->Draw("p same");

  gr_flipped5->SetMarkerStyle(24);
  gr_flipped5->SetMarkerSize(1.6);
  gr_flipped5->SetMarkerColor(38);
  gr_flipped5->Draw("p same");

  TF1* f1_lines2 = new TF1("lines2", "[0] + [1]*x", 0., xMax );
  gr_respSimul->Fit(f1_lines2,"RN");
  f1_lines2->SetLineWidth(1.);
  f1_lines2->SetLineColor(46);
  f1_lines2->Draw("L same");

  TF1* f1_line2 = new TF1("line2", "[0] + [1]*x", 0., xMax );
  gr_respIdeal->Fit(f1_line2, "RN");
  f1_line2->SetLineColor(38);
  f1_line2->SetLineWidth(1.);
  f1_line2->Draw("L same");

  TF1* f1_linef2 = new TF1("linef2", "[0] + [1]*x", 0., xMax );
  gr_flipped5->Fit(f1_linef2, "RN");
  f1_linef2->SetLineColor(38);
  f1_linef2->SetLineWidth(1.);
  f1_linef2->Draw("L same");

  TLegend* leg9 = new TLegend(0.45, 0.9, 0.2, 0.55);
  leg9->SetTextSize(0.038);
  leg9->AddEntry(gr_respSimul,"MC (1x1)","p");
  leg9->AddEntry(f1_lines2,Form("Offset = %.1f\n #pm %.1f\n",f1_lines2->GetParameter(0), f1_lines2->GetParError(0) ),"L");
  leg9->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.0f\n / %d",f1_lines->GetChisquare(), f1_lines2->GetNDF() ), "");

  leg9->AddEntry(gr_respIdeal,"MC (5x5)","p");
  leg9->AddEntry(f1_line2,Form("Offset = %.0f\n #pm %.1f\n",f1_line2->GetParameter(0), f1_line2->GetParError(0) ),"L");
  leg9->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.0f\n / %d",f1_line2->GetChisquare(), f1_line2->GetNDF() ), "");

  leg9->AddEntry(gr_flipped5,"MC (5x5) Reversed","p");
  leg9->AddEntry(f1_linef2,Form("Offset = %.1f\n #pm %.1f\n",f1_linef2->GetParameter(0), f1_linef2->GetParError(0) ),"L");
  leg9->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.0f\n / %d",f1_linef2->GetChisquare(), f1_linef2->GetNDF() ), "");
  leg9->SetFillColor(0);
  leg9->Draw("same");
  label_top->Draw("same");

  if(withHodo==0){  c1->SaveAs( Form( "%s/reverseWOHodo.pdf", outputdir.c_str() ) );}
  else if(withHodo==1){  c1->SaveAs( Form( "%s/reverseWHodo.pdf", outputdir.c_str() ) );}

  c1->Clear();
  */


  ///////////////////////RESOLUTION///////////////////////////////
  TH2D* h2_axes2 = new TH2D( "axes", "", 10, 0., xMax, 10, 0., 60. );
  //TH2D* h2_axes2 = new TH2D( "axes", "", 10, 0., xMax, 10, 0., 1000. );
  h2_axes2->SetXTitle("Electron Beam Energy [MeV]");
  h2_axes2->SetYTitle("CeF_{3} Resolution [%]");
  //h2_axes2->SetYTitle("CeF3 Sigma [ADC]");
  h2_axes2->Draw("");

  gr_reso_vs_energy->SetMarkerStyle(20);
  gr_reso_vs_energy->SetMarkerSize(1.4);
  gr_reso_vs_energy->SetMarkerColor(46);
  gr_reso_vs_energy->Draw("p same");


  TF1* fun= new TF1("fun", "sqrt([0]*[0]/x+[1]*[1]/(x*x)) ",50, xMax);
  fun->SetParameter(1,50.);
  TF1* fun2= new TF1("fun2", "sqrt([0]*[0]/x+[1]*[1]/(x*x))",50, xMax);
  fun2->SetParameter(1,500.);

  gr_reso_vs_energy->Fit(fun,"RN");
  fun->SetLineWidth(1.);
  fun->SetLineColor(46);
  fun->Draw("L same");

  gr_reso_vs_energy_simul->SetMarkerStyle(20);
  gr_reso_vs_energy_simul->SetMarkerSize(1.4);
  gr_reso_vs_energy_simul->SetMarkerColor(38);
  gr_reso_vs_energy_simul->Draw("p same");

  gr_reso_vs_energy_simul->Fit(fun2,"RN");
  fun2->SetLineWidth(1.);
  fun2->SetLineColor(38);
  fun2->Draw("L same");

  TLegend* leg4 = new TLegend(0.55, 0.55, 0.78, 0.9);
  leg4->SetTextSize(0.038);
  leg4->AddEntry(gr_reso_vs_energy,"Data","p");
  leg4->AddEntry(fun,Form("S = %.2f\n #pm %.2f\n ",1./sqrt(1000)*(fun->GetParameter(0)), 1./sqrt(1000)*(fun->GetParError(0)) ),"L");
  leg4->AddEntry((TObject*)0,Form("N =   %.2f\n #pm %.2f\n ",1./(1000)*(fun->GetParameter(1)), 1./(1000)*(fun->GetParError(1)) ),"");
  leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun->GetChisquare(), fun->GetNDF() ), "");

  leg4->AddEntry(gr_reso_vs_energy_simul,"MC","p");
  leg4->AddEntry(fun2,Form("S = %.2f\n #pm %.2f\n ",1./sqrt(1000)*(fun2->GetParameter(0)), 1./sqrt(1000)*(fun2->GetParError(0)) ),"L");
  leg4->AddEntry((TObject*)0,Form("N =   %.2f\n #pm %.2f\n ",1./(1000)*(fun2->GetParameter(1)), 1./(1000)*(fun2->GetParError(1)) ),"");
  leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.1f\n / %d",fun2->GetChisquare(), fun2->GetNDF() ), "");
  //leg4->AddEntry(gr_simideal,"simulation PS","p");

  leg4->SetFillColor(0);
  leg4->Draw("same");

  label_top->Draw("same");

  if(withHodo==1) c1->SaveAs( Form( "%s/reso_vs_energySWHodo.pdf", outputdir.c_str() ) );
  else  c1->SaveAs( Form( "%s/reso_vs_energySWOHodo.pdf", outputdir.c_str() ) );


  /*

  c1->Clear();
  ////////////////////RESO SIMULATION TO SIMULATION COMPARISON //////////////
  h2_axes2->SetXTitle("Electron Beam Energy [MeV]");
  h2_axes2->SetYTitle("CeF3 Resolution [%]");
  h2_axes2->Draw("");

  gr_reso_vs_energy_simul->SetMarkerStyle(20);
  gr_reso_vs_energy_simul->SetMarkerSize(1.6);
  gr_reso_vs_energy_simul->SetMarkerColor(46);
  gr_reso_vs_energy_simul->Draw("p same");

  //  TF1 *fun3= new TF1("fun3", "sqrt([0]*[0]/(x))",50, xMax);
  TF1 *fun3= new TF1("fun3", "sqrt([0]*[0]/(x)+[1]*[1]/(x*x))",50, xMax);
  gr_reso_vs_energy_simul->Fit(fun3,"RN");
  fun3->SetLineWidth(1.);
  fun3->SetLineColor(46);
  fun3->Draw("L same");

  gr_simideal->SetMarkerStyle(20);
  gr_simideal->SetMarkerSize(1.6);
  gr_simideal->SetMarkerColor(38);
  gr_simideal->Draw("p same");

  //  TF1 *fun22= new TF1("fun22", "sqrt([0]*[0]/(x))",50, xMax);
  TF1 *fun22= new TF1("fun22", "sqrt([0]*[0]/(x)+[1]*[1]/(x*x))",50, xMax);
  gr_simideal->Fit(fun22,"RN");
  fun22->SetLineWidth(1.);
  fun22->SetLineColor(38);
  fun22->Draw("L same");

  TLegend* legb = new TLegend(0.45, 0.65, 0.78, 0.9);
  legb->SetTextSize(0.038);
  legb->AddEntry(gr_reso_vs_energy_simul,"MC (1x1)","p");
  legb->AddEntry(fun,Form("S = %.2f\n #pm %.2f\n ",1./sqrt(1000)*(fun3->GetParameter(0)), 1./sqrt(1000)*(fun3->GetParError(0)) ),"L");
  legb->AddEntry( (TObject*)0,Form("N = %.2f\n #pm %.2f\n ",1./sqrt(1000)*(fun3->GetParameter(1)), 1./sqrt(1000)*(fun3->GetParError(1)) ),"L");
  legb->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun3->GetChisquare(), fun3->GetNDF() ), "");
  legb->AddEntry(gr_simideal,"MC (5x5)","p");
  legb->AddEntry(fun22,Form("S = %.2f\n #pm %.2f\n ",1./sqrt(1000)*(fun22->GetParameter(0)), 1./sqrt(1000)*(fun22->GetParError(0)) ),"L");
  legb->AddEntry((TObject*)0,Form("S = %.2f\n #pm %.2f\n ",1./sqrt(1000)*(fun22->GetParameter(1)), 1./sqrt(1000)*(fun22->GetParError(1)) ),"L");
  legb->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun22->GetChisquare(), fun22->GetNDF() ), "");
  legb->SetFillColor(0);
  legb->Draw("same");

  label_top->Draw("same");

  c1->SaveAs( Form( "%s/SimToSimComparison.pdf", outputdir.c_str() ) );

  */


  c1->Clear();

  ///////////////////////SIGMA  ///////////////////////////////
  TH2D* h2_axes3 = new TH2D( "axes3", "", 10, 0., xMax, 10, 0., 600. );
  h2_axes3->SetXTitle("Electron Beam Energy [MeV]");
  h2_axes3->SetYTitle("CeF_{3} Response Width [ADC]");
  h2_axes3->Draw("");

  gr_sigmaData->SetMarkerStyle(20);
  gr_sigmaData->SetMarkerSize(1.4);
  gr_sigmaData->SetMarkerColor(46);
  gr_sigmaData->Draw("p same");

  gr_sigmaSimul->SetMarkerStyle(20);
  gr_sigmaSimul->SetMarkerSize(1.4);
  gr_sigmaSimul->SetMarkerColor(38);
  gr_sigmaSimul->Draw("p same");

  gr_sigmaSimulPS->SetMarkerStyle(24);
  gr_sigmaSimulPS->SetMarkerSize(1.4);
  gr_sigmaSimulPS->SetMarkerColor(38);
  gr_sigmaSimulPS->Draw("p same");
  
  TLegend* leg = new TLegend(0.60, 0.4, 0.75, 0.2);
  leg->SetTextSize(0.038);
  leg->AddEntry(gr_sigmaData,"Data","p");
  leg->AddEntry(gr_sigmaSimul,"MC Sampling","p");
  leg->AddEntry(gr_sigmaSimulPS,"+ PhotoStat","p");

  leg->SetFillColor(0);
  leg->Draw("same");
  label_top->Draw("same");

  if(withHodo==1) c1->SaveAs( Form( "%s/sigmasWHodo.pdf", outputdir.c_str() ) );
  else  c1->SaveAs( Form( "%s/sigmasWOHodo.pdf", outputdir.c_str() ) );



  /*
  //////////////BGO ASYMMETRY//////////////////////////////
c1->Clear();

 TH2D* h2_axesb = new TH2D( "axesb", "", 10, 0., xMax, 10, -1,1);
  h2_axesb->SetXTitle("Electron Beam Energy [MeV]");
  h2_axesb->SetYTitle("BGO Asymmetry");
  h2_axesb->Draw("");

  gr_bgo_asymm->SetMarkerStyle(20);
  gr_bgo_asymm->SetMarkerSize(1.6);
  gr_bgo_asymm->SetMarkerColor(46);
  gr_bgo_asymm->Draw("p same");

  
  gr_resp_vs_energy_bgo_simul->SetMarkerStyle(20);
  gr_resp_vs_energy_bgo_simul->SetMarkerSize(1.6);
  gr_resp_vs_energy_bgo_simul->SetMarkerColor(38);
  gr_resp_vs_energy_bgo_simul->Draw("p same");


  TLegend* leg3 = new TLegend(0.60, 0.6, 0.855, 0.9);
  leg3->SetTextSize(0.038);
  leg3->AddEntry(gr_bgo_asymm,"data bgo","p");
  // leg3->AddEntry(f1_line_bgo,Form("offset=%.0f\n #pm %.0f\n",f1_line_bgo->GetParameter(0), f1_line_bgo->GetParError(0) ),"L");
  //leg3->AddEntry( (TObject*)0, Form("chi2/ndf=%.2f\n / %d",f1_line_bgo->GetChisquare(), f1_line_bgo->GetNDF() ), "");
   ;
    //leg3->AddEntry(gr_resp_vs_energy_bgo_simul,"simulation bgo","p");
    //leg3->AddEntry(f1_line_bgo_simul,Form("offset=%.0f\n #pm %.0f\n",f1_line_bgo_simul->GetParameter(0), f1_line_bgo_simul->GetParError(0) ),"L");
    //leg3->AddEntry( (TObject*)0, Form("chi2/ndf=%.2f\n / %d",f1_line_bgo_simul->GetChisquare(), f1_line_bgo_simul->GetNDF() ), "");
  leg3->SetFillColor(0);
  leg3->Draw("same");

   c1->SaveAs( Form( "%s/BGO_asym.pdf", outputdir.c_str() ) );
  */





  c1->Clear();

  //////////Linearity Deviation Plots ///////////////////////////////
  TH2D* h2_axes5 = new TH2D( "axes5", "", 10, 0., xMax, 10, 0.,0.5 );
  h2_axes5->SetXTitle("Electron Beam Energy [MeV]");
  h2_axes5->SetYTitle("E_{rec}/E_{beam}");
  h2_axes5->Draw("");

  gr_resp_vs_energy_dev->SetMarkerStyle(20);
  gr_resp_vs_energy_dev->SetMarkerSize(1.4);
  gr_resp_vs_energy_dev->SetMarkerColor(46);
  gr_resp_vs_energy_dev->Draw("p same");

  gr_resp_vs_energy_simul_dev->SetMarkerStyle(20);
  gr_resp_vs_energy_simul_dev->SetMarkerSize(1.4);
  gr_resp_vs_energy_simul_dev->SetMarkerColor(38);
  gr_resp_vs_energy_simul_dev->Draw("p same");

  
  TLegend* legd = new TLegend(0.60, 0.4, 0.75, 0.2);
  legd->SetTextSize(0.038);
  legd->AddEntry(gr_resp_vs_energy_dev,"Data","p");
  legd->AddEntry(gr_resp_vs_energy_simul_dev,"MC","p");

  legd->SetFillColor(0);
  legd->Draw("same");
  label_top->Draw("same");

  if(withHodo==1) c1->SaveAs( Form( "%s/deviationWHodo.pdf", outputdir.c_str() ) );
  else  c1->SaveAs( Form( "%s/deviationWOHodo.pdf", outputdir.c_str() ) );





  return 0;
}




/*
TF1* fitSingleElectronPeakBGO( const std::string& outputdir, const std::string& name, TTree* tree, int niter, float nSigma, std::string whattoProject, std::string cut ) {

  std::string histoName(Form("h1_%s", name.c_str()));
  TH1D* h1 = new TH1D(histoName.c_str(), "", 200, 10., 500.);

  tree->Project( histoName.c_str(), whattoProject.c_str(), cut.c_str());

  //tree->Project( histoName.c_str(),"bgo_corr[0]+bgo_corr[1]+bgo_corr[2]+bgo_corr[3]+bgo_corr[4]+bgo_corr[5]+bgo_corr[6]+bgo_corr[7]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 )");

  //  h1->Rebin(2);

  TF1* f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 10., 500.);
  f1->SetParameter(0, h1->Integral() );
  f1->SetParameter(1, h1->GetMean() );
  f1->SetParameter(2, h1->GetRMS() );
  f1->SetParError(1, h1->GetMeanError() );
  f1->SetParError(2, h1->GetRMSError() );

  doSingleFitBGO( h1, f1, outputdir, name, niter, nSigma );
  return f1;
}


void doSingleFitBGO( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name, int niter, float nSigma ) {

  h1->Fit( f1, "RQN" );
  f1->SetLineColor(kRed);
  f1->SetParLimits(1,1.,1000);
  f1->SetParLimits(2,1.,1000.);

  for( int iter=0; iter<niter; iter++ ) {

    float mean  = f1->GetParameter(1);
    float sigma = f1->GetParameter(2);
    float fitMin = mean - nSigma*sigma;
    //float fitMin = mean - nSigma*sigma/5.;
    if((mean-nSigma*sigma) <0.) {fitMin=20.;}
    float fitMax = mean + nSigma*sigma*1.5;
    f1->SetRange( fitMin, fitMax );
    if( iter==(niter-1) ){
     h1->Fit( f1, "RQ+" );
      } else{
      h1->Fit( f1, "RQN" ); }
  }
  TCanvas* c1 = new TCanvas("cX", "", 600, 600);
  c1->cd();
  h1->Draw();
  c1->SaveAs( Form("%s/fit_BGO_%s.pdf", outputdir.c_str(), name.c_str()) );
  delete c1;
}


/*
TF1* fitSingleElectronPeakBGOSimul( const std::string& outputdir, const std::string& name, TTree* tree, int niter, float nSigma ) {

  std::string histoName(Form("h1_%s", name.c_str()));
  TH1D* h1 = new TH1D(histoName.c_str(), "", 200, 0., 200.);

tree->Project( histoName.c_str(),"bgo_corr[0]+bgo_corr[1]+bgo_corr[2]+bgo_corr[3]+bgo_corr[4]+bgo_corr[5]+bgo_corr[6]+bgo_corr[7]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 )");
  TF1* f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 0., 200.);
  f1->SetParameter(0, h1->Integral() );
  f1->SetParameter(1, h1->GetMean() );
  f1->SetParameter(2, h1->GetRMS() );
  f1->SetParError(1, h1->GetMeanError() );
  f1->SetParError(2, h1->GetRMSError() );

  doSingleFitBGOSimul( h1, f1, outputdir, name, niter, nSigma );
  return f1;
}



void doSingleFitBGOSimul( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name, int niter, float nSigma ) {
  gStyle->SetOptFit(1);
  h1->Fit( f1, "RQN" );
   f1->SetLineColor(kRed);

   f1->SetParLimits(1,4.,60);
   f1->SetParLimits(2,2.,20.);

  for( int iter=0; iter<niter; iter++ ) {

    float mean  = f1->GetParameter(1); float sigma = f1->GetParameter(2);  float fitMin = mean - nSigma*sigma;
   if((mean-nSigma*sigma) <0. ) {fitMin=2.;}
    float fitMax = mean + nSigma*sigma;  f1->SetRange( fitMin, fitMax );
    if( iter==(niter-1) ){
     h1->Fit( f1, "RQ+" );
      }
    else{
      h1->Fit( f1, "RQN" ); }
  }
  //  TCanvas* c1 = new TCanvas("cX", "", 600, 600); c1->cd();
  //  h1->Draw();
  //  c1->SaveAs( Form("%s/fit_simulationBGO_%s.pdf", outputdir.c_str(), name.c_str()) );
  //  delete c1;
}

*/





ResoStruct getRespAndReso( TF1* f1, float energyS) {

  // float energy = energyS;
  float energyErrS = 0;

  float meanS = f1->GetParameter(1);
  float meanErrS = f1->GetParError(1);

  float rmsS = f1->GetParameter(2);
  float rmsErrS = f1->GetParError(2);

  // float resoS = 100.* rmsS/energyS; //in percent
  float resoS = 100.* rmsS/meanS; //in percent
  // float resoS = rmsS /173.097*3235.76; //in percent
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

ResoStruct addPhotoStatistics( ResoStruct rs, bool withHodo ) {

  // MC response is already in MeV, 0.49 is the number of p.e. per MeV
  // RM = 90% of energy, so assume 90% of energy (0.9*491) is deposited in central channel
  float nADC;

 if(withHodo==1){
       nADC = rs.resp/171.624*3235.37; //with hodo cut
 }  else{
       nADC = rs.resp/168.055*3223.23;
}


  float nPhotoElectrons = nADC/27.3;
  float poissonError = 100./sqrt( nPhotoElectrons ); // in percent
 
  float resoUnsmeared = rs.reso;
  float resoSmeared =  sqrt( rs.reso*rs.reso + poissonError*poissonError );

  rs.reso = resoSmeared;
  rs.reso_error = rs.reso_error * resoSmeared / resoUnsmeared; 

  return rs;
} 




TF1* fitSingleElectronPeak( const std::string& outputdir, const std::string& name, TTree* tree, int niter, float nSigma, bool withHodo ) {

  TH1D* h1 = new TH1D(name.c_str(), "", 200, 0., 6000.);
  h1 =  projectHodo(tree, name.c_str(), withHodo);

  gStyle->SetOptFit(2);
  TF1* f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 300., 6000.);
  f1->SetParameter(0, h1->Integral() );
  f1->SetParameter(1, h1->GetMean() );
  f1->SetParameter(2, h1->GetRMS() );

  f1->SetParError(1, h1->GetMeanError() );
  f1->SetParError(2, h1->GetRMSError() );

  FitTools::doSingleFit( h1, f1, outputdir, name, niter, nSigma );

  return f1;
}




TH1D* projectHodo ( TTree* tree, const std::string& name , bool withHodo){

  TH1D* h1 = new TH1D(name.c_str(), "", 200, 100., 6000.);

  if(withHodo==1){
    tree->Project(

		  //	  name.c_str() ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1&& pos_hodoClustY<2.5 && -2.5<pos_hodoClustY && pos_hodoClustX<2.5 && -2.5<pos_hodoClustX   )");
		 
		  //		  name.c_str() ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1&&  nFibres_hodoClustX > 1 && nFibres_hodoClustY>1 )");
		     name.c_str() ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 )");
}
  else{
    tree->Project( name.c_str() ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]","isSingleEle_scintFront");}


  return h1;
}




ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, const std::string& name, float energyS, bool withHodo ) {

  TH1D* h1 = new TH1D( name.c_str(), "", 500, 0., 500. );

  h1 =  projectHodoMC(tree, name.c_str(), withHodo);

  TF1* f1 = new TF1( Form("f1_%s", name.c_str()), "gaus", 5., 500.);

  FitTools::doSingleFit( h1, f1, outputdir, name,5.,1.5 );

  ResoStruct rs = getRespAndReso( f1, energyS );

  return rs;

}


TH1D* projectHodoMC ( TTree* tree, const std::string& name , bool withHodo){

  TH1D* h1 = new TH1D(name.c_str(), "", 200, 5., 500.);

  if(withHodo==1){
    tree->Project( name.c_str() ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 )");}
  else{
    tree->Project( name.c_str() ,"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]","isSingleEle_scintFront");}


  return h1;
}











  /*
  //////////////BGO //////////////////////////////
c1->Clear();

  TH2D* h2_axesb = new TH2D( "axesb", "", 10, 0., xMax, 10, 0., 600. );
  h2_axesb->SetXTitle("Electron Beam Energy [MeV]");
  h2_axesb->SetYTitle("BGO Response [ADC Counts]");
  h2_axesb->Draw("");

  gr_resp_vs_energy_bgo->SetMarkerStyle(20);
  gr_resp_vs_energy_bgo->SetMarkerSize(1.6);
  gr_resp_vs_energy_bgo->SetMarkerColor(46);
  gr_resp_vs_energy_bgo->Draw("p same");

  f1_line_bgo->Draw("L same");

  gr_resp_vs_energy_bgo_simul->SetMarkerStyle(20);
  gr_resp_vs_energy_bgo_simul->SetMarkerSize(1.6);
  gr_resp_vs_energy_bgo_simul->SetMarkerColor(38);
  gr_resp_vs_energy_bgo_simul->Draw("p same");

  TF1* f1_line_bgo_simul = new TF1("f1_line_bgo_simul", "[0] + [1]*x", 0., xMax );
  gr_resp_vs_energy_bgo_simul->Fit(f1_line_bgo_simul, "RN");
  f1_line_bgo_simul->SetLineWidth(1.);
  f1_line_bgo_simul->SetLineColor(38);
  f1_line_bgo_simul->Draw("L same");

  TLegend* leg3 = new TLegend(0.50, 0.2, 0.85, 0.4);
  leg3->SetTextSize(0.038);
  leg3->AddEntry(gr_resp_vs_energy_bgo,"data bgo","p");
  leg3->AddEntry(f1_line_bgo,Form("offset=%.0f\n #pm %.0f\n",f1_line_bgo->GetParameter(0), f1_line_bgo->GetParError(0) ),"L");
  leg3->AddEntry( (TObject*)0, Form("chi2/ndf=%.2f\n / %d",f1_line_bgo->GetChisquare(), f1_line_bgo->GetNDF() ), "");
  leg3->SetFillColor(0);
  leg3->AddEntry(gr_resp_vs_energy_bgo_simul,"simulation bgo","p");
  leg3->AddEntry(f1_line_bgo_simul,Form("offset=%.0f\n #pm %.0f\n",f1_line_bgo_simul->GetParameter(0), f1_line_bgo_simul->GetParError(0) ),"L");
  leg3->AddEntry( (TObject*)0, Form("chi2/ndf=%.2f\n / %d",f1_line_bgo_simul->GetChisquare(), f1_line_bgo_simul->GetNDF() ), "");
  leg3->SetFillColor(0);
  leg3->Draw("same");

  if(isIdeal==0){  c1->SaveAs( Form( "%s/resp_vs_energy_bgo.pdf", outputdir.c_str() ) );}
  */
