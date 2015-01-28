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

struct ResoStruct {

  float resp;
  float resp_error;

  float reso;
  float reso_error;

};

ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, const std::string& name, float energyS );
ResoStruct getRespAndReso( TF1* f1, float energyS );
float getRatioError( float num, float denom, float numErr, float denomErr );
ResoStruct addPhotoStatistics( ResoStruct rs );


int main( int argc, char* argv[] ) {

  TApplication* a = new TApplication("a",0,0);

  std::string outputdir = "H4Test/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  DrawTools::setStyle();
  TGaxis::SetMaxDigits(3);

  std::vector<std::string> runs; 
  std::vector<float> beamEnergy;
 
  std::vector<std::string> simulationBTF; 
  std::vector<float> beamEnergySimulationBTF; 

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
  
  
  simulationBTF.push_back("Simulation98");
  beamEnergySimulationBTF.push_back(98.);
  simulationBTF.push_back("Simulation147");
  beamEnergySimulationBTF.push_back(147.);
  simulationBTF.push_back("Simulation196");
  beamEnergySimulationBTF.push_back(196.);
  simulationBTF.push_back("Simulation295");
  beamEnergySimulationBTF.push_back(295.);
  simulationBTF.push_back("Simulation491");
  beamEnergySimulationBTF.push_back(491.);
  
  /*
  simulation.push_back("Simulation0.5");
  beamEnergySimulation.push_back(500.);
  simulation.push_back("Simulation1");
  beamEnergySimulation.push_back(1000.);
  simulation.push_back("Simulation5");
  beamEnergySimulation.push_back(5000.);
  */

  simulation.push_back("Simulation10");
  beamEnergySimulation.push_back(10000.);
  simulation.push_back("Simulation15");
  beamEnergySimulation.push_back(15000.);
  simulation.push_back("Simulation20");
  beamEnergySimulation.push_back(20000.);
  simulation.push_back("Simulation25");
  beamEnergySimulation.push_back(25000.);
  simulation.push_back("Simulation30");
  beamEnergySimulation.push_back(30000.);
  simulation.push_back("Simulation35");
  beamEnergySimulation.push_back(35000.);
  simulation.push_back("Simulation40");
  beamEnergySimulation.push_back(40000.);
  simulation.push_back("Simulation45");
  beamEnergySimulation.push_back(45000.);
  simulation.push_back("Simulation50");
  beamEnergySimulation.push_back(50000.);


  TGraphErrors* gr_resp_vs_energy = new TGraphErrors(0);
  TGraphErrors* gr_reso_vs_energy = new TGraphErrors(0);

  TGraphErrors* gr_resp_vs_energy_simul = new TGraphErrors(0);
  TGraphErrors* gr_reso_vs_energy_simul = new TGraphErrors(0);

  TGraphErrors* gr_resp_BTF_simul = new TGraphErrors(0);

  TGraphErrors* gr_resp_total_simul = new TGraphErrors(0);

  float xMax = 55000;

  TFile* fileMean = TFile::Open(Form("analysisTrees_V03/Reco_%s.root",runs[4].c_str()));
  TTree* treeMean = (TTree*)fileMean->Get("recoTree");
 
  //and for CeF3:
  TF1* energyfuncC = FitTools::fitSingleElectronPeak( outputdir, runs[4], treeMean,5,1.5 );

  float energy491C= beamEnergy[4];
  float adcEnergyC = energyfuncC->GetParameter(1);

  //for the simulated Cef3:
  TFile* energyfileS = TFile::Open(Form("OriginalSimulationData/H4/Reco_%s.root", simulation[4].c_str()));
  //TFile* energyfileS = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(), simulation[4].c_str()));
  TTree* energytreeS = (TTree*)energyfileS->Get("recoTree");
  
  ResoStruct energyrsS = getResponseResolutionMC( outputdir, energytreeS, simulation[4], beamEnergySimulation[4]);
  
  float energy491S = energyrsS.resp;



 for( unsigned i=0; i<runs.size(); ++i ) {

    TFile* file = TFile::Open(Form("analysisTrees_V03/Reco_%s.root",  runs[i].c_str()));
    TTree* tree = (TTree*)file->Get("recoTree");

    TF1* thisFunc = FitTools::fitSingleElectronPeak( outputdir, runs[i], tree,5, 1.5 );

    float energy = beamEnergy[i];
    float energyErr = 5.;
 
    float mean = thisFunc->GetParameter(1);
    float meanErr = thisFunc->GetParError(1);

    float rms = thisFunc->GetParameter(2);
    float rmsErr = thisFunc->GetParError(2);

    float reso = 100.* rms/mean ;
    float resoErr = 100* getRatioError( rms, mean, rmsErr, meanErr);

    gr_resp_vs_energy->SetPoint( i, energy/1000., mean/1000. );
    gr_resp_vs_energy->SetPointError( i, energyErr/100., meanErr/1000. );

    gr_reso_vs_energy->SetPoint( i, energy, reso );
    gr_reso_vs_energy->SetPointError( i, energyErr, resoErr );

    std::cout << "Resolution data = " << reso << std::endl;

    //for the simulated Cef3 of BTF:
    TFile* fileS = TFile::Open(Form("OriginalSimulationData/H4/Reco_%s.root", simulationBTF[i].c_str()));
    TTree* treeS = (TTree*)fileS->Get("recoTree");
    
    ResoStruct rsS = getResponseResolutionMC( outputdir, treeS, simulationBTF[i], beamEnergySimulationBTF[i]);
    
   
    gr_resp_BTF_simul->SetPoint( i, beamEnergySimulationBTF[i]/1000., rsS.resp/1000. );
    gr_resp_BTF_simul->SetPointError( i, 0., rsS.resp_error/1000. );
    

    //all together now//
    gr_resp_total_simul->SetPoint( i, beamEnergySimulationBTF[i]/1000.,  rsS.resp/1000.);
    gr_resp_total_simul->SetPointError( i, 0., rsS.resp_error/1000. );

  }



 for( unsigned i=0; i<simulation.size(); ++i ) {


    ///////////////// (1x1) Shashlik ("real setup") //////////////////////////////
    TFile* fileS = TFile::Open(Form("OriginalSimulationData/H4/Reco_%s.root", simulation[i].c_str()));


    TTree* treeS = (TTree*)fileS->Get("recoTree");

    ResoStruct rs = getResponseResolutionMC( outputdir, treeS, simulation[i], beamEnergySimulation[i] );
    ResoStruct rs_ps = addPhotoStatistics( rs );

    gr_resp_vs_energy_simul->SetPoint( i, beamEnergySimulation[i]/1000., rs.resp/1000. );
    //   gr_resp_vs_energy_simul->SetPoint( i, beamEnergySimulation[i]/1000., rs.resp * adcEnergyC/energy491S/1000. );
    gr_resp_vs_energy_simul->SetPointError( i,0, rs.resp_error/1000.);
    //    gr_resp_vs_energy_simul->SetPointError( i,0, rs.resp_error/1000.);

    gr_reso_vs_energy_simul->SetPoint( i, beamEnergySimulation[i]/1000., rs_ps.reso );
    //    gr_reso_vs_energy_simul->SetPoint( i, beamEnergySimulation[i], rs.reso );
    gr_reso_vs_energy_simul->SetPointError( i,0,  rs_ps.reso_error );

    std::cout << "Resolution = " << rs.reso << std::endl;
    std::cout << "Resolution.PE = " << rs_ps.reso << std::endl;


    gr_resp_total_simul->SetPoint( i+5, beamEnergySimulation[i]/1000. ,  rs.resp/1000.);
    gr_resp_total_simul->SetPointError( i+5, 0., rs.resp_error/1000. );

  }

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();


  TH2D* h2_axes = new TH2D( "axes", "", 10, 0., 51, 10, 0., 20. );
  h2_axes->SetXTitle("Electron Beam Energy [GeV]");
  h2_axes->SetYTitle("CeF3 Response [GeV]");
  h2_axes->Draw("");




  gr_resp_BTF_simul->SetMarkerStyle(20);
  gr_resp_BTF_simul->SetMarkerSize(1.6);
  gr_resp_BTF_simul->SetMarkerColor(46);
  gr_resp_BTF_simul->Draw("p same");

  gr_resp_vs_energy_simul->SetMarkerStyle(20);
  gr_resp_vs_energy_simul->SetMarkerSize(1.6);
  gr_resp_vs_energy_simul->SetMarkerColor(38);
  gr_resp_vs_energy_simul->Draw("p same");


  gr_resp_total_simul->SetMarkerStyle(25);
  gr_resp_total_simul->SetMarkerSize(1.6);
  gr_resp_total_simul->SetMarkerColor(kBlack);
  gr_resp_total_simul->Draw("p same");

  TF1* f1_line = new TF1("line", "[0] + [1]*x", 0., xMax/1000. );
  gr_resp_BTF_simul->Fit(f1_line, "RN");
  f1_line->SetLineColor(46);
  f1_line->SetLineWidth(1.);
  f1_line->Draw("L same");

  TF1* f1_lines = new TF1("lines", "[0] + [1]*x", 0., xMax/1000. );
  gr_resp_vs_energy_simul->Fit(f1_lines,"RN");
  f1_lines->SetLineWidth(1.);
  f1_lines->SetLineColor(38);
  f1_lines->Draw("L same");

  TF1* f1_lines2 = new TF1("lines", "[0] + [1]*x", 0., xMax/1000. );
  gr_resp_total_simul->Fit(f1_lines2,"RN");
  f1_lines2->SetLineWidth(1.);
  f1_lines2->SetLineColor(kBlack);
  f1_lines2->Draw("L same");

    TLegend* leg0 = new TLegend(0.18, 0.93, 0.6, 0.6);
  // TLegend* leg0 = new TLegend(0.45, 0.2, 0.9, 0.4);  
leg0->SetTextSize(0.038);
  /*
  leg0->AddEntry(gr_resp_vs_energy,"uncorrected data","p");
  leg0->AddEntry(f1_line,Form("offset=%.0f\n #pm %.0f\n",f1_line->GetParameter(0)*1000, f1_line->GetParError(0)*100 ),"L");
  leg0->AddEntry( (TObject*)0, Form("chi2/ndf=%.2f\n / %d",f1_line->GetChisquare(), f1_line->GetNDF() ), "");
  */

  leg0->AddEntry(gr_resp_BTF_simul,"MC (1x1) BTF","p");
  leg0->AddEntry(f1_line,Form("Offset = %.1f\n #pm %.2f\n MeV",f1_line->GetParameter(0)*1000., f1_line->GetParError(0)*1000 ),"L");
  leg0->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.0f\n / %d",f1_line->GetChisquare(), f1_line->GetNDF() ), "");

  leg0->AddEntry(gr_resp_vs_energy_simul,"MC (1x1) H4","p");
  leg0->AddEntry(f1_lines2,Form("Offset = %.1f\n #pm %.1f\n MeV",f1_lines->GetParameter(0)*1000., f1_lines->GetParError(0)*1000 ),"L");
  leg0->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.0f\n / %d",f1_lines->GetChisquare(), f1_lines->GetNDF() ), "");


  leg0->AddEntry(gr_resp_total_simul,"MC (1x1) BTF + H4","p");
  leg0->AddEntry(f1_lines2,Form("Offset = %.1f\n #pm %.1f\n MeV",f1_lines2->GetParameter(0)*1000., f1_lines2->GetParError(0)*1000 ),"L");
  leg0->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.0f\n / %d",f1_lines2->GetChisquare(), f1_lines2->GetNDF() ), "");


  leg0->SetFillColor(0);
  leg0->Draw("same");

  TPaveText* label_top2 = new TPaveText();
  label_top2 = DrawTools::getLabelTop("Single Electron Beam");
  label_top2->Draw("same");

  c1->SaveAs( Form( "%s/resp_vs_energy.pdf", outputdir.c_str() ) );
  
  c1->Clear();








 ///////////////////////RESOLUTION///////////////////////////////
 TH2D* h2_axes2 = new TH2D( "axes", "", 10, 0., xMax/1000.-3, 10, 0., 20. );
 h2_axes2->SetXTitle("Electron Beam Energy [GeV]");
 h2_axes2->SetYTitle("Energy Resolution [%]");
 h2_axes2->Draw("");
 
 /*
 // Data (1x1)
 gr_reso_data->SetMarkerStyle(20);
 gr_reso_data->SetMarkerSize(1.6);
 gr_reso_data->SetMarkerColor(46);
 gr_reso_data->Draw("p same");
 
 //TF1 *fun= new TF1("fun", "sqrt([0]*[0]/(x))",50, xMax);
 TF1 *fun= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]/(x*x))",50, xMax);
 fun->SetParameter(1,50.);
 gr_reso_data->Fit(fun,"RN");
 fun->SetLineWidth(1.);
 fun->SetLineColor(46);
 // fun->Draw("L same");
 */
 
 // MC (1x1)
 gr_reso_vs_energy_simul->SetMarkerStyle(20);
 gr_reso_vs_energy_simul->SetMarkerSize(1.6);
 gr_reso_vs_energy_simul->SetMarkerColor(38);
 gr_reso_vs_energy_simul->Draw("p same");
 
 TF1 *fun1= new TF1("fun1",  "sqrt([0]*[0]/x+[1]*[1])",0.2, xMax/1000.);
 // TF1 *fun1= new TF1("fun1",  "sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",0.2, xMax/1000.);
 //fun1->SetParameter(0,50.);
 fun1->SetParameter(2,50.);
 gr_reso_vs_energy_simul->Fit(fun1,"RN");
 fun1->SetLineWidth(1.);
 fun1->SetLineColor(38);
 fun1->Draw("L same");


 
 TLegend* leg4 = new TLegend(0.3, 0.55, 0.78, 0.9);
 leg4->SetTextSize(0.038);
 leg4->AddEntry(gr_reso_vs_energy_simul,"MC (1x1) H4","p");
 leg4->AddEntry(fun1,Form("S = %.2f\n #pm %.2f\n %s / #sqrt{E [GeV]}",fun1->GetParameter(0), (fun1->GetParError(0)),"%" ),"L");
 leg4->AddEntry( (TObject*)0,Form("C =   %.2f\n #pm %.2f\n %s",(fun1->GetParameter(1)), (fun1->GetParError(1)),"%" ),"");
 // leg4->AddEntry( (TObject*)0,Form("N =   %.2f\n #pm %.2f\n %s / E [GeV]",(fun1->GetParameter(2)), (fun1->GetParError(2)),"%" ),"");
 leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun1->GetChisquare(), fun1->GetNDF() ), "");

 leg4->SetFillColor(0);
 leg4->Draw("same");
 
 TPaveText* label_low = new TPaveText(0.165,0.175,0.5,0.21, "brNDC");
 label_low->SetFillColor(kWhite);
 label_low->SetTextSize(0.038);
 label_low->SetTextAlign(11); // align right
 label_low->SetTextFont(62);
 label_low->AddText( "W-CeF_{3} Single Tower");
 // label_low->Draw("same");
 
  label_top2->Draw("same");

 
 c1->SaveAs( Form( "%s/resolution.pdf", outputdir.c_str() ) );








  return 0;

}















ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, const std::string& name, float energyS ) {

  TH1D* h1 = new TH1D( name.c_str(), "", 500, 0., 25000. );

   tree->Project( name.c_str(),"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]","isSingleEle_scintFront" );

  // tree->Project( name.c_str(),"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]","(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 )" );

  TF1* f1 = new TF1( Form("f1_%s", name.c_str()), "gaus", 50., 25000.);

  FitTools::doSingleFit( h1, f1, outputdir, name,5.,1.5 );

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


ResoStruct addPhotoStatistics( ResoStruct rs ) {

  // MC response is already in MeV, 0.49 is the number of p.e. per MeV
  // RM = 90% of energy, so assume 90% of energy (0.9*491) is deposited in central channel
  //   float nADC = rs.resp/169.445*3224.36;
  //  float nADC = rs.resp/173.097*3235.76;
  float  nADC = rs.resp/169.445*3224.36;
   float nPhotoElectrons = nADC/27.3;
  //float nADC = rs.resp/185.*3200.;
  //float nPhotoElectrons = nADC/35.;

  float poissonError = 100./sqrt( nPhotoElectrons ); // in percent
 
  rs.reso = sqrt( rs.reso*rs.reso + poissonError*poissonError );
 
  return rs;

} 
