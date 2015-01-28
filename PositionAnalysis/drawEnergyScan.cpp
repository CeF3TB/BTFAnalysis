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

#include "interface/DrawTools.h"
#include "interface/FitTools.h"

#include "TApplication.h"

//Draw Linearty and Resolution for conditions of Up vs Down, Left vs Right


TF1* fitSingleElectronPeakUoDLoR( const std::string& outputdir, const std::string& name, TTree* tree, int niter, float nSigma, int isUp, bool isLeft );
TF1* fitSingleElectronPeakBGO( const std::string& outputdir, const std::string& name, TTree* tree, int niter, float nSigma );
TH1* meanOfSingleElectronPeakBGO( const std::string& outputdir, const std::string& name, TTree* tree );
//left and right are switched, only correct in legends

bool checkUpDown=true;
bool checkLeftRight=true;

int main( int argc, char* argv[] ) {

  TApplication* a = new TApplication("a",0,0);


  std::string tag="V00";
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }

  std::cout << "-> Using tag: " << tag << std::endl;

  std::string outputdir = "EnergyScanPlots/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  std::string outputdir2 = "EnergyScanPlots2/";
  std::string mkdir_command2 = "mkdir -p " + outputdir2;
  system( mkdir_command2.c_str() );

  DrawTools::setStyle();
 
  std::vector<std::string> runs; 
  std::vector<float> beamEnergy; 


  runs.push_back("BTF_314_20140503-024715_beam");
  beamEnergy.push_back(98.3);

  runs.push_back("BTF_308_20140503-002534_beam");
  beamEnergy.push_back(147.4);

  runs.push_back("BTF_293_20140502-180258_beam");
  beamEnergy.push_back(196.5);

  //runs.push_back("BTF_286_287");
  runs.push_back("BTF_286_20140502-153528_beam");
  beamEnergy.push_back(294.8);

  runs.push_back("BTF_259_20140502-012847_beam");
  beamEnergy.push_back(491.4);

  //runs.push_back("BTF_246_20140501-212512_beam");
  //beamEnergy.push_back(491.4);

  TGraphErrors* gr_resp_vs_energy = new TGraphErrors(0);
  TGraphErrors* gr_reso_vs_energy = new TGraphErrors(0);


  TGraphErrors* gr_resp_vs_energy_up = new TGraphErrors(0);
  TGraphErrors* gr_resp_vs_energy_down = new TGraphErrors(0);

  TGraphErrors* gr_resp_vs_energy_left = new TGraphErrors(0);
  TGraphErrors* gr_resp_vs_energy_right = new TGraphErrors(0);

  TGraphErrors* gr_resp_vs_energy_bgo = new TGraphErrors(0);

  TGraphErrors* gr_resp_vs_energy_bgo_m = new TGraphErrors(0);

  //491MeV=peak pos of run 259 in adc

  TFile* fileMean = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(), runs[4].c_str()));
    TTree* treeMean = (TTree*)fileMean->Get("recoTree");
    TF1* energyfunc = FitTools::fitSingleElectronPeak( outputdir, runs[4], treeMean );

    float energy491= beamEnergy[4];
    float adcEnergy = energyfunc->GetParameter(1);
    float adcEnergyErr = energyfunc->GetParError(1);



  for( unsigned i=0; i<runs.size(); ++i ) {

    TFile* file = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(), runs[i].c_str()));
    TTree* tree = (TTree*)file->Get("recoTree");

    TF1* thisFunc = FitTools::fitSingleElectronPeak( outputdir, runs[i], tree );

    float energy = beamEnergy[i];
    float energyErr = (energy>200.) ? 0.01*energy : 0.05*energy;
    //float energyErr = (energy>200.) ? 0.005*energy : 0.05*energy;

    float mean = thisFunc->GetParameter(1);
    float meanErr = thisFunc->GetParError(1);

    float rms = thisFunc->GetParameter(2);
    float rmsErr = thisFunc->GetParError(2);

    float reso = 100.* rms/mean; //in percent
    float resoErr = 100.* sqrt( rmsErr*rmsErr/(mean*mean) + rms*rms*meanErr*meanErr/(mean*mean*mean*mean) );

    // float reso = 100.* rms/energy* energy491 / adcEnergy ;
    //float resoErr = 100.* sqrt( rmsErr*rmsErr* energy491/adcEnergy*energy491/adcEnergy/(energy*energy) + rms*rms* energy491 / adcEnergy* energy491 / adcEnergy* energyErr*energyErr/(energy*energy*energy*energy)+    (rms*rms* energy491*0.01 / adcEnergy* energy491*0.01 / adcEnergy)/(energy*energy) + 	(rms*rms* energy491  *energy491*adcEnergyErr*adcEnergyErr)/(adcEnergy*adcEnergy*adcEnergy*adcEnergy*energy*energy) );

    gr_resp_vs_energy->SetPoint( i, energy, mean );
    gr_resp_vs_energy->SetPointError( i, energyErr, meanErr );

    gr_reso_vs_energy->SetPoint( i, energy, reso );
    gr_reso_vs_energy->SetPointError( i, energyErr, resoErr );



    //////BGO /////////////
    TF1* funcBGO = fitSingleElectronPeakBGO(outputdir2, runs[i],tree,4, 0.9);

    float meanBGO = funcBGO->GetParameter(1);
    float meanBGOerr = funcBGO->GetParError(1);

    gr_resp_vs_energy_bgo->SetPoint(i, energy, meanBGO);
    gr_resp_vs_energy_bgo->SetPointError(i, energyErr, meanBGOerr);


    /////BGO mean of histo not fitted ////////////////
    TH1* histoBGO = meanOfSingleElectronPeakBGO(outputdir2, runs[i], tree );

    float meanBGOm = histoBGO->GetMean();

    float meanBGOerrm = histoBGO->GetMeanError();

    gr_resp_vs_energy_bgo_m->SetPoint(i, energy, meanBGOm);
    gr_resp_vs_energy_bgo_m->SetPointError(i, energyErr, meanBGOerrm);




if(checkUpDown==1){
    TF1* funcUp =fitSingleElectronPeakUoDLoR(outputdir, runs[i], tree, 4, 1.5, 1,0 );
    TF1* funcDown =fitSingleElectronPeakUoDLoR(outputdir, runs[i], tree, 4, 1.5,  0 ,0);

    float meanUp = funcUp->GetParameter(1);
    float meanDown = funcDown->GetParameter(1);
    float meanUperr = funcUp->GetParError(1);
    float meanDownerr = funcDown->GetParError(1);
    std::cout << "MeanUP = " << meanUp << std::endl;
    std::cout << "MeanDown = " << meanDown << std::endl;

    gr_resp_vs_energy_up->SetPoint(i, energy, meanUp);
    gr_resp_vs_energy_down->SetPoint(i, energy, meanDown);
    gr_resp_vs_energy_up->SetPointError(i, energyErr, meanUperr);
    gr_resp_vs_energy_down->SetPointError(i, energyErr, meanDownerr);
    }

if(checkLeftRight==1){
    TF1* funcLeft =fitSingleElectronPeakUoDLoR( outputdir, runs[i], tree, 4, 1.5, -1,1 );
    TF1* funcRight =fitSingleElectronPeakUoDLoR( outputdir, runs[i], tree, 4, 1.5, -1 ,0);

    float meanLeft = funcLeft->GetParameter(1);
    float meanRight = funcRight->GetParameter(1);
    float meanLefterr = funcLeft->GetParError(1);
    float meanRighterr = funcRight->GetParError(1);

    gr_resp_vs_energy_left->SetPoint(i, energy, meanLeft);
    gr_resp_vs_energy_right->SetPoint(i, energy, meanRight);
    gr_resp_vs_energy_left->SetPointError(i, energyErr, meanLefterr);
    gr_resp_vs_energy_right->SetPointError(i, energyErr, meanRighterr);
    }



  }


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  float xMax = 550.;

  TH2D* h2_axes = new TH2D( "axes", "", 10, 0., xMax, 10, 0., 4000. );
  h2_axes->SetXTitle("Electron Beam Energy [MeV]");
  h2_axes->SetYTitle("CeF3 Response [ADC Counts]");
  h2_axes->Draw("");

  gr_resp_vs_energy->SetMarkerStyle(20);
  gr_resp_vs_energy->SetMarkerSize(1.6);
  gr_resp_vs_energy->SetMarkerColor(46);
  gr_resp_vs_energy->Draw("p same");

  TF1* f1_line = new TF1("line", "[0] + [1]*x", 0., xMax );
  //f1_line->FixParameter(0, 0.);
  gr_resp_vs_energy->Fit(f1_line, "RN");
  f1_line->SetLineWidth(1.);
  f1_line->Draw("L same");

  std::cout << "Chi2: " << f1_line->GetChisquare() << std::endl;
  std::cout << "Prob: " << TMath::Prob(f1_line->GetChisquare(), 2 ) << std::endl;

  TLegend* leg0 = new TLegend(0.6, 0.2, 0.9, 0.4);
  leg0->AddEntry(gr_resp_vs_energy,"data","p");
  leg0->AddEntry(f1_line,Form("offset=%d #pm %d",int(f1_line->GetParameter(0)),int( f1_line->GetParError(0)) ),"L");
  //leg3->AddEntry( (TObject*)0 ,Form("slope=%f",f1_line_bgo->GetParameter(1)),"");
  leg0->AddEntry( (TObject*)0, Form("chi2/ndf=%.2f\n / %d",f1_line->GetChisquare(), f1_line->GetNDF() ), "");
  leg0->SetFillColor(0);
  leg0->Draw("same");


  //c1->SaveAs( Form( "%s/resp_vs_energy.eps", outputdir.c_str() ) );
  //c1->SaveAs( Form( "%s/resp_vs_energy.png", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/resp_vs_energy.pdf", outputdir.c_str() ) );


  c1->Clear();


  //////////////////////////////////////////////////////////
  //////////////UP and DOWN comparison /////////////////////
  if(checkUpDown==1){
c1->Clear();
  h2_axes->Draw("");

  gr_resp_vs_energy_up->SetMarkerStyle(21);
  gr_resp_vs_energy_up->SetMarkerSize(1.6);
  gr_resp_vs_energy_up->SetMarkerColor(46);
  gr_resp_vs_energy_up->Draw("p same");

  gr_resp_vs_energy_down->SetMarkerStyle(20);
  gr_resp_vs_energy_down->SetMarkerSize(1.6);
  gr_resp_vs_energy_down->SetMarkerColor(38);
  gr_resp_vs_energy_down->Draw("p same");


  TF1* f1_line_up = new TF1("line", "[0] + [1]*x", 0., xMax );
  //f1_line->FixParameter(0, 0.);
  gr_resp_vs_energy_up->Fit(f1_line_up, "RN");
  f1_line_up->SetLineWidth(1.);
  f1_line_up->SetLineColor(46);
  f1_line_up->Draw("L same");

  TF1* f1_line_down = new TF1("line", "[0] + [1]*x", 0., xMax );
  //f1_line->FixParameter(0, 0.);
  gr_resp_vs_energy_down->Fit(f1_line_down, "RN");
  f1_line_down->SetLineWidth(1.);
  f1_line_down->SetLineColor(38);
  f1_line_down->Draw("L same");


  TLegend* leg2 = new TLegend(0.6, 0.2, 0.9, 0.4);
leg2->SetTextSize(0.038);
  leg2->AddEntry(gr_resp_vs_energy_up,"data up","p");  
  leg2->AddEntry(f1_line_up,Form("offset=%d #pm %d",int(f1_line_up->GetParameter(0)),int( f1_line_up->GetParError(0))),"L");
  // leg2->AddEntry( (TObject*)0 ,Form("slope=%f",f1_line_up->GetParameter(1)),"");
  leg2->AddEntry( (TObject*)0, Form("chi2/ndf=%.2f\n / %d",f1_line_up->GetChisquare(), f1_line_up->GetNDF() ), "");
  leg2->AddEntry(gr_resp_vs_energy_down,"data down","p");
  leg2->AddEntry(f1_line_down,Form("offset=%d #pm %d",int(f1_line_down->GetParameter(0)),int( f1_line_down->GetParError(0))),"L");
  // leg2->AddEntry( (TObject*)0 ,Form("slope=%f",f1_line_down->GetParameter(1)),"");
  leg2->AddEntry( (TObject*)0, Form("chi2/ndf=%.2f\n / %d",f1_line_down->GetChisquare() , f1_line_down->GetNDF()), "");

  leg2->SetFillColor(0);
  leg2->Draw("same");

  c1->SaveAs( Form( "%s/resp_vs_energy_UpDown.pdf", outputdir.c_str()) );

  }

  //////////////////////////////////////////////////////////
  //////////////LEFT and RIGHT comparison /////////////////////
  if(checkLeftRight==1){

c1->Clear();

  h2_axes->Draw("");

  gr_resp_vs_energy_left->SetMarkerStyle(21);
  gr_resp_vs_energy_left->SetMarkerSize(1.6);
  gr_resp_vs_energy_left->SetMarkerColor(46);
  gr_resp_vs_energy_left->Draw("p same");

  gr_resp_vs_energy_right->SetMarkerStyle(20);
  gr_resp_vs_energy_right->SetMarkerSize(1.6);
  gr_resp_vs_energy_right->SetMarkerColor(38);
  gr_resp_vs_energy_right->Draw("p same");


  TF1* f1_line_left = new TF1("line", "[0] + [1]*x", 0., xMax );
  //f1_line->FixParameter(0, 0.);
  gr_resp_vs_energy_left->Fit(f1_line_left, "RN");
  f1_line_left->SetLineWidth(1.);
  f1_line_left->SetLineColor(46);
  f1_line_left->Draw("L same");

  TF1* f1_line_right = new TF1("line", "[0] + [1]*x", 0., xMax );
  //f1_line->FixParameter(0, 0.);
  gr_resp_vs_energy_right->Fit(f1_line_right, "RN");
  f1_line_right->SetLineWidth(1.);
  f1_line_right->SetLineColor(38);
  f1_line_right->Draw("L same");


  TLegend* leg3 = new TLegend(0.55, 0.2, 0.9, 0.45);
leg3->SetTextSize(0.038);
  leg3->AddEntry(gr_resp_vs_energy_left,"data right","p");
  leg3->AddEntry(f1_line_left,Form("offset=%d #pm %d",int(f1_line_left->GetParameter(0)),int( f1_line_left->GetParError(0)) ),"L");
  //leg3->AddEntry( (TObject*)0 ,Form("slope=%f",f1_line_left->GetParameter(1)),"");
  leg3->AddEntry( (TObject*)0, Form("chi2/ndf=%.2f\n / %d",f1_line_left->GetChisquare(), f1_line_left->GetNDF() ), "");
  leg3->AddEntry(gr_resp_vs_energy_right,"data left","p");
  leg3->AddEntry(f1_line_right,Form("offset=%d #pm %d",int(f1_line_right->GetParameter(0)),int( f1_line_right->GetParError(0))),"L");
  //leg3->AddEntry( (TObject*)0 ,Form("slope=%f",f1_line_right->GetParameter(1)),"");
  leg3->AddEntry( (TObject*)0, Form("chi2/ndf=%.2f\n / %d",f1_line_right->GetChisquare(), f1_line_right->GetNDF() ), "");
  leg3->SetFillColor(0);
  leg3->Draw("same");

  c1->SaveAs( Form( "%s/resp_vs_energy_LeftRight.pdf", outputdir.c_str() ) );

  }

  //////////////BGO //////////////////////////////
c1->Clear();

  TH2D* h2_axesb = new TH2D( "axesb", "", 10, 0., xMax, 10, 0., 700. );
  h2_axesb->SetXTitle("Electron Beam Energy [MeV]");
  h2_axesb->SetYTitle("BGO Response [ADC Counts]");
  h2_axesb->Draw("");

  gr_resp_vs_energy_bgo->SetMarkerStyle(20);
  gr_resp_vs_energy_bgo->SetMarkerSize(1.6);
  gr_resp_vs_energy_bgo->SetMarkerColor(46);
  gr_resp_vs_energy_bgo->Draw("p same");

  TF1* f1_line_bgo = new TF1("line", "[0] + [1]*x", 0., xMax );
  //f1_line->FixParameter(0, 0.);
  gr_resp_vs_energy_bgo->Fit(f1_line_bgo, "RN");
  f1_line_bgo->SetLineWidth(1.);
  f1_line_bgo->SetLineColor(46);
  f1_line_bgo->Draw("L same");

  TLegend* leg3 = new TLegend(0.6, 0.2, 0.9, 0.4);
  leg3->AddEntry(gr_resp_vs_energy_bgo,"data bgo","p");
  leg3->AddEntry(f1_line_bgo,Form("offset=%d #pm %d",int(f1_line_bgo->GetParameter(0)),int( f1_line_bgo->GetParError(0)) ),"L");
  //leg3->AddEntry( (TObject*)0 ,Form("slope=%f",f1_line_bgo->GetParameter(1)),"");
  leg3->AddEntry( (TObject*)0, Form("chi2/ndf=%.2f\n / %d",f1_line_bgo->GetChisquare(), f1_line_bgo->GetNDF() ), "");
  leg3->SetFillColor(0);
  leg3->Draw("same");

  c1->SaveAs( Form( "%s/resp_vs_energy_bgo.pdf", outputdir.c_str() ) );





  //////////////BGO mean not fitted//////////////////////////////
  if(false){
c1->Clear();

  h2_axesb->Draw("");

  gr_resp_vs_energy_bgo_m->SetMarkerStyle(20);
  gr_resp_vs_energy_bgo_m->SetMarkerSize(1.6);
  gr_resp_vs_energy_bgo_m->SetMarkerColor(46);
  gr_resp_vs_energy_bgo_m->Draw("p same");

  TF1* f1_line_bgo2 = new TF1("line", "[0] + [1]*x", 0., xMax );
  gr_resp_vs_energy_bgo_m->Fit(f1_line_bgo2, "RN");
  f1_line_bgo2->SetLineWidth(1.);
  f1_line_bgo2->SetLineColor(46);
  f1_line_bgo2->Draw("L same");

  TLegend* legb = new TLegend(0.6, 0.2, 0.9, 0.4);
legb->SetTextSize(0.038);
  legb->AddEntry(gr_resp_vs_energy_bgo_m,"data bgo","p");
  legb->AddEntry(f1_line_bgo2,Form("offset=%d #pm %d",int(f1_line_bgo2->GetParameter(0)),int( f1_line_bgo2->GetParError(0)) ),"L");
  //leg3->AddEntry( (TObject*)0 ,Form("slope=%f",f1_line_bgo->GetParameter(1)),"");
  legb->AddEntry( (TObject*)0, Form("chi2/ndf=%.2f\n / %d",f1_line_bgo2->GetChisquare(), f1_line_bgo2->GetNDF() ), "");
  legb->SetFillColor(0);
  legb->Draw("same");

  c1->SaveAs( Form( "%s/resp_vs_energy_bgo2.pdf", outputdir2.c_str() ) );
  }

  ///////////////////////RESOLUTION///////////////////////////////
  TH2D* h2_axes2 = new TH2D( "axes", "", 10, 0., xMax, 10, 0., 55. );
  h2_axes2->SetXTitle("Electron Beam Energy [MeV]");
  h2_axes2->SetYTitle("CeF3 Resolution [%]");
  h2_axes2->Draw("");

  gr_reso_vs_energy->SetMarkerStyle(20);
  gr_reso_vs_energy->SetMarkerSize(1.6);
  gr_reso_vs_energy->SetMarkerColor(46);
  gr_reso_vs_energy->Draw("p same");

  TF1 *fun= new TF1("fun", "sqrt([0]*[0]+[1]*[1]/(x))",50, xMax);
  //TF1 *fun= new TF1("fun", "sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))",50, xMax);

  //fun->SetParameter(1,100);
  //fun->SetParameter(0,5);
   gr_reso_vs_energy->Fit(fun,"RN");
  fun->SetLineWidth(1.);
  fun->Draw("L same");

  TLegend* leg4 = new TLegend(0.35, 0.65, 0.9, 0.9);
leg4->SetTextSize(0.038);
leg4->AddEntry(gr_reso_vs_energy,"data (nominal beam energy)","p");
//leg4->AddEntry(gr_reso_vs_energy,"data (CeF3 energy)","p");

  leg4->AddEntry(fun,Form("S = %.2f\n #pm %.2f\n ",1./sqrt(1000)*(fun->GetParameter(1)), 1./sqrt(1000)*(fun->GetParError(1)) ),"L");
  leg4->AddEntry((TObject*)0,Form("C = %.2f\n #pm %.2f\n  ", fun->GetParameter(0), fun->GetParError(0) ), "");
  //leg4->AddEntry((TObject*)0,Form("c=%.2f\n #pm %.2f\n ",fun->GetParameter(2), fun->GetParError(2) ),"");
  leg4->AddEntry( (TObject*)0, Form("chi2/ndf=%.2f\n / %d",fun->GetChisquare(), fun->GetNDF() ), "");
  leg4->SetFillColor(0);
  leg4->Draw("same");

  //c1->SaveAs( Form( "%s/reso_vs_energy.eps", outputdir.c_str() ) );
  //c1->SaveAs( Form( "%s/reso_vs_energy.png", outputdir.c_str() ) );
  //c1->SaveAs( Form( "%s/reso_vs_energy_cef3peak.pdf", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/reso_vs_energy.pdf", outputdir.c_str() ) );


  return 0;

}







TF1* fitSingleElectronPeakUoDLoR( const std::string& outputdir, const std::string& name, TTree* tree, int niter, float nSigma, int isUp, bool isLeft ) {

  std::string histoName(Form("h1_%s", name.c_str()));
  TH1D* h1 = new TH1D(histoName.c_str(), "", 200, 0., 6000.);

  //////////////////////left right with -2 < and > 2:///////////////
  /*
  if(isUp==1){tree->Project( histoName.c_str(),"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 && pos_hodoClustY>0.)");
  } else if(isUp==0){
      tree->Project( histoName.c_str(),"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 && pos_hodoClustY<0. )");
  }  else if(isUp==-1){
  if(isLeft==1){tree->Project( histoName.c_str(),"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 && pos_hodoClustX<-2.)");
	       } else{
    tree->Project( histoName.c_str(),"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 && pos_hodoClustX>2. )");}}
  */
  
  if(isUp==1){tree->Project( histoName.c_str(),"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 && pos_hodoClustY>0.)");
  } else if(isUp==0){
      tree->Project( histoName.c_str(),"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 && pos_hodoClustY<0. )");
  }  else if(isUp==-1){
  if(isLeft==1){tree->Project( histoName.c_str(),"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 && pos_hodoClustX<0.)");
	       } else{
    tree->Project( histoName.c_str(),"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 && pos_hodoClustX>0. )");}}
  

  TF1* f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 0., 6000.);
  f1->SetParameter(0, h1->Integral() );
  f1->SetParameter(1, h1->GetMean() );
  f1->SetParameter(2, h1->GetRMS() );

  f1->SetParError(1, h1->GetMeanError() );
  f1->SetParError(2, h1->GetRMSError() );

  FitTools::doSingleFit( h1, f1, outputdir, name, niter, nSigma );

  return f1;

}



TH1* meanOfSingleElectronPeakBGO( const std::string& outputdir, const std::string& name, TTree* tree ) {

  std::string histoName(Form("h1_%s", name.c_str()));
  TH1D* h1 = new TH1D(histoName.c_str(), "", 200, 5., 3000.);

tree->Project( histoName.c_str(),"bgo_corr[0]+bgo_corr[1]+bgo_corr[2]+bgo_corr[3]+bgo_corr[4]+bgo_corr[5]+bgo_corr[6]+bgo_corr[7]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1)");

 TCanvas* can = new TCanvas("can","", 600,600);
 can->cd();

 h1->Draw();
 can->SaveAs( Form("%s/histo_%s.pdf", outputdir.c_str(), name.c_str() ) );

 delete can;


  return h1;

}



TF1* fitSingleElectronPeakBGO( const std::string& outputdir, const std::string& name, TTree* tree, int niter, float nSigma ) {

  std::string histoName(Form("h1_%s", name.c_str()));
  TH1D* h1 = new TH1D(histoName.c_str(), "", 200, 0., 6000.);

tree->Project( histoName.c_str(),"bgo_corr[0]+bgo_corr[1]+bgo_corr[2]+bgo_corr[3]+bgo_corr[4]+bgo_corr[5]+bgo_corr[6]+bgo_corr[7]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 )");
 

  TF1* f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 0., 6000.);
  f1->SetParameter(0, h1->Integral() );
  f1->SetParameter(1, h1->GetMean() );
  f1->SetParameter(2, h1->GetRMS() );

  f1->SetParError(1, h1->GetMeanError() );
  f1->SetParError(2, h1->GetRMSError() );

  FitTools::doSingleFit( h1, f1, outputdir, name, niter, nSigma );

  return f1;

}
