#include "../interface/FitTools.h"

#include "TCanvas.h"





TF1* FitTools::fitSingleElectronPeak( const std::string& outputdir, const std::string& name, TTree* tree, int niter, float nSigma ) {

  std::string histoName(Form("h1_%s", name.c_str()));
  TH1D* h1 = new TH1D(histoName.c_str(), "", 200, 0., 6000.);
  tree->Project( histoName.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1)");

  // TF1* f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 400., 6000.);
  TF1* f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 0., 6000.);
  f1->SetParameter(0, h1->Integral() );
  f1->SetParameter(1, h1->GetMean() );
  f1->SetParameter(2, h1->GetRMS() );

  f1->SetParError(1, h1->GetMeanError() );
  f1->SetParError(2, h1->GetRMSError() );

  doSingleFit( h1, f1, outputdir, name, niter, nSigma );

  return f1;

}




void FitTools::doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name, int niter, float nSigma ) {

  h1->Fit( f1, "RQN" );

  for( unsigned iter=0; iter<niter; iter++ ) {

    float mean  = f1->GetParameter(1);
    float sigma = f1->GetParameter(2);
    float fitMin = mean - nSigma*sigma;
    float fitMax = mean + nSigma*sigma;
    f1->SetRange( fitMin, fitMax );
    if( iter==(niter-1) )
      h1->Fit( f1, "RQN" );
    else
      h1->Fit( f1, "RQ+" );
  }


  TCanvas* c1 = new TCanvas("cX", "", 600, 600);
  c1->cd();

  h1->Draw();

  c1->SaveAs( Form("%s/fit_%s.eps", outputdir.c_str(), name.c_str()) );
  c1->SaveAs( Form("%s/fit_%s.png", outputdir.c_str(), name.c_str()) );

  delete c1;

}

