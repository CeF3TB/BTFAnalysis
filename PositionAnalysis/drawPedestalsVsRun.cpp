#include <cstdlib>
#include "interface/DrawTools.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TH2D.h"
#include "TF1.h"

#include "fastDQM_CeF3_BTF.h"



Double_t positiveConst(Double_t *x, Double_t *par) {
  if (x[0] < 0. ) {
     TF1::RejectPoint();
     return 0;
  }
  return par[0];
}


void drawPedestals( TFile* file, const std::string& name );


int main() {

  DrawTools::setStyle();

  TFile* file = TFile::Open("pedestalFile.root");

  std::system( "mkdir -p pedestalPlots" );

  drawPedestals( file, "cef3" );
  drawPedestals( file, "bgo" );
  drawPedestals( file, "hodox" );
  drawPedestals( file, "hodoy" );

  return 0;

}




void drawPedestals( TFile* file, const std::string& name ) {

  std::vector<int> colors;
  colors.push_back( kRed );
  colors.push_back( 38 );
  colors.push_back( 30 );
  colors.push_back( 42 );
  colors.push_back( 29 );
  colors.push_back( kBlack );
  colors.push_back( kGreen );
  colors.push_back( kBlue  );

  int nChannels;
  if( name=="cef3" )
    nChannels=CEF3_CHANNELS;
  else if( name=="bgo" )
    nChannels=BGO_CHANNELS;
  else if( name=="hodox" )
    nChannels=HODOX_CHANNELS;
  else if( name=="hodoy" )
    nChannels=HODOY_CHANNELS;
  else {
    std::cout << "Unknown type: " << name << std::endl;
    std::cout << "Exiting." << std::endl;
    return;
  }

  
  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  float typicalPed = (name=="bgo") ? 140. : 100.; 
  float yMax = ((float)nChannels+1.)*typicalPed*1.2;

  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 950., 10, 0., yMax );
  h2_axes->SetXTitle( "Run Number");
  h2_axes->SetYTitle( "Pedestal Mean [ADC Counts]" );

  h2_axes->Draw();

  TLegend* legend = new TLegend( 0.68, 0.5-0.025*nChannels, 0.9, 0.5+0.025*nChannels );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);

  std::vector<TH1D*> histos;


  float scalingFactor = (float)(nChannels);

  for( unsigned i=0; i<nChannels; ++i ) {

    TH1D* h1 = (TH1D*)file->Get( Form("%s_%d", name.c_str(), i) );

    h1->Scale( scalingFactor );
    h1->SetMarkerStyle( 20+i );
    h1->SetMarkerColor( colors[i] );
    h1->SetLineColor( colors[i] );

    TF1* f1 = new TF1("f1", positiveConst, 90., 600., 1 );
    f1->SetLineColor(colors[i]);
    h1->Fit(f1, "RQ");

    //h1->Draw("p same");

    legend->AddEntry( h1, Form("%s %d (#times %d)", name.c_str(), i, (int)scalingFactor), "P" );

    histos.push_back(h1);

    scalingFactor -= 1.;

    delete f1;

  }

  h2_axes->Draw();
  for(unsigned i=0; i<histos.size(); ++i ) histos[i]->Draw("p same");

  legend->Draw("same");

  gPad->RedrawAxis();
 
  c1->SaveAs( Form("pedestalPlots/%s.eps", name.c_str()) );
  c1->SaveAs( Form("pedestalPlots/%s.png", name.c_str()) );
  c1->SaveAs( Form("pedestalPlots/%s.pdf", name.c_str()) );

  delete c1;
  delete h2_axes;
  delete legend;

}
