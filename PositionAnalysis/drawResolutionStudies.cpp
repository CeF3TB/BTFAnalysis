#include <cmath>

#include "TFile.h"
#include "TTree.h"


#include "interface/DrawTools.h"
#include "interface/FitTools.h"



struct ResoStruct {

  float resp;
  float resp_error;

  float reso;
  float reso_error;

  float Sres;
  float Sres_error;

};




ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, float LYSF[], const std::string& name );
ResoStruct getRespAndReso( TF1* f1, float energyErrorPercent );


int main() {

  DrawTools::setStyle();

  TFile* file_mc = TFile::Open("EEShash_491MeV_10000ev.root");

  TFile* file_data = TFile::Open("analysisTrees_V00/Reco_BTF_259_20140502-012847_beam.root");

  TTree* tree_data = (TTree*)file_data->Get("recoTree");
  TTree* tree_mc = (TTree*)file_mc->Get("EEShash");

  std::string outputdir = "ResolutionStudiesPlots";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );
  
  TF1* f1_data = FitTools::fitSingleElectronPeak( outputdir, "data", tree_data );

  ResoStruct rs_data = getRespAndReso(f1_data, 0.01);



  float LYSF_ideal[] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
  float LYSF_real [] = {0.85, 0.94, 0.95, 0.98, 1.00, 1.02, 1.05, 1.05, 1.05, 0.74};
  float LYSF_hole [] = {0.85, 0.94, 0.95, 0.98, 1.00, 0.00, 1.05, 1.05, 1.05, 0.74};

  ResoStruct rs_ideal = getResponseResolutionMC( outputdir, tree_mc, LYSF_ideal, "ideal" );
  ResoStruct rs_real  = getResponseResolutionMC( outputdir, tree_mc, LYSF_real, "real" );
  ResoStruct rs_hole  = getResponseResolutionMC( outputdir, tree_mc, LYSF_hole, "hole" );
  
  std::cout << std::endl;
  std::cout << "* DATA: " << std::endl;
  std::cout << "reso: " << rs_data.reso << " +/- " << rs_data.reso_error << std::endl;
  std::cout << "Sres: " << rs_data.Sres << " +/- " << rs_data.Sres_error << std::endl;
  
  std::cout << "MC ideal: " << rs_ideal.Sres << " +/- " << rs_ideal.Sres_error << std::endl;
  std::cout << "MC real: " << rs_real.Sres << " +/- " << rs_real.Sres_error << std::endl;
  std::cout << "MC hole: " << rs_hole.Sres << " +/- " << rs_hole.Sres_error << std::endl;


  return 0;

}




ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, float LYSF[], const std::string& name ) {

  std::string fullVarName = "";
  for( unsigned i=0; i<10; ++i ) {
    std::string plusSign = (i==0) ? "" : " + ";
    std::string thisPiece(Form("%s%f*Eact_%d", plusSign.c_str(), LYSF[i], i));
    fullVarName += thisPiece;
  }

  TH1D* h1 = new TH1D( name.c_str(), "", 500, 0., 500. );

  tree->Project( name.c_str(), fullVarName.c_str() );

  TF1* f1 = new TF1( Form("f1_%s", name.c_str()), "gaus", 100., 500.);

  FitTools::doSingleFit( h1, f1, outputdir, name );


  ResoStruct rs = getRespAndReso( f1, 0. );

  return rs;

}



ResoStruct getRespAndReso( TF1* f1, float energyErrorPercent ) {

  float energy = 491.;
  float energyErr = energyErrorPercent*energy;

  float mean = f1->GetParameter(1);
  float meanErr = f1->GetParError(1);

  float rms = f1->GetParameter(2);
  float rmsErr = f1->GetParError(2);

  float reso = 100.* rms/mean; //in percent
  float resoErr = 100.* sqrt( rmsErr*rmsErr/(mean*mean) + rms*rms*meanErr*meanErr/(mean*mean*mean*mean) );

  float energyGeV = energy/1000.;
  float energyGeVErr = energyErr/1000.;
  float Sres = reso*sqrt(energyGeV);
  float Sres_error = sqrt( resoErr*resoErr*energyGeV + reso*reso*0.5*0.5*energyGeVErr*energyGeVErr/energyGeV );

  ResoStruct rs;
  rs.resp = mean;
  rs.resp_error = meanErr;
  rs.reso = reso;
  rs.reso_error = resoErr;
  rs.Sres = Sres;
  rs.Sres_error = Sres_error;

  return rs;

}
