#include <iostream>
#include <fstream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"

#include "interface/DrawTools.h"

#include "fastDQM_CeF3_BTF.h"




struct HistoAndName {

  HistoAndName( TH1D* h, const std::string& n ) {
    histo = h;
    name = n;
  }

  TH1D* histo;
  std::string name;

};



struct Histo2AndName {

  Histo2AndName( TH2D* h, const std::string& n ) {
    histo = h;
    name = n;
  }

  TH2D* histo;
  std::string name;

};




TH2D* getSingleHisto( TTree* tree, const std::string& name, const std::string& varName, const std::string& varExpr, const std::string& cuts );
std::pair<TH1D*,TH1D*> getBiasAndResoHistos( TH2D* h2 );
void checkLateralScan( const std::string& outputdir, const std::string& name, TTree* tree, const std::string& refVar, const std::string& axisName, const std::string& cut );
void drawPerformancePlot( const std::string& outputdir, const std::string& name, const std::string& var, std::vector< HistoAndName > hn, const std::string& bias_reso, const std::string& axisName, float xMin, float xMax );
void drawProjections( const std::string& outputdir, const std::string& name, std::vector<Histo2AndName> hn, const std::string& var );
void checkBGOScan( const std::string& outputdir, const std::string& name, TTree* tree );


int main( int argc, char* argv[] ) {

 // if( argv==1 ) {
 //   std::cout << "USAGE: ./checkBGOPositioning [runSet] [tag=\"V02\"]" << std::endl;
 //   exit(11);
 // }



  std::string tag = "V02";
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }


  TChain* tree = new TChain("posTree");
  tree->Add( Form("PosAnTrees_%s/crossScanFile.root/posTree", tag.c_str()) );
  tree->Add( Form("PosAnTrees_%s/diagonalScanFile.root/posTree", tag.c_str()) );


  DrawTools::setStyle();

  std::string outputdir = "BGOPositioningPlots_" + tag;
  system(Form("mkdir -p %s", outputdir.c_str()));

  checkLateralScan( outputdir, "all"    , tree, "Beam", "Beam", "(yBeam==0. && abs(xBeam)>0.2) || (xBeam==0. && abs(yBeam)>0.2) || (xBeam==yBeam && xBeam>0.) || (xBeam==-yBeam && xBeam>0.)" );
  checkLateralScan( outputdir, "horiz"  , tree, "Beam", "Beam", "yBeam==0. && abs(xBeam)>0.2" );
  checkLateralScan( outputdir, "vert"   , tree, "Beam", "Beam", "xBeam==0. && abs(yBeam)>0.2" );
  checkLateralScan( outputdir, "diag13" , tree, "Beam", "Beam", "xBeam==yBeam && xBeam>0." );
  checkLateralScan( outputdir, "diag02" , tree, "Beam", "Beam", "xBeam==-yBeam && xBeam>0." );

  TFile* file_hodo = TFile::Open(Form("PosAnTrees_%s/PosAn_BTF_246_beam.root", tag.c_str()));
  TTree* tree_hodo = (TTree*)file_hodo->Get("posTree");
  checkLateralScan( outputdir, "hodo" , tree_hodo, "Pos_hodo", "Hodoscope", "nHodoClustersX==1 && nHodoClustersY==1 && nFibres_hodoClustX[0]<=2 && nFibres_hodoClustY[0]<=2" );

  //TFile* file_bgo = TFile::Open(Form("PosAnTrees_%s/BGORuns.root", tag.c_str()));
  //TTree* tree_bgo = (TTree*)file_bgo->Get("posTree");
  //checkBGOScan( outputdir, "bgo", tree_bgo );

  return 0;

}



void checkBGOScan( const std::string& outputdir, const std::string& name, TTree* tree ) {

  std::vector<int> runs;
  runs.push_back(229); //0
  runs.push_back(231); //1
  runs.push_back(233); //2
  runs.push_back(237); //3
  runs.push_back(235); //4
  runs.push_back(239); //5
  runs.push_back(241); //6
  runs.push_back(243); //7


  float xMin = -0.5;
  float xMax = (-0.5+(float)runs.size());

  TH1D* h1_bias_asymm_x = new TH1D( Form("%s_bias_asymm_x_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );
  TH1D* h1_bias_asymm_log_x = new TH1D( Form("%s_bias_asymm_log_x_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );
  TH1D* h1_bias_wa_x = new TH1D( Form("%s_bias_wa_x_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );
  TH1D* h1_bias_wa_log_x = new TH1D( Form("%s_bias_wa_log_x_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );

  TH1D* h1_reso_asymm_x = new TH1D( Form("%s_reso_asymm_x_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );
  TH1D* h1_reso_asymm_log_x = new TH1D( Form("%s_reso_asymm_log_x_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );
  TH1D* h1_reso_wa_x = new TH1D( Form("%s_reso_wa_x_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );
  TH1D* h1_reso_wa_log_x = new TH1D( Form("%s_reso_wa_log_x_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );

  TH1D* h1_bias_asymm_y = new TH1D( Form("%s_bias_asymm_y_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );
  TH1D* h1_bias_asymm_log_y = new TH1D( Form("%s_bias_asymm_log_y_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );
  TH1D* h1_bias_wa_y = new TH1D( Form("%s_bias_wa_y_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );
  TH1D* h1_bias_wa_log_y = new TH1D( Form("%s_bias_wa_log_y_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );

  TH1D* h1_reso_asymm_y = new TH1D( Form("%s_reso_asymm_y_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );
  TH1D* h1_reso_asymm_log_y = new TH1D( Form("%s_reso_asymm_log_y_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );
  TH1D* h1_reso_wa_y = new TH1D( Form("%s_reso_wa_y_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );
  TH1D* h1_reso_wa_log_y = new TH1D( Form("%s_reso_wa_log_y_vs_chan", name.c_str()), "", runs.size(), xMin, xMax );



  for( unsigned i=0; i<runs.size(); ++i ) {

    std::string cut(Form("run==%d", runs[i]));

    TH2D* h2_asymm_x = getSingleHisto( tree, name, "asymm_x", "(xPos_bgo_asymm-xBeam):xBeam", cut );
    TH2D* h2_asymm_y = getSingleHisto( tree, name, "asymm_y", "(yPos_bgo_asymm-yBeam):yBeam", cut );
    TH2D* h2_asymm_log_x = getSingleHisto( tree, name, "asymm_log_x", "(xPos_bgo_asymmlog-xBeam):xBeam", cut );
    TH2D* h2_asymm_log_y = getSingleHisto( tree, name, "asymm_log_y", "(yPos_bgo_asymmlog-yBeam):yBeam", cut );
    TH2D* h2_wa_x = getSingleHisto( tree, name, "wa_x", "(xPos_bgo_wa-xBeam):xBeam", cut );
    TH2D* h2_wa_y = getSingleHisto( tree, name, "wa_y", "(yPos_bgo_wa-yBeam):yBeam", cut );
    TH2D* h2_wa_log_x = getSingleHisto( tree, name, "wa_log_x", "(xPos_bgo_walog-xBeam):xBeam", cut );
    TH2D* h2_wa_log_y = getSingleHisto( tree, name, "wa_log_y", "(yPos_bgo_walog-yBeam):yBeam", cut );

    std::vector< Histo2AndName > hn_x;
    hn_x.push_back( Histo2AndName(h2_asymm_x, "Asymmetry") );
    hn_x.push_back( Histo2AndName(h2_asymm_log_x, "Log-Asymmetry") );
    hn_x.push_back( Histo2AndName(h2_wa_x, "Weighted Average") );
    hn_x.push_back( Histo2AndName(h2_wa_log_x, "Log-Weighted Average") );

    std::vector< Histo2AndName > hn_y;
    hn_y.push_back( Histo2AndName(h2_asymm_y, "Asymmetry") );
    hn_y.push_back( Histo2AndName(h2_asymm_log_y, "Log-Asymmetry") );
    hn_y.push_back( Histo2AndName(h2_wa_y, "Weighted Average") );
    hn_y.push_back( Histo2AndName(h2_wa_log_y, "Log-Weighted Average") );

    drawProjections( outputdir, Form("%s_%d", name.c_str(), i), hn_x, "x" );
    drawProjections( outputdir, Form("%s_%d", name.c_str(), i), hn_y, "y" );



    std::pair<TH1D*,TH1D*> p_asymm_x = getBiasAndResoHistos(h2_asymm_x);
    std::pair<TH1D*,TH1D*> p_asymm_y = getBiasAndResoHistos(h2_asymm_y);
    std::pair<TH1D*,TH1D*> p_asymm_log_x = getBiasAndResoHistos(h2_asymm_log_x);
    std::pair<TH1D*,TH1D*> p_asymm_log_y = getBiasAndResoHistos(h2_asymm_log_y);
    std::pair<TH1D*,TH1D*> p_wa_x = getBiasAndResoHistos(h2_wa_x);
    std::pair<TH1D*,TH1D*> p_wa_y = getBiasAndResoHistos(h2_wa_y);
    std::pair<TH1D*,TH1D*> p_wa_log_x = getBiasAndResoHistos(h2_wa_log_x);
    std::pair<TH1D*,TH1D*> p_wa_log_y = getBiasAndResoHistos(h2_wa_log_y);

    TH1D* refHisto_x = p_asymm_x.first;
    TH1D* refHisto_y = p_asymm_y.first;


    for( unsigned ibin=1; ibin<refHisto_x->GetNbinsX()+1; ibin++ ) {

      if( refHisto_x->GetBinError(ibin)>0. ) {
 
        h1_bias_asymm_x->SetBinContent(i+1, p_asymm_x.first->GetBinContent(ibin));
        h1_bias_asymm_log_x->SetBinContent(i+1, p_asymm_log_x.first->GetBinContent(ibin));
        h1_bias_wa_x->SetBinContent(i+1, p_wa_x.first->GetBinContent(ibin));
        h1_bias_wa_log_x->SetBinContent(i+1, p_wa_log_x.first->GetBinContent(ibin));

        h1_bias_asymm_x->SetBinError(i+1, p_asymm_x.first->GetBinError(ibin));
        h1_bias_asymm_log_x->SetBinError(i+1, p_asymm_log_x.first->GetBinError(ibin));
        h1_bias_wa_x->SetBinError(i+1, p_wa_x.first->GetBinError(ibin));
        h1_bias_wa_log_x->SetBinError(i+1, p_wa_log_x.first->GetBinError(ibin));

        h1_reso_asymm_x->SetBinContent(i+1, p_asymm_x.second->GetBinContent(ibin));
        h1_reso_asymm_log_x->SetBinContent(i+1, p_asymm_log_x.second->GetBinContent(ibin));
        h1_reso_wa_x->SetBinContent(i+1, p_wa_x.second->GetBinContent(ibin));
        h1_reso_wa_log_x->SetBinContent(i+1, p_wa_log_x.second->GetBinContent(ibin));

        h1_reso_asymm_x->SetBinError(i+1, p_asymm_x.second->GetBinError(ibin));
        h1_reso_asymm_log_x->SetBinError(i+1, p_asymm_log_x.second->GetBinError(ibin));
        h1_reso_wa_x->SetBinError(i+1, p_wa_x.second->GetBinError(ibin));
        h1_reso_wa_log_x->SetBinError(i+1, p_wa_log_x.second->GetBinError(ibin));
       
      }


      if( refHisto_y->GetBinError(ibin)>0. ) {
 
        h1_bias_asymm_y->SetBinContent(i+1, p_asymm_y.first->GetBinContent(ibin));
        h1_bias_asymm_log_y->SetBinContent(i+1, p_asymm_log_y.first->GetBinContent(ibin));
        h1_bias_wa_y->SetBinContent(i+1, p_wa_y.first->GetBinContent(ibin));
        h1_bias_wa_log_y->SetBinContent(i+1, p_wa_log_y.first->GetBinContent(ibin));

        h1_bias_asymm_y->SetBinError(i+1, p_asymm_y.first->GetBinError(ibin));
        h1_bias_asymm_log_y->SetBinError(i+1, p_asymm_log_y.first->GetBinError(ibin));
        h1_bias_wa_y->SetBinError(i+1, p_wa_y.first->GetBinError(ibin));
        h1_bias_wa_log_y->SetBinError(i+1, p_wa_log_y.first->GetBinError(ibin));

        h1_reso_asymm_y->SetBinContent(i+1, p_asymm_y.second->GetBinContent(ibin));
        h1_reso_asymm_log_y->SetBinContent(i+1, p_asymm_log_y.second->GetBinContent(ibin));
        h1_reso_wa_y->SetBinContent(i+1, p_wa_y.second->GetBinContent(ibin));
        h1_reso_wa_log_y->SetBinContent(i+1, p_wa_log_y.second->GetBinContent(ibin));

        h1_reso_asymm_y->SetBinError(i+1, p_asymm_y.second->GetBinError(ibin));
        h1_reso_asymm_log_y->SetBinError(i+1, p_asymm_log_y.second->GetBinError(ibin));
        h1_reso_wa_y->SetBinError(i+1, p_wa_y.second->GetBinError(ibin));
        h1_reso_wa_log_y->SetBinError(i+1, p_wa_log_y.second->GetBinError(ibin));

      }

    }  // for ibin

  }  // for i run

        
  std::vector< HistoAndName > hn_bias_x;
  hn_bias_x.push_back( HistoAndName(h1_bias_asymm_x, "Asymmetry") );
  hn_bias_x.push_back( HistoAndName(h1_bias_asymm_log_x, "Log-Asymmetry") );
  hn_bias_x.push_back( HistoAndName(h1_bias_wa_x, "Weighted Average") );
  hn_bias_x.push_back( HistoAndName(h1_bias_wa_log_x, "Log-Weighted Average") );

  std::vector< HistoAndName > hn_reso_x;
  hn_reso_x.push_back( HistoAndName(h1_reso_asymm_x, "Asymmetry") );
  hn_reso_x.push_back( HistoAndName(h1_reso_asymm_log_x, "Log-Asymmetry") );
  hn_reso_x.push_back( HistoAndName(h1_reso_wa_x, "Weighted Average") );
  hn_reso_x.push_back( HistoAndName(h1_reso_wa_log_x, "Log-Weighted Average") );
        
  std::vector< HistoAndName > hn_bias_y;
  hn_bias_y.push_back( HistoAndName(h1_bias_asymm_y, "Asymmetry") );
  hn_bias_y.push_back( HistoAndName(h1_bias_asymm_log_y, "Log-Asymmetry") );
  hn_bias_y.push_back( HistoAndName(h1_bias_wa_y, "Weighted Average") );
  hn_bias_y.push_back( HistoAndName(h1_bias_wa_log_y, "Log-Weighted Average") );

  std::vector< HistoAndName > hn_reso_y;
  hn_reso_y.push_back( HistoAndName(h1_reso_asymm_y, "Asymmetry") );
  hn_reso_y.push_back( HistoAndName(h1_reso_asymm_log_y, "Log-Asymmetry") );
  hn_reso_y.push_back( HistoAndName(h1_reso_wa_y, "Weighted Average") );
  hn_reso_y.push_back( HistoAndName(h1_reso_wa_log_y, "Log-Weighted Average") );

   
  drawPerformancePlot( outputdir, name, "x", hn_bias_x, "bias", "BGO Channel", xMin, xMax );
  drawPerformancePlot( outputdir, name, "x", hn_reso_x, "reso", "BGO Channel", xMin, xMax );
  drawPerformancePlot( outputdir, name, "y", hn_bias_y, "bias", "BGO Channel", xMin, xMax );
  drawPerformancePlot( outputdir, name, "y", hn_reso_y, "reso", "BGO Channel", xMin, xMax );

    

}




void checkLateralScan( const std::string& outputdir, const std::string& name, TTree* tree, const std::string& refVar, const std::string& axisName, const std::string& cut ) {

  TH2D* h2_asymm_x = getSingleHisto( tree, name, "asymm_x", Form("(xPos_bgo_asymm-x%s):x%s", refVar.c_str(), refVar.c_str()), cut );
  TH2D* h2_asymm_y = getSingleHisto( tree, name, "asymm_y", Form("(yPos_bgo_asymm-y%s):y%s", refVar.c_str(), refVar.c_str()), cut );
  TH2D* h2_asymm_log_x = getSingleHisto( tree, name, "asymm_log_x", Form("(xPos_bgo_asymmlog-x%s):x%s", refVar.c_str(), refVar.c_str()), cut );
  TH2D* h2_asymm_log_y = getSingleHisto( tree, name, "asymm_log_y", Form("(yPos_bgo_asymmlog-y%s):y%s", refVar.c_str(), refVar.c_str()), cut );
  TH2D* h2_wa_x = getSingleHisto( tree, name, "wa_x", Form("(xPos_bgo_wa-x%s):x%s", refVar.c_str(), refVar.c_str()), cut );
  TH2D* h2_wa_y = getSingleHisto( tree, name, "wa_y", Form("(yPos_bgo_wa-y%s):y%s", refVar.c_str(), refVar.c_str()), cut );
  TH2D* h2_wa_log_x = getSingleHisto( tree, name, "wa_log_x", Form("(xPos_bgo_walog-x%s):x%s", refVar.c_str(), refVar.c_str()), cut );
  TH2D* h2_wa_log_y = getSingleHisto( tree, name, "wa_log_y", Form("(yPos_bgo_walog-y%s):y%s", refVar.c_str(), refVar.c_str()), cut );


  std::vector< Histo2AndName > hn_x;
  hn_x.push_back( Histo2AndName(h2_asymm_x, "Asymmetry") );
  hn_x.push_back( Histo2AndName(h2_asymm_log_x, "Log-Asymmetry") );
  hn_x.push_back( Histo2AndName(h2_wa_x, "Weighted Average") );
  hn_x.push_back( Histo2AndName(h2_wa_log_x, "Log-Weighted Average") );

  std::vector< Histo2AndName > hn_y;
  hn_y.push_back( Histo2AndName(h2_asymm_y, "Asymmetry") );
  hn_y.push_back( Histo2AndName(h2_asymm_log_y, "Log-Asymmetry") );
  hn_y.push_back( Histo2AndName(h2_wa_y, "Weighted Average") );
  hn_y.push_back( Histo2AndName(h2_wa_log_y, "Log-Weighted Average") );

  drawProjections( outputdir, name, hn_x, "x" );
  drawProjections( outputdir, name, hn_y, "y" );

  std::pair<TH1D*,TH1D*> p_asymm_x = getBiasAndResoHistos(h2_asymm_x);
  std::pair<TH1D*,TH1D*> p_asymm_y = getBiasAndResoHistos(h2_asymm_y);
  std::pair<TH1D*,TH1D*> p_asymm_log_x = getBiasAndResoHistos(h2_asymm_log_x);
  std::pair<TH1D*,TH1D*> p_asymm_log_y = getBiasAndResoHistos(h2_asymm_log_y);
  std::pair<TH1D*,TH1D*> p_wa_x = getBiasAndResoHistos(h2_wa_x);
  std::pair<TH1D*,TH1D*> p_wa_y = getBiasAndResoHistos(h2_wa_y);
  std::pair<TH1D*,TH1D*> p_wa_log_x = getBiasAndResoHistos(h2_wa_log_x);
  std::pair<TH1D*,TH1D*> p_wa_log_y = getBiasAndResoHistos(h2_wa_log_y);

  std::vector< HistoAndName > hn_bias_x;
  hn_bias_x.push_back( HistoAndName(p_asymm_x.first, "Asymmetry") );
  hn_bias_x.push_back( HistoAndName(p_asymm_log_x.first, "Log-Asymmetry") );
  hn_bias_x.push_back( HistoAndName(p_wa_x.first, "Weighted Average") );
  hn_bias_x.push_back( HistoAndName(p_wa_log_x.first, "Log-Weighted Average") );

  std::vector< HistoAndName > hn_bias_y;
  hn_bias_y.push_back( HistoAndName(p_asymm_y.first, "Asymmetry") );
  hn_bias_y.push_back( HistoAndName(p_asymm_log_y.first, "Log-Asymmetry") );
  hn_bias_y.push_back( HistoAndName(p_wa_y.first, "Weighted Average") );
  hn_bias_y.push_back( HistoAndName(p_wa_log_y.first, "Log-Weighted Average") );

  std::vector< HistoAndName > hn_reso_x;
  hn_reso_x.push_back( HistoAndName(p_asymm_x.second, "Asymmetry") );
  hn_reso_x.push_back( HistoAndName(p_asymm_log_x.second, "Log-Asymmetry") );
  hn_reso_x.push_back( HistoAndName(p_wa_x.second, "Weighted Average") );
  hn_reso_x.push_back( HistoAndName(p_wa_log_x.second, "Log-Weighted Average") );

  std::vector< HistoAndName > hn_reso_y;
  hn_reso_y.push_back( HistoAndName(p_asymm_y.second, "Asymmetry") );
  hn_reso_y.push_back( HistoAndName(p_asymm_log_y.second, "Log-Asymmetry") );
  hn_reso_y.push_back( HistoAndName(p_wa_y.second, "Weighted Average") );
  hn_reso_y.push_back( HistoAndName(p_wa_log_y.second, "Log-Weighted Average") );

  if( name!="vert" ) drawPerformancePlot( outputdir, name, "x", hn_bias_x, "bias", Form("%s x Position [mm]", axisName.c_str()), -10., 10. );
  if( name!="horiz") drawPerformancePlot( outputdir, name, "y", hn_bias_y, "bias", Form("%s x Position [mm]", axisName.c_str()), -10., 10. );
  if( name!="vert" ) drawPerformancePlot( outputdir, name, "x", hn_reso_x, "reso", Form("%s y Position [mm]", axisName.c_str()), -10., 10. );
  if( name!="horiz") drawPerformancePlot( outputdir, name, "y", hn_reso_y, "reso", Form("%s y Position [mm]", axisName.c_str()), -10., 10. );

}



std::pair<TH1D*,TH1D*> getBiasAndResoHistos( TH2D* h2 ) {

  int  nBins = h2->GetXaxis()->GetNbins();
  float xMin = h2->GetXaxis()->GetXmin();
  float xMax = h2->GetXaxis()->GetXmax();

  TH1D* h1_bias = new TH1D( Form("bias_%s", h2->GetName()), "", nBins, xMin, xMax );
  TH1D* h1_reso = new TH1D( Form("reso_%s", h2->GetName()), "", nBins, xMin, xMax );

  for( unsigned i=0; i<nBins; ++i ) {

    TH1D* h1_proj = h2->ProjectionY(Form("%s_%d", h2->GetName(), i), i+1, i+1);

    if( h1_proj->GetEntries()>0 ) {
      h1_bias->SetBinContent( i+1, h1_proj->GetMean() );
      h1_bias->SetBinError  ( i+1, h1_proj->GetMeanError() );
      h1_reso->SetBinContent( i+1, h1_proj->GetRMS() );
      h1_reso->SetBinError  ( i+1, h1_proj->GetRMSError() );
    }

  }

  std::pair<TH1D*,TH1D*> returnPair;
  returnPair.first  = h1_bias;
  returnPair.second = h1_reso;

  return returnPair;

}




TH2D* getSingleHisto( TTree* tree, const std::string& name, const std::string& varName, const std::string& varExpr, const std::string& cuts ) {

  std::string histoName = name + "_" + varName;
  TH2D* h2 = new TH2D( histoName.c_str(), "", 81, -40.5, 40.5, 40, -40., 40.);

  std::string fullCuts = cuts + " && isSingleEle_scintFront";
  tree->Project( histoName.c_str(), varExpr.c_str(), fullCuts.c_str() );

  return h2;

}



void drawProjections( const std::string& outputdir, const std::string& name, std::vector<Histo2AndName> hn, const std::string& var ) {

  std::vector<int> colors;
  colors.push_back( 38 );
  colors.push_back( 46 );
  colors.push_back( 29 );
  colors.push_back( 42 );
  colors.push_back( 40 );
  colors.push_back( 41 );
  colors.push_back( 45 );
  colors.push_back( 30 );
  colors.push_back( 16 );  


  int  nBins = hn[0].histo->GetXaxis()->GetNbins();
  float xMin = hn[0].histo->GetXaxis()->GetXmin();
  float xMax = hn[0].histo->GetXaxis()->GetXmax();
  float binSize = (xMax-xMin)/((float)nBins);

  for( unsigned i=0; i<nBins; ++i ) {

    TCanvas* c1 = new TCanvas("cxx", "", 600, 600);
    c1->cd();

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 0.2);
    h2_axes->SetXTitle("Estimate - Beam [mm]");
    h2_axes->SetYTitle("Normalized to Unity");
    h2_axes->Draw();

    float xMin_i = xMin + (float)i*binSize;

    TLegend* legend = new TLegend( 0.2, 0.6, 0.45, 0.9, Form("%.1f < %s < %.1f mm", xMin_i, var.c_str(), xMin_i+binSize) );
    legend->SetTextSize(0.033);
    legend->SetFillColor(0);

    bool saveIt = true;

    for( unsigned j=0; j<hn.size(); j++ ) { 

      TH1D* h1_proj = hn[j].histo->ProjectionY(Form("%s_%d", hn[j].histo->GetName(), i), i+1, i+1);
  
      if( h1_proj->GetEntries()==0 ) 
        saveIt = false;

      if( saveIt ) {
        h1_proj->SetLineColor( colors[j] );
        h1_proj->SetLineWidth(2);
        h1_proj->DrawNormalized("same");
        legend->AddEntry( h1_proj, hn[j].name.c_str(), "L" );
      }
      
    }

    if( saveIt ) {

      legend->Draw("same");

      TPaveText* labelTop = DrawTools::getLabelTop();
      labelTop->Draw("same");

      gPad->RedrawAxis();

      std::string subdir = "proj_" + name;
      system( Form("mkdir -p %s/%s", outputdir.c_str(), subdir.c_str()) );


      c1->SaveAs(Form("%s/%s/proj%s_%s_%d.eps", outputdir.c_str(), subdir.c_str(), var.c_str(), name.c_str(), i) );
      c1->SaveAs(Form("%s/%s/proj%s_%s_%d.png", outputdir.c_str(), subdir.c_str(), var.c_str(), name.c_str(), i) );
      c1->SaveAs(Form("%s/%s/proj%s_%s_%d.pdf", outputdir.c_str(), subdir.c_str(), var.c_str(), name.c_str(), i) );

    }

    delete c1;
    delete h2_axes;

  } // i

}


void drawPerformancePlot( const std::string& outputdir, const std::string& name, const std::string& var, std::vector< HistoAndName > hn, const std::string& bias_reso, const std::string& axisName, float xMin, float xMax ) {

  bool isReso;
  if( bias_reso=="bias" ) isReso = false;
  else if( bias_reso=="reso" ) isReso = true;
  else {
    std::cout << "ERROR! bias_reso has to be either 'bias' or 'reso'. Exiting." << std::endl;
    exit(91);
  }


  std::vector<int> colors;
  colors.push_back( 38 );
  colors.push_back( 46 );
  colors.push_back( 29 );
  colors.push_back( 42 );
  colors.push_back( 40 );
  colors.push_back( 41 );
  colors.push_back( 45 );
  colors.push_back( 30 );
  colors.push_back( 16 );  




  TCanvas* c1 = new TCanvas( "cx", "", 600, 600 );
  c1->Clear();
  c1->cd();


  // first: bias
  //float xMax = (name=="bgo") ? 40. : 15.;
  float yMin = (isReso) ? 0.  : -1.;
  float yMax = (isReso) ? 15. : 3.;
  TString name_tstr(name);
  if( name_tstr.Contains("bgo") ) {
    if(!isReso )
      yMin = -12.;
    yMax = 12.;
  }

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, yMin, yMax );
  h2_axes->SetXTitle( axisName.c_str() );
  //h2_axes->SetXTitle( Form( "Beam %s Position [mm]", var.c_str() ) );
  if(isReso)
    h2_axes->SetYTitle( "RMS(Estimate - Beam) [mm]" );
  else
    h2_axes->SetYTitle( "Estimate - Beam [mm]" );

  

  TLegend* legend = new TLegend( 0.2, 0.69, 0.4, 0.92 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.033);

  for( unsigned i=0; i<hn.size(); ++i ) {

    hn[i].histo->SetLineColor( colors[i] );
    hn[i].histo->SetMarkerColor( colors[i] );
    hn[i].histo->SetMarkerSize( 1.6 );
    hn[i].histo->SetMarkerStyle( 20+i );

    TF1* f1 = new TF1( Form("f1_%s_%s_%d", name.c_str(), var.c_str(), i), "[0]", -xMax, xMax);
    f1->SetRange(-xMax, xMax);
    f1->SetLineWidth(2);
    f1->SetLineColor(colors[i]);
    hn[i].histo->Fit(f1, "R+");

    legend->AddEntry( hn[i].histo, hn[i].name.c_str(), "P" );
  
  }

  h2_axes->Draw();

  for( unsigned i=0; i<hn.size(); ++i )
    hn[i].histo->Draw("p same");

  legend->Draw("same");

  TLine* lineZero = new TLine(xMin, 0., xMax, 0.);
  lineZero->Draw("same");

  TPaveText* labelTop = DrawTools::getLabelTop();
  labelTop->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/%s_%s_%s.eps", outputdir.c_str(), name.c_str(), bias_reso.c_str(), var.c_str()) );
  c1->SaveAs( Form("%s/%s_%s_%s.png", outputdir.c_str(), name.c_str(), bias_reso.c_str(), var.c_str()) );
  c1->SaveAs( Form("%s/%s_%s_%s.pdf", outputdir.c_str(), name.c_str(), bias_reso.c_str(), var.c_str()) );

  delete c1;
  delete h2_axes;


}


