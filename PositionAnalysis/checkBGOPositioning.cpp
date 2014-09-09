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
std::pair<TH1D*,TH1D*> getBiasAndResoHistosHodo( const std::string& outputdir, TH2D* h2 );
//void checkLateralScan( const std::string& outputdir, const std::string& name, TTree* tree, const std::string& refVar, const std::string& axisName, const std::string& cut );
void checkLateralScan( const std::string& outputdir, const std::string& name, TTree* tree, std::vector< std::pair< std::string, std::string> > names, const std::string& refVar, const std::string& axisName, const std::string& cut );
void drawPerformancePlot( const std::string& outputdir, const std::string& name, const std::string& var, std::vector< HistoAndName > hn, const std::string& bias_reso, const std::string& axisName, float xMin, float xMax );
void drawProjections( const std::string& outputdir, const std::string& name, std::vector<Histo2AndName> hn, const std::string& var );
void checkBGOScan( const std::string& outputdir, const std::string& name, TTree* tree );
void checkHodo( const std::string& outputdir, const std::string& name, TTree* tree, const std::string& var, const std::string& cut );
TF1* fitHodoWithBeam( const std::string& outputdir, const std::string& suffix, TH1D* h1, float bias, float r );


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


  std::vector< std::pair< std::string, std::string> > names; // first is varName, second is legendName
  names.push_back( std::pair<std::string, std::string>("bgo_asymm", "Asymmetry") );
  names.push_back( std::pair<std::string, std::string>("bgo_asymmlog", "Log-Asymmetry") );
  names.push_back( std::pair<std::string, std::string>("bgo_wa", "Weighted Average") );
  names.push_back( std::pair<std::string, std::string>("bgo_walog", "Log-Weighted Average") );

  checkLateralScan( outputdir, "all"    , tree, names, "Beam", "Beam", "(yBeam==0. && abs(xBeam)>0.2) || (xBeam==0. && abs(yBeam)>0.2) || (xBeam==yBeam && xBeam>0.) || (xBeam==-yBeam && xBeam>0.)" );
  checkLateralScan( outputdir, "horiz"  , tree, names, "Beam", "Beam", "yBeam==0. && abs(xBeam)>0.2" );
  checkLateralScan( outputdir, "vert"   , tree, names, "Beam", "Beam", "xBeam==0. && abs(yBeam)>0.2" );
  checkLateralScan( outputdir, "diag13" , tree, names, "Beam", "Beam", "xBeam==yBeam && xBeam>0." );
  checkLateralScan( outputdir, "diag02" , tree, names, "Beam", "Beam", "xBeam==-yBeam && xBeam>0." );

  TFile* file_hodo = TFile::Open(Form("PosAnTrees_%s/PosAn_BTF_246_beam.root", tag.c_str()));
  TTree* tree_hodo = (TTree*)file_hodo->Get("posTree");
  checkLateralScan( outputdir, "hodo" , tree_hodo, names, "Pos_hodo", "Hodoscope", "nHodoClustersX==1 && nHodoClustersY==1 && nFibres_hodoClustX[0]<=2 && nFibres_hodoClustY[0]<=2" );

  checkHodo( outputdir, "horizHodo", tree, "x", "run==149 || run==150 || run==151" );
  checkHodo( outputdir, "vertHodo" , tree, "y", "run==160 || run==161 || run==162" );

  checkHodo( outputdir, "diag02Hodo" , tree, "x", "run==176 || run==177 || run==94" );
  checkHodo( outputdir, "diag02Hodo" , tree, "y", "run==176 || run==177 || run==94" );

  checkHodo( outputdir, "diag13Hodo" , tree, "x", "run==136 || run==141" );
  checkHodo( outputdir, "diag13Hodo" , tree, "y", "run==136 || run==141" );


  TFile* file_hodo_defoc = TFile::Open(Form("PosAnTrees_%s/PosAn_BTF_248_beam.root", tag.c_str()));
  TTree* tree_hodo_defoc = (TTree*)file_hodo_defoc->Get("posTree");

  checkHodo( outputdir, "hodoDefoc" , tree_hodo_defoc, "x", "run>0" );
  checkHodo( outputdir, "hodoDefoc" , tree_hodo_defoc, "y", "run>0" );

  //checkLateralScan( outputdir, "hodoDefoc" , tree_hodo_defoc, "Pos_hodo", "Hodoscope", "nHodoClustersX==1 && nHodoClustersY==1 && nFibres_hodoClustX[0]<=2 && nFibres_hodoClustY[0]<=2" );


  //TFile* file_bgo = TFile::Open(Form("PosAnTrees_%s/BGORuns.root", tag.c_str()));
  //TTree* tree_bgo = (TTree*)file_bgo->Get("posTree");
  //checkBGOScan( outputdir, "bgo", tree_bgo );

  std::vector< std::pair< std::string, std::string> > names_calo; // first is varName, second is legendName
  names_calo.push_back( std::pair<std::string, std::string>("bgo_wa", "BGO Weighted Average") );
  names_calo.push_back( std::pair<std::string, std::string>("calo_wa", "3x3 Weighted Average") );
  names_calo.push_back( std::pair<std::string, std::string>("calo_walog", "3x3 Log-Weighted Average") );

  checkLateralScan( outputdir, "horizCalo"  , tree, names_calo, "Beam", "Beam", "yBeam==0. && abs(xBeam)>0.2 " );
  checkLateralScan( outputdir, "vertCalo"   , tree, names_calo, "Beam", "Beam", "xBeam==0. && abs(yBeam)>0.2 " );
  checkLateralScan( outputdir, "diag13Calo" , tree, names_calo, "Beam", "Beam", "xBeam==yBeam  && xBeam>0. " );
  checkLateralScan( outputdir, "diag02Calo" , tree, names_calo, "Beam", "Beam", "xBeam==-yBeam && xBeam>0. " );

  return 0;

}



void checkHodo( const std::string& outputdir, const std::string& name, TTree* tree, const std::string& var, const std::string& cut ) {

  std::string fullCuts = cut + " && isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 && nFibres_hodoClustX[0]<=2 && nFibres_hodoClustY[0]<=2";
    
  TH2D* h2_hodo = getSingleHisto( tree, name, Form("hodo_%s", var.c_str()), Form("(%sPos_hodo-%sBeam):%sBeam", var.c_str(), var.c_str(), var.c_str()), fullCuts );

  std::pair<TH1D*,TH1D*> p_hodo = getBiasAndResoHistosHodo(outputdir, h2_hodo);

  std::vector< HistoAndName > hn_bias_hodo;
  hn_bias_hodo.push_back( HistoAndName(p_hodo.first, "Hodoscope") );

  std::vector< HistoAndName > hn_reso_hodo;
  hn_reso_hodo.push_back( HistoAndName(p_hodo.second, "Hodoscope") );
   
  drawPerformancePlot( outputdir, name, var, hn_bias_hodo, "bias", Form("Beam %s Position [mm]", var.c_str()), -4., 4. );
  drawPerformancePlot( outputdir, name, var, hn_reso_hodo, "reso", Form("Beam %s Position [mm]", var.c_str()), -4., 4. );

}





void checkLateralScan( const std::string& outputdir, const std::string& name, TTree* tree, std::vector< std::pair< std::string, std::string> > names, const std::string& refVar, const std::string& axisName, const std::string& cut ) {


  std::vector< Histo2AndName > hn_x;
  std::vector< Histo2AndName > hn_y;

  std::vector< HistoAndName > hn_bias_x;
  std::vector< HistoAndName > hn_bias_y;

  std::vector< HistoAndName > hn_reso_x;
  std::vector< HistoAndName > hn_reso_y;



  for( unsigned i=0; i<names.size(); ++i ) {

    TH2D* h2_x = getSingleHisto( tree, name, Form("%s_x", names[i].first.c_str()), Form("(xPos_%s-x%s):x%s", names[i].first.c_str(), refVar.c_str(), refVar.c_str()), cut );
    TH2D* h2_y = getSingleHisto( tree, name, Form("%s_y", names[i].first.c_str()), Form("(yPos_%s-y%s):y%s", names[i].first.c_str(), refVar.c_str(), refVar.c_str()), cut );

    hn_x.push_back( Histo2AndName(h2_x, names[i].second) );
    hn_y.push_back( Histo2AndName(h2_y, names[i].second) );

    std::pair<TH1D*,TH1D*> p_x = getBiasAndResoHistos(h2_x);
    std::pair<TH1D*,TH1D*> p_y = getBiasAndResoHistos(h2_y);
    
    hn_bias_x.push_back( HistoAndName(p_x.first, names[i].second) );
    hn_bias_y.push_back( HistoAndName(p_y.first, names[i].second) );

    hn_reso_x.push_back( HistoAndName(p_x.second, names[i].second) );
    hn_reso_y.push_back( HistoAndName(p_y.second, names[i].second) );

  }


  drawProjections( outputdir, name, hn_x, "x" );
  drawProjections( outputdir, name, hn_y, "y" );


  if( name!="vert" ) drawPerformancePlot( outputdir, name, "x", hn_bias_x, "bias", Form("%s x Position [mm]", axisName.c_str()), -15., 15. );
  if( name!="horiz") drawPerformancePlot( outputdir, name, "y", hn_bias_y, "bias", Form("%s y Position [mm]", axisName.c_str()), -15., 15. );
  if( name!="vert" ) drawPerformancePlot( outputdir, name, "x", hn_reso_x, "reso", Form("%s x Position [mm]", axisName.c_str()), -15., 15. );
  if( name!="horiz") drawPerformancePlot( outputdir, name, "y", hn_reso_y, "reso", Form("%s y Position [mm]", axisName.c_str()), -15., 15. );

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



std::pair<TH1D*,TH1D*> getBiasAndResoHistosHodo( const std::string& outputdir, TH2D* h2 ) {

  int  nBins = h2->GetXaxis()->GetNbins();
  float xMin = h2->GetXaxis()->GetXmin();
  float xMax = h2->GetXaxis()->GetXmax();

  TH1D* h1_bias = new TH1D( Form("bias_%s", h2->GetName()), "", nBins, xMin, xMax );
  TH1D* h1_reso = new TH1D( Form("reso_%s", h2->GetName()), "", nBins, xMin, xMax );

  for( unsigned i=0; i<nBins; ++i ) {

    int iBin = i+1;
    TH1D* h1_proj = h2->ProjectionY(Form("%s_%d", h2->GetName(), iBin), iBin, iBin);

    if( h1_proj->GetEntries()>0 ) {

      TF1* f1 = fitHodoWithBeam( outputdir, h1_proj->GetName(), h1_proj, h2->GetXaxis()->GetBinCenter(iBin), 3. );
 
      h1_bias->SetBinContent( iBin, f1->GetParameter(1) );
      h1_bias->SetBinError  ( iBin, f1->GetParError(1) );
      h1_reso->SetBinContent( iBin, f1->GetParameter(2) );
      h1_reso->SetBinError  ( iBin, f1->GetParError(2) );

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
  colors.push_back( 46 );
  colors.push_back( 38 );
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
  colors.push_back( 46 );
  colors.push_back( 38 );
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

  float yMax_leg = 0.91;
  float yMin_leg = yMax_leg - (float)(hn.size())*0.06;

  TLegend* legend = new TLegend( 0.2, yMin_leg, 0.4, yMax_leg );
  legend->SetFillColor(0);
  legend->SetTextSize(0.033);

  float yMax = -999.;
  float yMin = +999.;

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

    if( hn[i].histo->GetMaximum() > yMax ) yMax = hn[i].histo->GetMaximum();
    if( hn[i].histo->GetMinimum() < yMin ) yMin = hn[i].histo->GetMinimum();
  
  }


  if( isReso ) yMin = 0.;
  else yMin -= 1.;
  if( isReso ) yMax += 6.;
  else yMax += 2.;
 

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, yMin, yMax );
  h2_axes->SetXTitle( axisName.c_str() );
  if(isReso)
    h2_axes->SetYTitle( "RMS(Estimate - Reference) [mm]" );
  else
    h2_axes->SetYTitle( "Estimate - Reference [mm]" );

  

  h2_axes->Draw();

  legend->Draw("same");

  if( isReso ) {

    float nominalReso = (var=="x") ? 4. : 3.;
  
    TLine* line_nominalReso = new TLine( xMin, nominalReso, xMax, nominalReso );
    line_nominalReso->SetLineStyle(2);
    line_nominalReso->Draw("same");

  }


  for( unsigned i=0; i<hn.size(); ++i )
    hn[i].histo->Draw("p same");


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



TF1* fitHodoWithBeam( const std::string& outputdir, const std::string& suffix, TH1D* h1, float bias, float r ) {

  //h1->Rebin(5);
  float nentries =  h1->GetEntries();

  float xFitMin = -bias-4.5;
  float xFitMax = -bias+4.5;

  TF1* f1_gaus = new TF1( Form("f1_%s", suffix.c_str()), "gaus" );
  f1_gaus->SetRange(xFitMin, xFitMax);
  //f1_gaus->SetParameter(0, nentries);
  f1_gaus->FixParameter(0, nentries);
  f1_gaus->SetParameter(1, 0.);

  //f1_gaus->SetParLimits(0, 0.01*nentries, nentries);
  f1_gaus->SetParLimits(1, -40., 40.);
  f1_gaus->SetParameter(2, r );

  h1->Fit(f1_gaus, "RL+" );
  TCanvas* c1 = new TCanvas("c1_temp", "", 600, 600);
  c1->cd();
  h1->SetLineColor(kRed);
  h1->SetLineWidth(2);
  h1->SetXTitle("Hodo Position - Beam [mm]");
  h1->SetYTitle("Entries");
  h1->Draw("");

  //f1_gaus->SetRange(-40., 40.);
  //f1_gaus->Draw("same");
  c1->SaveAs(Form("%s/tmpFit%s.eps", outputdir.c_str(), suffix.c_str()));
  c1->SaveAs(Form("%s/tmpFit%s.pdf", outputdir.c_str(), suffix.c_str()));
  c1->SaveAs(Form("%s/tmpFit%s.png", outputdir.c_str(), suffix.c_str()));
  delete c1;


  return f1_gaus;

}

