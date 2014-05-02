#include <TStyle.h>
#include <TCanvas.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"

#include "fastDQM_CeF3_BTF.h"

std::string runName;

int main( int argc, char* argv[] ) {

  runName = "test_10";
  if( argc>1 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
  }


  std::string fileName = "PosAn_" + runName + ".root";
  TFile* file = TFile::Open( fileName.c_str() );
  std::cout << "-> Opened file: " << fileName << std::endl;

  TTree* tree = (TTree*)file->Get("tree_passedEvents");

  float hodox_corr[HODOX_CHANNELS], hodoy_corr[HODOY_CHANNELS], cef3_corr[CEF3_CHANNELS];
  int nHodoFibersX,nHodoFibersY;

  tree->SetBranchAddress( "nHodoFibersX", &nHodoFibersX );
  tree->SetBranchAddress( "nHodoFibersY", &nHodoFibersY );
  tree->SetBranchAddress( "hodox_corr", hodox_corr );
  tree->SetBranchAddress( "hodoy_corr", hodoy_corr );
  tree->SetBranchAddress( "cef3_corr", cef3_corr );

  int nentries = tree->GetEntries();

  std::string outfileName = "HodoscopeAn_" + runName + ".root";
  TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );

  TH1D* h1_totalEnergyDistribution = new TH1D("totalEnergyDistribution", "",1600 , 0 , 16000);
  TH1D* h1_nFibersSingleXDistribution = new TH1D("nFibersSingleXDistribution", "",  10 , -0.5 , 10.5 );
  TH1D* h1_nFibersSingleYDistribution = new TH1D("nFibersSingleYDistribution", "",  10 , -0.5 , 10.5 );
  TH2D* h2_nFibersSingleXYDistribution = new TH2D("nFibersSingleXYDistribution", "",  10 , -0.5 , 10.5 ,  10 , -0.5 , 10.5);

  TH1D* h1_nFibersDoubleXDistribution = new TH1D("nFibersDoubleXDistribution", "",  10 , -0.5 , 10.5 );
  TH1D* h1_nFibersDoubleYDistribution = new TH1D("nFibersDoubleYDistribution", "",  10 , -0.5 , 10.5 );
  TH2D* h2_nFibersDoubleXYDistribution = new TH2D("nFibersDoubleXYDistribution", "",  10 , -0.5 , 10.5 ,  10 , -0.5 , 10.5);

  TH1D* h1_nFibersTripleXDistribution = new TH1D("nFibersTripleXDistribution", "",  10 , -0.5 , 10.5 );
  TH1D* h1_nFibersTripleYDistribution = new TH1D("nFibersTripleYDistribution", "",  10 , -0.5 , 10.5 );
  TH2D* h2_nFibersTripleXYDistribution = new TH2D("nFibersTripleXYDistribution", "",  10 , -0.5 , 10.5 ,  10 , -0.5 , 10.5);


  float cef3_corrLowerCutSingle=0;
  float cef3_corrUpperCutSingle=0;

  float cef3_corrLowerCutDouble=0;
  float cef3_corrUpperCutDouble=0;

  float cef3_corrLowerCutTriple=0;
  float cef3_corrUpperCutTriple=0;

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {
    
    tree->GetEntry(iEntry);

    if( iEntry % 1000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;

    float cef3_corrTotal=0.;    
    for (int i=0;i<CEF3_CHANNELS;i++){
      cef3_corrTotal+=cef3_corr[i];
    }
    //    std::cout<<cef3_corrTotal<<" ";
    h1_totalEnergyDistribution->Fill(cef3_corrTotal);
  }

  TF1* f1= new TF1("gausSum", "gaus(0)+gaus(3)+gaus(6)", 2000., 11000.);
  f1->SetParameter(0,1000);
  f1->SetParameter(1,3000.);
  f1->SetParameter(2,500);
  f1->SetParameter(3,100);
  f1->SetParameter(4,6500.);
  f1->SetParameter(5,700);
  f1->SetParameter(6,50);
  f1->SetParameter(7,9500.);
  f1->SetParameter(8,800);
  h1_totalEnergyDistribution->Fit("gausSum","LR+");

  float meanSingle= f1->GetParameter(1);
  float sigmaSingle= f1->GetParameter(2);
  std::cout<<"resolution: "<<sigmaSingle/meanSingle<<std::endl;
  std::cout<<"window single goes from "<<meanSingle-2*sigmaSingle<<"+/-"<<meanSingle+2*sigmaSingle<<std::endl;

  cef3_corrLowerCutSingle=meanSingle-2*sigmaSingle;
  cef3_corrUpperCutSingle=meanSingle+2*sigmaSingle;


  float meanDouble= f1->GetParameter(4);
  float sigmaDouble= f1->GetParameter(5);
  std::cout<<"resolution: "<<sigmaDouble/meanDouble<<std::endl;
  std::cout<<"window double goes from "<<meanDouble-2*sigmaDouble<<"+/-"<<meanDouble+2*sigmaDouble<<std::endl;

  cef3_corrLowerCutDouble=meanDouble-2*sigmaDouble;
  cef3_corrUpperCutDouble=meanDouble+2*sigmaDouble;

  float meanTriple= f1->GetParameter(7);
  float sigmaTriple= f1->GetParameter(8);
  std::cout<<"resolution: "<<sigmaTriple/meanTriple<<std::endl;
  std::cout<<"window triple goes from "<<meanTriple-2*sigmaTriple<<"+/-"<<meanTriple+2*sigmaTriple<<std::endl;

  cef3_corrLowerCutTriple=meanTriple-2*sigmaTriple;
  cef3_corrUpperCutTriple=meanTriple+2*sigmaTriple;


  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {
    
    tree->GetEntry(iEntry);

    if( iEntry % 1000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;

    float cef3_corrTotal=0.;    
    for (int i=0;i<CEF3_CHANNELS;i++){
      cef3_corrTotal+=cef3_corr[i];
    }

    if(cef3_corrTotal>cef3_corrLowerCutSingle && cef3_corrTotal<cef3_corrUpperCutSingle){
      h1_nFibersSingleXDistribution->Fill(nHodoFibersX);
      h1_nFibersSingleYDistribution->Fill(nHodoFibersY);
      h2_nFibersSingleXYDistribution->Fill(nHodoFibersX,nHodoFibersY);
    }else if(cef3_corrTotal>cef3_corrLowerCutDouble && cef3_corrTotal<cef3_corrUpperCutDouble){
      h1_nFibersDoubleXDistribution->Fill(nHodoFibersX);
      h1_nFibersDoubleYDistribution->Fill(nHodoFibersY);
      h2_nFibersDoubleXYDistribution->Fill(nHodoFibersX,nHodoFibersY);
    }else if(cef3_corrTotal>cef3_corrLowerCutTriple && cef3_corrTotal<cef3_corrUpperCutTriple){
      h1_nFibersTripleXDistribution->Fill(nHodoFibersX);
      h1_nFibersTripleYDistribution->Fill(nHodoFibersY);
      h2_nFibersTripleXYDistribution->Fill(nHodoFibersX,nHodoFibersY);
    }
  }




  outfile->cd();
  h1_totalEnergyDistribution->Write();
  h1_nFibersSingleXDistribution->Write();
  h1_nFibersSingleYDistribution->Write();  
  h2_nFibersSingleXYDistribution->Write();

  h1_nFibersDoubleXDistribution->Write();
  h1_nFibersDoubleYDistribution->Write();  
  h2_nFibersDoubleXYDistribution->Write();


  h1_nFibersTripleXDistribution->Write();
  h1_nFibersTripleYDistribution->Write();  
  h2_nFibersTripleXYDistribution->Write();


  outfile->Close();
  
  return 0;


}



