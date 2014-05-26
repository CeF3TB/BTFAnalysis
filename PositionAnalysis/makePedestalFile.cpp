#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"


#include "fastDQM_CeF3_BTF.h"






std::vector<TH1D*> getHistoVector(const std::string name, int nChannels);

void setPedestal( const std::string& fileName, int runNumber, const std::string& type, std::vector<TH1D*> vh_peds );

void savePedestalsToFile( TFile* file, std::vector<TH1D*> vh_peds );




int main( int argc, char* argv[] ) {


  ifstream pedList("pedestalFilesList.txt");

  TFile* outfile = TFile::Open("pedestalFile.root", "recreate");
 
  std::vector<TH1D*> ped_hodox = getHistoVector("hodox", HODOX_CHANNELS);
  std::vector<TH1D*> ped_hodoy = getHistoVector("hodoy", HODOY_CHANNELS);
  std::vector<TH1D*> ped_bgo   = getHistoVector("bgo",   BGO_CHANNELS);
  std::vector<TH1D*> ped_cef3  = getHistoVector("cef3" , CEF3_CHANNELS);

  

  while( pedList.good() ) {

    std::string fileName;
    pedList >> fileName;

    char runName[200];
    sscanf( fileName.c_str(), "data/run_%s_pedestal.root", runName );

    int runNumber, date, time;
    sscanf( runName, "BTF_%d_%d-%d", &runNumber, &date, &time );

    setPedestal( fileName, runNumber, "hodox", ped_hodox );
    setPedestal( fileName, runNumber, "hodoy", ped_hodoy );
    setPedestal( fileName, runNumber, "bgo", ped_bgo );
    setPedestal( fileName, runNumber, "cef3", ped_cef3 );

  }

  outfile->cd();

  savePedestalsToFile( outfile, ped_hodox );  
  savePedestalsToFile( outfile, ped_hodoy );  
  savePedestalsToFile( outfile, ped_bgo );  
  savePedestalsToFile( outfile, ped_cef3 );  

  outfile->Close();

  std::cout << "-> Saved pedestals in: " << outfile->GetName() << std::endl;

  return 0;

}









void setPedestal( const std::string& fileName, int runNumber, const std::string& type, std::vector<TH1D*> vh_peds ) {

  int nBoard=-1;
  int nChannels=-1;
  int firstChannel=-1;
  if( type=="cef3" ) {
    nBoard       = CEF3_ADC_BOARD;
    nChannels    = CEF3_CHANNELS;
    firstChannel = CEF3_ADC_START_CHANNEL;
  } else if( type=="bgo" ) {
    nBoard       = BGO_ADC_BOARD;
    nChannels    = BGO_CHANNELS;
    firstChannel = BGO_ADC_START_CHANNEL;
  } else if( type=="hodox" ) {
    nBoard       = HODOX_ADC_BOARD;
    nChannels    = HODOX_CHANNELS;
    firstChannel = HODOX_ADC_START_CHANNEL;
  } else if( type=="hodoy" ) {
    nBoard       = HODOY_ADC_BOARD;
    nChannels    = HODOY_CHANNELS;
    firstChannel = HODOY_ADC_START_CHANNEL;
  } else {
    std::cout << "ERROR! Unkown type '" << type << "'!" << std::endl;
    std::cout << "Don't know what pedestals you're looking for." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(77);
  }

    

  TFile* file = TFile::Open(fileName.c_str());
  if( file==0 ) return;
  TTree* tree = (TTree*)file->Get("eventRawData");

  std::vector< std::pair<float, float> > peds;
  for( unsigned i=0; i<nChannels; ++i ) {
    int iChannel = firstChannel+i;
    TH1D* h1_ped = new TH1D("ped", "", 500, 0., 500.);
    tree->Project( "ped", "adcData", Form("adcBoard==%d && adcChannel==%d", nBoard, iChannel) );
    vh_peds[i]->SetBinContent( runNumber, h1_ped->GetMean() );
    vh_peds[i]->SetBinError( runNumber, h1_ped->GetRMS() );
    delete h1_ped;
  }


}




std::vector<TH1D*> getHistoVector(const std::string name, int nChannels) {

  std::vector<TH1D*> returnVector;
  int nBins = 500;

  for( unsigned i=0; i<nChannels; ++i ) {
    TH1D* newHisto = new TH1D(Form("%s_%d", name.c_str(), i), "", nBins, 0., (float)nBins);
    for(unsigned ibin=1;ibin<nBins+1;++ibin ) {
      newHisto->SetBinContent( ibin, -1. );
      newHisto->SetBinError( ibin, 0.);
    }
    returnVector.push_back(newHisto);
  }

  return returnVector;

}



void savePedestalsToFile( TFile* file, std::vector<TH1D*> vh_peds ) {


  file->cd();

  for( unsigned i=0; i<vh_peds.size(); ++i )
    vh_peds[i]->Write();

}
