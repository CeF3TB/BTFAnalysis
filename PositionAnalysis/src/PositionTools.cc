#include "../interface/PositionTools.h"
#include "../fastDQM_CeF3_BTF.h"
#include <iostream>
#include <cmath>
#include <cstdlib>






void PositionTools::getBGOCoordinates( int iChannel, float& x , float& y ) {

  if( iChannel==0 ) {
    x = -532.0+511.; y = -178.0+201.;
  } else if( iChannel==1 ) {
    x = -510.0+511.; y = -178.0+201.;
  } else if( iChannel==2 ) {
    x = -488.0+511.; y = -180.0+201.;
  } else if( iChannel==3 ) {
    x = -534.0+511.; y = -200.1+201.;
  } else if( iChannel==4 ) {
    x = -488.0+511.; y = -202.1+201.;
  } else if( iChannel==5 ) {
    x = -534.0+511.; y = -222.1+201.;
  } else if( iChannel==6 ) {
    x = -512.0+511.; y = -224.1+201.;
  } else if( iChannel==7 ) {
    x = -490.0+511.; y = -224.1+201.;
  } else {
    std::cout << "[PositionTools::getBGOCoordinates] ERROR! BGO index must be between 0 and 7!" << std::endl;
    std::exit(99);
  } 


}



void PositionTools::getCef3Coordinates( int iChannel, float& x , float& y ) {

  float side = 12. - 0.696; // in mm

  if( iChannel==0 ) {
    x = -side; y = +side;
  } else if( iChannel==1 ) {
    x = +side; y = +side;
  } else if( iChannel==2 ) {
    x = +side; y = -side;
  } else if( iChannel==3 ) {
    x = -side; y = -side;
  } else {
    std::cout << "[PositionTools::getCef3Coordinates] ERROR! Cef3 index must be between 0 and 3!" << std::endl;
    exit(97);
  }

}



void PositionTools::getBGOPositionAsymm( std::vector<float> bgo, float& x, float& y, bool logWeights ) {

  x = 0.;
  y = 0.;

  if( bgo.size()!=BGO_CHANNELS ) {
    std::cout << "[PositionTools::getBGOPositionAsymm] ERROR! Expecting a BGO vector of size: " << BGO_CHANNELS << " (and not " << bgo.size() << ")!!" << std::endl;
    std::cout << "Returning (0,0)" << std::endl;
    return;
  }

  
  //  0 1 2 
  //  3   4
  //  5 6 7

  if( logWeights ) {
    bgo[0] = (bgo[0]>=1.) ? log( bgo[0] ) : 0.;
    bgo[1] = (bgo[1]>=1.) ? log( bgo[1] ) : 0.;
    bgo[2] = (bgo[2]>=1.) ? log( bgo[2] ) : 0.;
    bgo[3] = (bgo[3]>=1.) ? log( bgo[3] ) : 0.;
    bgo[4] = (bgo[4]>=1.) ? log( bgo[4] ) : 0.;
    bgo[5] = (bgo[5]>=1.) ? log( bgo[5] ) : 0.;
    bgo[6] = (bgo[6]>=1.) ? log( bgo[6] ) : 0.;
    bgo[7] = (bgo[7]>=1.) ? log( bgo[7] ) : 0.;
  }

  float left   = bgo[0]+bgo[3]+bgo[5];
  float right  = bgo[2]+bgo[4]+bgo[7];
  float top    = bgo[0]+bgo[1]+bgo[2];
  float bottom = bgo[5]+bgo[6]+bgo[7];

  //if( logWeights ) {
  //  left = log(left);
  //  right = log(right);
  //  top = log(top);
  //  bottom = log(bottom);
  //}

  float asymm_x = (right-left)/(right+left);
  float asymm_y = (top-bottom)/(top+bottom);

  //float p0_x = (logWeights) ? 0.0598468 : 0.153043 ;
  //float p1_x = (logWeights) ? 0.0348552 : 0.0711287;
  //
  //float p0_y = (logWeights) ? 0.0149668 : 0.0406251;
  //float p1_y = (logWeights) ? 0.0389653 : 0.0787041;
  
  float p0_x = (logWeights) ? 0.0373235 : 0.139943;
  float p1_x = (logWeights) ? 0.0130919 : 0.0532351;
  
  float p0_y = (logWeights) ? 0.0195486 : 0.0476113;
  float p1_y = (logWeights) ? 0.0136892 : 0.0556958;
  
  x = (asymm_x - p0_x)/p1_x;
  y = (asymm_y - p0_y)/p1_y;

}




void PositionTools::getBGOPositionWA( std::vector<float> bgo, float& x, float& y, bool logWeights ) {

  x = 0.;
  y = 0.;


  if( bgo.size()!=BGO_CHANNELS ) {
    std::cout << "[PositionTools::getBGOPositionWA] ERROR! Expecting a BGO vector of size: " << BGO_CHANNELS << " (and not " << bgo.size() << ")!!" << std::endl;
    std::cout << "Returning (0,0)" << std::endl;
    return;
  }


  float wa_x=0.;
  float wa_y=0.;
  float sumw = 0.;

  for( unsigned i=0; i<bgo.size(); ++i ) {

    float xbgo, ybgo;
    getBGOCoordinates( i, xbgo, ybgo );

    float logbgo = (bgo[i]>=1.) ? log(bgo[i]) : 0.;
    float w = (logWeights) ? logbgo : bgo[i];

    wa_x += w*xbgo;
    wa_y += w*ybgo;
    sumw += w;

  }

  wa_x /= sumw;
  wa_y /= sumw;

  //float p0_x = (logWeights) ? 0.988866 : 2.70895;
  //float p1_x = (logWeights) ? 0.555462 : 1.30758;
  //
  //float p0_y = (logWeights) ? 0.107338 : 0.240749;
  //float p1_y = (logWeights) ? 0.613692 : 1.34084;

  float p0_x = (logWeights) ? 0.637256 : 2.42003;
  float p1_x = (logWeights) ? 0.221399 : 0.967115;
  
  float p0_y = (logWeights) ? 0.213944 : 0.409758;
  float p1_y = (logWeights) ? 0.228218 : 0.941134;

  x = (wa_x - p0_x)/p1_x;
  y = (wa_y - p0_y)/p1_y;

}






void PositionTools::getCaloPositionWA( std::vector<float> cef3, std::vector<float> bgo, float& x, float& y, bool logWeights ) {

  x = 0.;
  y = 0.;


  if( cef3.size()!=CEF3_CHANNELS ) {
    std::cout << "[PositionTools::getBGOPositionWA] ERROR! Expecting a CEF3 vector of size: " << CEF3_CHANNELS << " (and not " << cef3.size() << ")!!" << std::endl;
    std::cout << "Returning (0,0)" << std::endl;
    return;
  }

  if( bgo.size()!=BGO_CHANNELS ) {
    std::cout << "[PositionTools::getBGOPositionWA] ERROR! Expecting a BGO vector of size: " << BGO_CHANNELS << " (and not " << bgo.size() << ")!!" << std::endl;
    std::cout << "Returning (0,0)" << std::endl;
    return;
  }


  float wa_x=0.;
  float wa_y=0.;
  float sumw = 0.;

  for( unsigned i=0; i<bgo.size(); ++i ) {

    float xbgo, ybgo;
    getBGOCoordinates( i, xbgo, ybgo );

    float logbgo = (bgo[i]>=1.) ? log(bgo[i]) : 0.;
    float w = (logWeights) ? logbgo : bgo[i];

    wa_x += w*xbgo;
    wa_y += w*ybgo;
    sumw += w;

  }

  float etot_cef3 = 0.;
  for( unsigned i=0; i<cef3.size(); ++i )  etot_cef3 += cef3[i];
  etot_cef3 *= 0.79; // calibrate to BGO

  if( logWeights ) {
    if( etot_cef3>=1. ) etot_cef3 = log(etot_cef3);
    else etot_cef3 = 0.;
  }
  
  sumw += etot_cef3;

  wa_x /= sumw;
  wa_y /= sumw;




  float p0_x = (logWeights) ? 0.783777 : 0.   ;
  float p1_x = (logWeights) ? 0.435678 : 1. ;
  
  float p0_y = (logWeights) ? 0.083535 : 0.;
  float p1_y = (logWeights) ? 0.479405 : 1. ;

  //float p0_x = (logWeights) ? 0.783777 : 0.9611   ;
  //float p1_x = (logWeights) ? 0.435678 : 0.432683 ;
  //
  //float p0_y = (logWeights) ? 0.083535 : 0.0307057;
  //float p1_y = (logWeights) ? 0.479405 : 0.376184 ;


  x = (wa_x - p0_x)/p1_x;
  y = (wa_y - p0_y)/p1_y;

}

#include "TLine.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TPaveText.h"

void PositionTools::getCeF3Position( std::vector<float> cef3, float& x, float& y, int event ) {


  float etot=0.;
  for( unsigned i=0; i<cef3.size(); ++i ) etot += cef3[i];


  std::vector<float> v_d;

  for( unsigned i=0; i<cef3.size(); ++i ) {
  
    float r_i = cef3[i]/etot;
std::cout << "cef3[" <<  i << "]: " << cef3[i] << " (" << 100.*r_i << "%%)" << std::endl;

    if( r_i>0.25 ) {
      float d_i = -sqrt( (r_i - 0.25)/0.000949544 ) + 16.;
      v_d.push_back( d_i );
    } else {
      v_d.push_back(-1.);
    }

  }  //for

  TCanvas* c1 = new TCanvas( "C11", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D( "axes", "", 10, -40., 40., 10, -40., 40. );
  h2_axes->Draw();

  float pos = 12. - 0.696; // in mm

  TLine* line_up = new TLine( -pos, pos, pos, pos );
  line_up->Draw("same");
  TLine* line_down = new TLine( -pos, -pos, pos, -pos );
  line_down->Draw("same");
  TLine* line_left = new TLine( -pos, -pos, -pos, pos );
  line_left->Draw("same");
  TLine* line_right = new TLine( pos, -pos, pos, pos );
  line_right->Draw("same");


  for( unsigned i=0; i<v_d.size(); ++i ) {
    if( v_d[i]<=0. ) continue;

    float x, y;
    getCef3Coordinates( i, x, y );
    TEllipse* ell = new TEllipse( x, y, v_d[i] );
    ell->SetFillStyle(0);
    ell->Draw("same");

  }


  TEllipse* beam = new TEllipse( -12., 0., 4., 3. );
  beam->SetLineColor(kRed);
  beam->SetFillStyle(0);
  beam->Draw("same");

  TPaveText* text0 = new TPaveText( -pos-0.5, pos+1., -pos+0.5, pos+2. );
  text0->SetFillColor(0);
  text0->SetTextSize(0.027);
  text0->AddText( Form("%.1f%%", 100.*cef3[0]/(etot)));
  text0->Draw("same");

  TPaveText* text1 = new TPaveText( +pos-0.5, pos+1., +pos+0.5, pos+2. );
  text1->SetFillColor(0);
  text1->SetTextSize(0.027);
  text1->AddText( Form("%.1f%%", 100.*cef3[1]/(etot)));
  text1->Draw("same");

  TPaveText* text2 = new TPaveText( +pos-0.5, -pos-1., +pos+0.5, -pos-2. );
  text2->SetFillColor(0);
  text2->SetTextSize(0.027);
  text2->AddText( Form("%.1f%%", 100.*cef3[2]/(etot)));
  text2->Draw("same");

  TPaveText* text3 = new TPaveText( -pos-0.5, -pos-1., -pos+0.5, -pos-2. );
  text3->SetFillColor(0);
  text3->SetTextSize(0.027);
  text3->AddText( Form("%.1f%%", 100.*cef3[3]/(etot)));
  text3->Draw("same");


  c1->SaveAs(Form("prova_%d.png", event));
    
  delete c1;
  delete h2_axes;

}

