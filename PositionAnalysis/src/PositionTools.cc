#include "../interface/PositionTools.h"
#include "../fastDQM_CeF3_BTF.h"
#include <iostream>
#include <cmath>






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
    exit(99);
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

  float left   = bgo[0]+bgo[3]+bgo[5];
  float right  = bgo[2]+bgo[4]+bgo[7];
  float top    = bgo[0]+bgo[1]+bgo[2];
  float bottom = bgo[5]+bgo[6]+bgo[7];

  if( logWeights ) {
    left = log(left);
    right = log(right);
    top = log(top);
    bottom = log(bottom);
  }

  float asymm_x = (right-left)/(right+left);
  float asymm_y = (top-bottom)/(top+bottom);

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

    float w = (logWeights) ? log(bgo[i]) : bgo[i];

    wa_x += w*xbgo;
    wa_y += w*ybgo;
    sumw += w;

  }

  wa_x /= sumw;
  wa_y /= sumw;

  float p0_x = (logWeights) ? 0.637256 : 2.42003;
  float p1_x = (logWeights) ? 0.221399 : 0.967115;
  
  float p0_y = (logWeights) ? 0.213944 : 0.409758;
  float p1_y = (logWeights) ? 0.228218 : 0.941134;

  x = (wa_x - p0_x)/p1_x;
  y = (wa_y - p0_y)/p1_y;

}
