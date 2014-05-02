#ifndef DrawTools_h
#define DrawTools_h


#include "TStyle.h"
#include "TPaveText.h"



class DrawTools {

 public:

  static TStyle* setStyle();

  static TPaveText* getLabelTop( int beamEnergy=500 );
  static TPaveText* getLabelRun( const std::string& runName, bool top=true );


};


#endif
