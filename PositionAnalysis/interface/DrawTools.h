#ifndef DrawTools_h
#define DrawTools_h


#include "TStyle.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TROOT.h"

class DrawTools {

 public:

  static TStyle* setStyle();

  static TPaveText* getLabelTop( const std::string& text );
  static TPaveText* getLabelTop( int beamEnergy=491 );
  static TPaveText* getLabelTop_2D( const std::string& text );
  static TPaveText* getLabelTop_2D( int beamEnergy=491 );
  static TPaveText* getLabelRun( const std::string& runName, bool top=true );
  static TPaveText* getCef3LabelLeft();
  static TText* getAuthorsLabel();

};


#endif
