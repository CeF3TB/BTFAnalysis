#ifndef RunHelper_h
#define RunHelper_h




#include <string>

class RunHelper {

 public:

  static void getBeamPosition( const std::string& runNamei, float& xBeam, float& yBeam );
  static void getBGOCoordinates( int iChannel, float& x , float& y );


};


#endif
