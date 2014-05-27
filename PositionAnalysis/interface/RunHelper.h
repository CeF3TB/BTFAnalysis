#ifndef RunHelper_h
#define RunHelper_h



#include "TString.h"
#include <string>

class RunHelper {

 public:

  static void getBeamPosition( const std::string& runNamei, float& xBeam, float& yBeam );
  static void getBGOCoordinates( int iChannel, float& x , float& y );
  void setCalibrationFiles(std::string tagName);
  std::string getCalibrationFiles(bool isCef3);

  std::string cef3CalibrationVersion_;
  std::string bgoCalibrationVersion_;

};


#endif
