#ifndef PositionTools_h
#define PositionTools_h

#include <vector>




class PositionTools {

 public:

  static void getBGOCoordinates( int iChannel, float& x , float& y );
  static void getCef3Coordinates( int iChannel, float& x , float& y );

  static void getBGOPositionAsymm( std::vector<float> bgo, float& x, float& y, bool logWeights=false );
  static void getBGOPositionWA( std::vector<float> bgo, float& x, float& y, bool logWeights=false );
  static void getCaloPositionWA( std::vector<float> cef3, std::vector<float> bgo, float& x, float& y, bool logWeights=false );
  static void getCeF3Position( std::vector<float> cef3, float& x, float& y, int event );


 private:


};


#endif
