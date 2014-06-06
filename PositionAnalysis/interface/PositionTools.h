#ifndef PositionTools_h
#define PositionTools_h

#include <vector>




class PositionTools {

 public:

  static void getBGOCoordinates( int iChannel, float& x , float& y );

  static void getBGOPositionAsymm( std::vector<float> bgo, float& x, float& y, bool logWeights=false );
  static void getBGOPositionWA( std::vector<float> bgo, float& x, float& y, bool logWeights=false );


 private:


};


#endif
