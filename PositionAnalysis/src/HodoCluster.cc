#include "../interface/HodoCluster.h"




void HodoCluster::addFibre( int i ) {

  if( size_==0 ) {

    pos_ += -(i-3.5);
    size_+=1;

  } else {

    pos_ *= (float)size_;
    pos_ += -(i-3.5);
    size_+=1;
    pos_ /= size_;

  }

}

