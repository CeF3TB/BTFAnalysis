#ifndef HodoCluster_h
#define HodoCluster_h

class HodoCluster {

 public:

  HodoCluster() {
    size_ = 0;
    pos_ = 0.;
  }

  ~HodoCluster() {};

  float getSize() { return size_; };
  float getPosition() { return pos_; };

  void addFibre( int i );


 private:

  int size_;
  float pos_;

};


#endif
