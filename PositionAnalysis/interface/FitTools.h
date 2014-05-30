#ifndef FitTools_h
#define FitTools_h


#include "TF1.h"
#include "TH1D.h"
#include "TTree.h"

#include <string>



class FitTools {

 public:

  static void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name );

  static TF1* fitSingleElectronPeak( const std::string& outputdir, const std::string& name, TTree* tree );


 private:

};



#endif
