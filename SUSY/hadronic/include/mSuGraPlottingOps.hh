#ifndef MSUGRAPLOTTINGOPS
#define MSUGRAPLOTTINGOPS

#include "PlottingBase.hh"
#include "Utils.hh"
#include "Types.hh"
#include "Math/VectorUtil.h"
#include "ThrustStuff.hh"
#include "TH3.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

#include <strstream>
#include <iostream>
#include <fstream>

class TH1D;
class TH2D;



namespace Operation {

   class mSuGraPlottingOps : public PlottingBase {


  public:

    mSuGraPlottingOps( const Utils::ParameterSet& );
    ~mSuGraPlottingOps();

     void Start( Event::Data& );
     bool Process( Event::Data& );

   



     std::ostream& Description( std::ostream& );

   private:

     void BookHistos();
     
     std::string dirName_;
     std::string mSuGraFile_;
     
     bool verbose_;
     
     
     
     std::vector<TH3D *> H_M0_M12_mChi;   std::vector<TH3D *> H_M0_M12_mChi_noweight;
     std::vector<TH2D *> H_M0_M12_nn;     std::vector<TH2D *> H_M0_M12_nn_noweight;
     std::vector<TH2D *> H_M0_M12_ll;     std::vector<TH2D *> H_M0_M12_ll_noweight;
     std::vector<TH2D *> H_M0_M12_tb;     std::vector<TH2D *> H_M0_M12_tb_noweight;
     std::vector<TH2D *> H_M0_M12_bb;     std::vector<TH2D *> H_M0_M12_bb_noweight;
     std::vector<TH2D *> H_M0_M12_ns;     std::vector<TH2D *> H_M0_M12_ns_noweight;
     std::vector<TH2D *> H_M0_M12_ng;     std::vector<TH2D *> H_M0_M12_ng_noweight;
     std::vector<TH2D *> H_M0_M12_ss;     std::vector<TH2D *> H_M0_M12_ss_noweight;
     std::vector<TH2D *> H_M0_M12_sg;     std::vector<TH2D *> H_M0_M12_sg_noweight;
     std::vector<TH2D *> H_M0_M12_sb;     std::vector<TH2D *> H_M0_M12_sb_noweight;
     std::vector<TH2D *> H_M0_M12_gg;     std::vector<TH2D *> H_M0_M12_gg_noweight;


   
    int    xBins_;
    double xLow_;
    double xHigh_;
    int    yBins_;
    double yLow_;
    double yHigh_;
    int    zBins_;
    double zLow_;
    double zHigh_;


  }; //~mSuGra scan class
}
#endif 
