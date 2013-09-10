#include "mSuGraPlottingOps.hh"
#include "NLOTools.hh"
#include "AlphaT.hh"
#include "CommonOps.hh"
#include "EventData.hh"
#include "GenMatrixBin.hh"
#include "Jet.hh"
#include "KinSuite.hh"
#include "Lepton.hh"
#include "Math/VectorUtil.h"
#include "Photon.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "Types.hh"
#include <sstream>
#include <string>
#include <iomanip>

#include <strstream>
#include <iostream>
#include <fstream>


#include "JetData.hh"
#include "CommonOps.hh"


using namespace Operation;




//----------------------------------------------------------------------
mSuGraPlottingOps::mSuGraPlottingOps(  const Utils::ParameterSet& ps ) :
  dirName_( ps.Get<std::string>("DirectoryName")),
  mSuGraFile_(ps.Get<std::string>("mSUGRAFile")),
  xBins_(ps.Get<int>("xBins")),
  xLow_(ps.Get<double>("xLow")),
  xHigh_(ps.Get<double>("xHigh")),
  yBins_(ps.Get<int>("yBins")),
  yLow_(ps.Get<double>("yLow")),
  yHigh_(ps.Get<double>("yHigh")),
  zBins_(ps.Get<int>("zBins")),
  zLow_(ps.Get<double>("zLow")),
  zHigh_(ps.Get<double>("zHigh")),
  verbose_(ps.Get<bool>("verbose"))
{
  
}

mSuGraPlottingOps::~mSuGraPlottingOps(){}

//
void mSuGraPlottingOps::Start( Event::Data& ev ) {


  initDir( ev.OutputFile(), dirName_.c_str() );

  BookHistos();

  
 
}

//
void mSuGraPlottingOps::BookHistos() {

  BookHistArray( H_M0_M12_mChi,
		 "m0_m12_mChi",
		 ";m0;m12;mChi;",
		 xBins_,xLow_,xHigh_, //m0
		 yBins_,yLow_,yHigh_, //m12
		 zBins_,zLow_,zHigh_, //mChi
		 1,0,1,false);
  
  BookHistArray( H_M0_M12_mChi_noweight,
		 "m0_m12_mChi_noweight",
		 ";m0;m12;mChi;",
		 xBins_,xLow_,xHigh_, //m0
		 yBins_,yLow_,yHigh_, //m12
		 zBins_,zLow_,zHigh_, //mChi
		 1,0,1,false);
  
  BookHistArray( H_M0_M12_sb,
		 "m0_m12_sb",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );
  
  BookHistArray( H_M0_M12_sb_noweight,
		 "m0_m12_sb_noweight",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );
  
  BookHistArray( H_M0_M12_ss,
		 "m0_m12_ss",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );
  
  BookHistArray( H_M0_M12_ss_noweight,
		 "m0_m12_ss_noweight",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );
  
  BookHistArray( H_M0_M12_sg,
		 "m0_m12_sg",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );
  
  BookHistArray( H_M0_M12_sg_noweight,
		 "m0_m12_sg_noweight",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );

  BookHistArray( H_M0_M12_gg,
		 "m0_m12_gg",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );

  BookHistArray( H_M0_M12_gg_noweight,
		 "m0_m12_gg_noweight",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );


  BookHistArray( H_M0_M12_ll,
		 "m0_m12_ll",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );

  BookHistArray( H_M0_M12_ll_noweight,
		 "m0_m12_ll_noweight",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );

  BookHistArray( H_M0_M12_nn,
		 "m0_m12_nn",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );

  BookHistArray( H_M0_M12_nn_noweight,
		 "m0_m12_nn_noweight",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );

  BookHistArray( H_M0_M12_ng,
		 "m0_m12_ng",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );

  BookHistArray( H_M0_M12_ng_noweight,
		 "m0_m12_ng_noweight",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );

  BookHistArray( H_M0_M12_ns,
		 "m0_m12_ns",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );

  BookHistArray( H_M0_M12_ns_noweight,
		 "m0_m12_ns_noweight",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );


  BookHistArray( H_M0_M12_bb,
		 "m0_m12_bb",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );

   BookHistArray( H_M0_M12_bb_noweight,
		  "m0_m12_bb_noweight",
		  ";m0;m12",
		  xBins_,xLow_,xHigh_,//m0
		  yBins_,yLow_,yHigh_,//m12
		  1, 0, 1,false );

  BookHistArray( H_M0_M12_tb,
		 "m0_m12_tb",
		 ";m0;m12",
		 xBins_,xLow_,xHigh_,//m0
		 yBins_,yLow_,yHigh_,//m12
		 1, 0, 1,false );

   BookHistArray( H_M0_M12_tb_noweight,
		  "m0_m12_tb_noweight",
		  ";m0;m12",
		  xBins_,xLow_,xHigh_,//m0
		  yBins_,yLow_,yHigh_,//m12
		  1, 0, 1,false );


  



}




bool mSuGraPlottingOps::Process( Event::Data& ev ) {
  // Event weight
  Double_t weight = ev.GetEventWeight();
  Double_t nlocross = ev.GetEventWeight();

  string kFactorFile = mSuGraFile_;
  TString process = "";

  double M0 = 0.;
  double M12 = 0.;
  double MChi = 0.;

  if(ev.M0.enabled()){
    M0 = ev.M0();
  }
  if(ev.MG.enabled()){
    M0 = ev.MG();
  }
  if(ev.M12.enabled()){
    M12 = ev.M12();
  }
  if(ev.MLSP.enabled()){
    M12 = ev.MLSP();
  }
  if(ev.MChi.enabled()){
    MChi = ev.MChi();
  }

   //all jet mult
 
  H_M0_M12_mChi[0]->Fill(M0,M12,MChi,weight);
  H_M0_M12_mChi_noweight[0]->Fill(M0,M12,MChi,1);





  //NLO stuff: for the calculation of the NLO cross-section the processes are filled separately
  if((ev.M0.enabled() || ev.MG.enabled())){

     process = NLO::GetProcess(ev);
 
     Double_t NLOcrosssection = 1;
    
       
     NLOcrosssection = NLO::GetNLOCross(M0,M12,process,kFactorFile);     
     if(NLOcrosssection != -1) nlocross = NLOcrosssection;//if no valid kFactorFile has been given simply fill LO eventweight
    
     if(verbose_ == true)cout << " m0 " << M0 << " m12 " << M12 << " nlocross " << nlocross << endl;
             
     if(process == "nn"){H_M0_M12_nn[0]->Fill(M0,M12,nlocross); H_M0_M12_nn_noweight[0]->Fill(M0,M12,1);}
     if(process == "ns"){H_M0_M12_ns[0]->Fill(M0,M12,nlocross); H_M0_M12_ns_noweight[0]->Fill(M0,M12,1);}
     if(process == "ng"){H_M0_M12_ng[0]->Fill(M0,M12,nlocross); H_M0_M12_ng_noweight[0]->Fill(M0,M12,1);}
     if(process == "ss"){H_M0_M12_ss[0]->Fill(M0,M12,nlocross); H_M0_M12_ss_noweight[0]->Fill(M0,M12,1);}
     if(process == "ll"){H_M0_M12_ll[0]->Fill(M0,M12,nlocross); H_M0_M12_ll_noweight[0]->Fill(M0,M12,1);}
     if(process == "sb"){H_M0_M12_sb[0]->Fill(M0,M12,nlocross); H_M0_M12_sb_noweight[0]->Fill(M0,M12,1);}
     if(process == "tb"){H_M0_M12_tb[0]->Fill(M0,M12,nlocross); H_M0_M12_tb_noweight[0]->Fill(M0,M12,1);}
     if(process == "gg"){H_M0_M12_gg[0]->Fill(M0,M12,nlocross); H_M0_M12_gg_noweight[0]->Fill(M0,M12,1);}
     if(process == "bb"){H_M0_M12_bb[0]->Fill(M0,M12,nlocross); H_M0_M12_bb_noweight[0]->Fill(M0,M12,1);}
     if(process == "sg"){H_M0_M12_sg[0]->Fill(M0,M12,nlocross); H_M0_M12_sg_noweight[0]->Fill(M0,M12,1);}
     
  }
  return true;
    
}


//
std::ostream& mSuGraPlottingOps::Description( std::ostream& ostrm ) {
  ostrm << "mSuGra scan 2d Plots ";
  return ostrm;
}
