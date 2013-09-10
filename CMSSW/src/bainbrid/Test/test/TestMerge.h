#ifndef bainbrid_Test_TestMerge_h
#define bainbrid_Test_TestMerge_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <vector>
#include <string>
#include <map>

class TDirectory;
class TFile;
class TList;
class TH1D;
class TH2D;

class TestMerge : public edm::EDAnalyzer {

 public:
  
  explicit TestMerge( const edm::ParameterSet& );
  ~TestMerge() {;}

  typedef edm::ParameterSet PSet;
  typedef std::vector<PSet> VPSet;
  
 private:
  
  void beginJob( const edm::EventSetup& );
  void analyze( const edm::Event&, const edm::EventSetup& ) {;}
  void merge( TDirectory*, TList* );
  
  void files();
  void files( const PSet& );
  
  TH1D* histo( const std::string& histogram_name );
  TH2D* histo2d( const std::string& histogram_name );

 private:
  
  class Data {
  public:
    Data( std::string, float );
    Data();
    std::string name_;
    float xSec_;
  };
  
  std::map<TFile*,Data> files_;

  std::string file_;

  std::map<std::string,TH1D*> histos_; 
  std::map<std::string,TH2D*> histos2d_; 
  
  float lumi_;
  
};

#endif // bainbrid_Test_TestMerge_h
