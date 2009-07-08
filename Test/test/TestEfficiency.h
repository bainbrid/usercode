#ifndef bainbrid_Test_TestEfficiency_h
#define bainbrid_Test_TestEfficiency_h

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

class TestEfficiency : public edm::EDAnalyzer {

 public:
  
  explicit TestEfficiency( const edm::ParameterSet& );
  ~TestEfficiency() {;}

  typedef std::vector<std::string> vstring;
  
 private:
  
  void beginJob( const edm::EventSetup& );
  void analyze( const edm::Event&, const edm::EventSetup& ) {;}
  void efficiencyPlots( TFile&, const vstring& );
  void efficiencyCurves( TFile& );
  void significanceCurves( TFile& );

 private:
  
  std::string append_;
  std::string file_;
  vstring signal_;
  vstring bkgd_;
  vstring histos_;
  
};

#endif // bainbrid_Test_TestEfficiency_h
