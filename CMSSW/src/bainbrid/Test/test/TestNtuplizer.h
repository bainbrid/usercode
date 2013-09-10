#ifndef bainbrid_Test_TestNtuplizer_h
#define bainbrid_Test_TestNtuplizer_h

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include <vector>

class TTree;

/**
   
 */
class TestNtuplizer : public edm::EDAnalyzer {

 public:
  
  explicit TestNtuplizer( const edm::ParameterSet& );
  ~TestNtuplizer() {;}
  
 private:

  void beginJob( const edm::EventSetup& );
  void analyze( const edm::Event&, const edm::EventSetup& );
  
  TTree* tree_;
  int n_;
  double e_[50];
  
};


#endif // bainbrid_Test_TestNtuplizer_h
