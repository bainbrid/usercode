#ifndef bainbrid_Test_SimplePhotonIDAnalysis_h
#define bainbrid_Test_SimplePhotonIDAnalysis_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include <vector>
#include <string>
#include <map>

class TH1D;
class TH2D;

class SimplePhotonIDAnalysis : public edm::EDAnalyzer {

public:
  
  explicit SimplePhotonIDAnalysis( const edm::ParameterSet& );

  ~SimplePhotonIDAnalysis() {;}
  
private:
  
  virtual void beginJob( const edm::EventSetup& );

  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  
  // Histograms
  std::map<std::string,TH1D*> histos_; 
  std::map<std::string,TH2D*> histos2d_; 

  // Categories of photons
  std::vector<std::string> labels_;

  // InputTags for collections
  edm::InputTag photons_;
  edm::InputTag others_;

};

#endif // bainbrid_Test_SimplePhotonIDAnalysis_h
