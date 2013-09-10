#ifndef bainbrid_Test_SimplePhotonIDAnalysis_h
#define bainbrid_Test_SimplePhotonIDAnalysis_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "boost/cstdint.hpp"
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
  
  void beginJob( const edm::EventSetup& );
  void endJob();

  void analyze( const edm::Event&, const edm::EventSetup& );
  
  // Histograms
  std::map<std::string,TH1D*> histos_; 
  std::map<std::string,TH2D*> histos2d_; 

  // Categories of photons
  std::vector<std::string> labels_;

  // InputTags for collections
  edm::InputTag photons_;
  edm::InputTag others_;
  edm::InputTag gen_;

  typedef std::map<int32_t,uint32_t> IdMap;
  typedef std::vector<IdMap> IdMaps;
  
  IdMaps genDaughterPdgId_;
  IdMaps genMotherPdgId_;
  IdMaps matchedDaughterPdgId_;
  IdMaps matchedMotherPdgId_;

  std::vector<int32_t> phoIds_;
  std::vector<int32_t> eleIds_;
  std::vector<int32_t> allIds_;
  
  static int   nbins_;
  static float lower_;
  static float upper_;
  static int32_t invalidId_;

};

#endif // bainbrid_Test_SimplePhotonIDAnalysis_h
