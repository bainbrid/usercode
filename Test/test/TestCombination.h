#ifndef bainbrid_Test_TestCombination_h
#define bainbrid_Test_TestCombination_h

#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include <vector>
#include <map>

class TH1D;
class TH2D;

class TestCombination : public edm::EDAnalyzer {

 public:
  
  explicit TestCombination( const edm::ParameterSet& );

  ~TestCombination() {;}
  
 private:

  // Some init
  void beginJob( const edm::EventSetup& );

  // Top-level analysis method
  void analyze( const edm::Event&, const edm::EventSetup& );

  // Builds pseudo-jets
  float pseudoJets( const std::vector<reco::Particle>&,
		    std::vector<reco::Particle>&, 
		    std::vector<reco::Particle>&  );

  // Builds possible combinations
  void combinations( int, 
		     std::vector< std::vector<int> >&, 
		     std::vector< std::vector<int> >&  );

  // Finds object combination with minimum deltaPt
  float minDeltaPt( const std::vector<reco::Particle>&, 
		    const std::vector< std::vector<int> >&, 
		    const std::vector< std::vector<int> >&,
		    std::vector<reco::Particle>&, 
		    std::vector<reco::Particle>& );
  
  // Input tags for collections
  edm::InputTag jets_;
  edm::InputTag photons_;

  // Some thresholds
  float jetPt_;
  float jetEta_;
  float photonPt_;
  float photonEta_;
  float totalPt_;

  // Histograms
  std::map<std::string,TH1D*> histos_; 
  std::map<std::string,TH2D*> histos2d_; 
  
};
  
#endif // bainbrid_Test_TestCombination_h
