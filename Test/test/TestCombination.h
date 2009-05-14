#ifndef bainbrid_Test_TestCombination_h
#define bainbrid_Test_TestCombination_h

#include "DataFormats/Candidate/interface/Particle.h"
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

  // Calculate minimum difference in Et from all pseudo-dijet combinations
  float minDeltaEt( const std::vector<reco::Particle>&,
		    std::vector<reco::Particle>&, 
		    std::vector<reco::Particle>& );

  float sumEt( const std::vector<reco::Particle>& );
  
  float massT( const std::vector<reco::Particle>& );
  
  float alphaT( float minDeltaEt, float sumEt, float massT );

  float alphaT( const std::vector<reco::Particle>&,
		std::vector<reco::Particle>&, 
		std::vector<reco::Particle>& );
  
  // Retrieve objects
  bool getPhotons( const edm::Event&, std::vector<reco::Particle>& );
  bool getJets( const edm::Event&, std::vector<reco::Particle>& );
  bool getMuons( const edm::Event&, std::vector<reco::Particle>& );
  bool getElectrons( const edm::Event&, std::vector<reco::Particle>& );

 private:

  // Maximum number of objects handled
  int maximum_;

  // Number of objects to test
  int test_;
  
  // Input tags for collections
  edm::InputTag photons_;
  edm::InputTag jets_;
  edm::InputTag muons_;
  edm::InputTag electrons_;
  
  // Some kinematic thresholds
  float photonEt_;
  float photonEta_;
  float jetEt_;
  float jetEta_;
  float jetEMfrac_;
  float muonPt_;
  float muonEta_;
  float muonTrkIso_;
  float electronPt_;
  float electronEta_;
  float electronTrkIso_;

  // Event selection
  float totalEt_;

  // Histograms
  std::map<std::string,TH1D*> histos_; 
  std::map<std::string,TH2D*> histos2d_; 
  
};

#endif // bainbrid_Test_TestCombination_h



/*   // Builds possible combinations */
/*   void combinations( uint8_t,  */
/* 		     std::vector< std::vector<uint8_t> >&,  */
/* 		     std::vector< std::vector<uint8_t> >&  ); */
  
/*   // Finds object combination with minimum deltaPt */
/*   float minDeltaPt( const std::vector<reco::Particle>&,  */
/* 		    const std::vector< std::vector<uint8_t> >&,  */
/* 		    const std::vector< std::vector<uint8_t> >&, */
/* 		    std::vector<reco::Particle>&,  */
/* 		    std::vector<reco::Particle>& ); */
  
/* inline int TestCombination::factorial( int n ) { */
/*   int fact = 1; */
/*   for ( int i=n; i>=1; i-- ) fact = fact * i; */
/*   return fact; */
/* } */
