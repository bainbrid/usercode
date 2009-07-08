#ifndef bainbrid_Test_TestCombination_h
#define bainbrid_Test_TestCombination_h

#include "DataFormats/Candidate/interface/Particle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include <vector>
#include <string>
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
  double minDHT( const std::vector<reco::Particle>&,
		 std::vector<reco::Particle>&, 
		 std::vector<reco::Particle>& );

  double recoilDHT( const std::vector<reco::Particle>& photons,
		    const std::vector<reco::Particle>& jets,
		    std::vector<reco::Particle>&, 
		    std::vector<reco::Particle>& );
  
  double minSuperDHT( const std::vector<reco::Particle>&,
		      std::vector<reco::Particle>&, 
		      std::vector<reco::Particle>& );
  
  double superDPHI( const std::vector<reco::Particle>&,
		    std::vector<reco::Particle>&, 
		    std::vector<reco::Particle>& );
  
  double HT( const std::vector<reco::Particle>& );
  
  double MHT( const std::vector<reco::Particle>& );
  
  double MT( const std::vector<reco::Particle>& );
  
  double HT_MHT( const std::vector<reco::Particle>& );
  
  double alphaT( double minDeltaEt, double sumEt, double massT );
  
  double alphaT( const std::vector<reco::Particle>&,
		 std::vector<reco::Particle>&, 
		 std::vector<reco::Particle>& );
  
  // Retrieve objects
  bool getPhotons( const edm::Event&, std::vector<reco::Particle>&, double threshold = -1. );
  bool getJets( const edm::Event&, std::vector<reco::Particle>&, double threshold = -1. );
  bool getMuons( const edm::Event&, std::vector<reco::Particle>& );
  bool getElectrons( const edm::Event&, std::vector<reco::Particle>& );

  TH1D* histo( const std::string& histogram_name );
  TH2D* histo2d( const std::string& histogram_name );

 private:

  typedef reco::Particle::LorentzVector LorentzVector;

  double Et( const LorentzVector& );

  LorentzVector t4Mom( const LorentzVector& );

  // Maximum number of objects handled
  int maximum_;

  // Number of objects to test
  int test_;
  
  // Input tags for collections
  edm::InputTag photons_;
  edm::InputTag jets_;
  edm::InputTag muons_;
  edm::InputTag electrons_;
  edm::InputTag met_;  
  edm::InputTag ccMet_;  
  edm::InputTag genMet_;  
  edm::InputTag gen_;  

  // Some kinematic thresholds
  double photonEt_;
  double photonEta_;
  double jetEt_;
  double jetEta_;
  double jetEMfrac_;
  double muonPt_;
  double muonEta_;
  double muonTrkIso_;
  double electronPt_;
  double electronEta_;
  double electronTrkIso_;
  
  // Event selection
  double totalHt_;
  uint32_t nObjects_;
  uint32_t nPhotons_;
  uint32_t minJets_;
  uint32_t maxJets_;

  // Histograms
  std::map<std::string,TH1D*> histos_; 
  std::map<std::string,TH2D*> histos2d_; 
  
};

inline double TestCombination::Et( const LorentzVector& input ) {
  double et = input.E() * input.E() - input.Pz() * input.Pz();
  et = et < 0. ? -1.*sqrt(-1.*et) : sqrt(et);
  return et;
}

inline TestCombination::LorentzVector TestCombination::t4Mom( const LorentzVector& input ) {
  return LorentzVector( input.Px(), input.Py(), 0., Et(input) );
}

#endif // bainbrid_Test_TestCombination_h








/*   // Event weighting */
/*   bool  eventWeight_; */
/*   double normLumi_; */
/*   double xSec_; */
/*   int   nEvents_; */
/*   double weight_; */


/*   // Builds possible combinations */
/*   void combinations( uint8_t,  */
/* 		     std::vector< std::vector<uint8_t> >&,  */
/* 		     std::vector< std::vector<uint8_t> >&  ); */
  
/*   // Finds object combination with minimum deltaPt */
/*   double minDeltaPt( const std::vector<reco::Particle>&,  */
/* 		    const std::vector< std::vector<uint8_t> >&,  */
/* 		    const std::vector< std::vector<uint8_t> >&, */
/* 		    std::vector<reco::Particle>&,  */
/* 		    std::vector<reco::Particle>& ); */
  
/* inline int TestCombination::factorial( int n ) { */
/*   int fact = 1; */
/*   for ( int i=n; i>=1; i-- ) fact = fact * i; */
/*   return fact; */
/* } */
