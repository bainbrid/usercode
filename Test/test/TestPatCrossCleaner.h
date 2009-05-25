#ifndef bainbrid_Test_TestPatCrossCleaner_h
#define bainbrid_Test_TestPatCrossCleaner_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include <vector>
#include <map>

class TH1D;
class TH2D;
class TH1;

class TestPatCrossCleaner : public edm::EDAnalyzer {

 public:
  
  explicit TestPatCrossCleaner( const edm::ParameterSet& );

  ~TestPatCrossCleaner() {;}
  
 private:

  // Some init
  void beginJob( const edm::EventSetup& );

  // Top-level analysis method
  void analyze( const edm::Event&, const edm::EventSetup& );

  TH1* histo( const std::string& histogram_name );

 private:

  // Input tags for collections
  edm::InputTag rawPhotons_;
  edm::InputTag rawJets_;
  edm::InputTag rawMuons_;
  edm::InputTag rawElectrons_;

  edm::InputTag ccPhotons_;
  edm::InputTag ccJets_;
  edm::InputTag ccMuons_;
  edm::InputTag ccElectrons_;

  edm::InputTag droppedPhotons_;
  edm::InputTag droppedJets_;
  edm::InputTag droppedMuons_;
  edm::InputTag droppedElectrons_;

  edm::InputTag jetMatchedCCPhotons_;
  edm::InputTag jetMatchedCCElectrons_;
  edm::InputTag jetMatchedCCPartons_;

  edm::InputTag jetMatchedDroppedPhotons_;
  edm::InputTag jetMatchedDroppedElectrons_;
  edm::InputTag jetMatchedDroppedPartons_;

  edm::InputTag photonMatchedCCPhotons_;
  edm::InputTag photonMatchedDroppedPhotons_;

  // Histograms
  std::map<std::string,TH1D*> histos_; 
  std::map<std::string,TH2D*> histos2d_; 
  
};

#endif // bainbrid_Test_TestPatCrossCleaner_h
