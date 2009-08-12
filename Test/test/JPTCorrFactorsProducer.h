#ifndef bainbrid_Test_JPTCorrFactorsProducer_h
#define bainbrid_Test_JPTCorrFactorsProducer_h

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class JPTCorrFactorsProducer : public edm::EDProducer {

 public:
  
  explicit JPTCorrFactorsProducer( const edm::ParameterSet& );
  ~JPTCorrFactorsProducer();

  virtual void produce( edm::Event&, const edm::EventSetup& );
  
  static bool matchJetsByCaloTowers( const pat::Jet&, 
				     const pat::Jet& );
  
 private:
  
  edm::InputTag uncorrected_;
  edm::InputTag corrected_;
  
};

#endif // bainbrid_Test_JPTCorrFactorsProducer_h
