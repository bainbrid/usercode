#ifndef bainbrid_Test_RawPATJetProducer_h
#define bainbrid_Test_RawPATJetProducer_h

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class RawPATJetProducer : public edm::EDProducer {

 public:
  
  explicit RawPATJetProducer( const edm::ParameterSet& );
  ~RawPATJetProducer();

  virtual void produce( edm::Event&, const edm::EventSetup& );
  
 private:
  
  edm::InputTag jets_;
  
};

#endif // bainbrid_Test_RawPATJetProducer_h
