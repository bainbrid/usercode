#ifndef PhysicsTools_PatAlgos_PATPhotonIDProducer_h
#define PhysicsTools_PatAlgos_PATPhotonIDProducer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "RecoEgamma/PhotonIdentification/interface/CutBasedPhotonIDAlgo.h"

namespace pat {

  class PATPhotonIDProducer : public edm::EDProducer {
    
  public:
    
    explicit PATPhotonIDProducer( const edm::ParameterSet& );
    
    ~PATPhotonIDProducer() {;}
    
  private:
    
    void produce( edm::Event&, const edm::EventSetup& );
    
    CutBasedPhotonIDAlgo algo_;

    edm::InputTag phoLabel_;
    
  };
  
}

#endif // PhysicsTools_PatAlgos_PATPhotonIDProducer_h
