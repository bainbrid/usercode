#ifndef JetMETCorrections_JetPlusTrack_EnergyScaleAnalyzer_H
#define JetMETCorrections_JetPlusTrack_EnergyScaleAnalyzer_H

#include "FWCore/Framework/interface/EDAnalyzer.h"

namespace edm {
  class ActivityRegistry;
  class Event;
  class EventSetup;
}
class ObjectMatcherBase;

/**
   @brief Service 
*/
class EnergyScaleAnalyzer : public edm::EDAnalyzer {
  
 public:
  
  /// Service constructor called by framework
  explicit EnergyScaleAnalyzer( const edm::ParameterSet& );
  
  /// Destructor
  virtual ~EnergyScaleAnalyzer();
  
  /// Event loop
  void analyze( const edm::Event&, const edm::EventSetup& );
  
 private:

  ObjectMatcherBase* matcher_;
  
};

#endif // JetMETCorrections_JetPlusTrack_EnergyScaleAnalyzer_H
