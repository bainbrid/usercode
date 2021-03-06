#ifndef JetMETCorrections_JetPlusTrack_ObjectMatcher_H
#define JetMETCorrections_JetPlusTrack_ObjectMatcher_H

#include "JetMETCorrections/JetPlusTrack/interface/ObjectMatcherBase.h"
#include "CLHEP/Vector/LorentzVector.h"
#include <vector>

namespace edm  { 
  class Event; 
  class EventSetup; 
  class ParameterSet; 
}

template<class GEN, class RECO>
class ObjectMatcher : public ObjectMatcherBase {
  
 public:
  
  explicit ObjectMatcher( const edm::ParameterSet& );
  
  virtual ~ObjectMatcher();
  
 private:
  
  ObjectMatcher();
  
  void gen( const edm::Event&, 
	    const edm::EventSetup&, 
	    std::vector<HepLorentzVector>& );
  
  void reco( const edm::Event&, 
	     const edm::EventSetup&, 
	     std::vector<HepLorentzVector>& );
  
};

#endif // JetMETCorrections_JetPlusTrack_ObjectMatcher_H

