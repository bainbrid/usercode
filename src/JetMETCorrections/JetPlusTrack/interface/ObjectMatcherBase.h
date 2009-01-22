#ifndef JetMETCorrections_JetPlusTrack_ObjectMatcherBase_H
#define JetMETCorrections_JetPlusTrack_ObjectMatcherBase_H

#include "JetMETCorrections/JetPlusTrack/interface/ObjectTags.h"
#include "CLHEP/Vector/LorentzVector.h"
#include <vector>
#include <string>

namespace edm  { 
  class Event; 
  class EventSetup; 
  class ParameterSet; 
}

class ObjectMatcherBase {
  
 public:
  
  /// Destructor
  virtual ~ObjectMatcherBase();
  
  /// Creates instance of derived class 
  static ObjectMatcherBase* Instance( const edm::ParameterSet& );
  
  /// Public analyze method 
  void analyze( const edm::Event&, const edm::EventSetup& );
  
 protected:
  
  /// Constructor
  explicit ObjectMatcherBase( const edm::ParameterSet& );
  
  /// Default constructor
  ObjectMatcherBase();
  
  /// Extracts 4-vectors from generator objects
  virtual void gen( const edm::Event&, 
		    const edm::EventSetup&, 
		    std::vector<HepLorentzVector>& ) = 0;
  
  /// Extracts 4-vectors from reco objects
  virtual void reco( const edm::Event&, 
		     const edm::EventSetup&, 
		     std::vector<HepLorentzVector>& ) = 0;

  /// Returns string to identify specialised class  
  std::string id();
  
 protected:
  
  /// Tags used to identify objects
  ObjectTags tags_;
  
};

#endif // JetMETCorrections_JetPlusTrack_ObjectMatcherBase_H
