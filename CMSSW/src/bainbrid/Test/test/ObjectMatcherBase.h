#ifndef JetMETCorrections_JetPlusTrack_ObjectMatcherBase_H
#define JetMETCorrections_JetPlusTrack_ObjectMatcherBase_H

#include "bainbrid/Test/test/ObjectTags.h"
#include "bainbrid/Test/test/LorentzVectorPair.h"
#include <vector>
#include <string>

#include <TFile.h>
#include <TTree.h>

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
  
  /// Set verbosity
  void verbose( bool );
  
 protected:
  
  /// Constructor
  explicit ObjectMatcherBase( const edm::ParameterSet& );
  
  /// Default constructor
  ObjectMatcherBase();
  
  /// Extracts 4-vectors from generator objects
  virtual void gen( const edm::Event&, 
		    const edm::EventSetup&, 
		    std::vector<math::XYZTLorentzVector>& ) = 0;
  
  /// Extracts 4-vectors from reco objects
  virtual void reco( const edm::Event&, 
		     const edm::EventSetup&, 
		     std::vector<math::XYZTLorentzVector>& ) = 0;
  
  /// Returns string to identify specialised class  
  std::string id();
  
 protected:

  TFile* hOutputFile;
  TTree* t1;
  double  EtaGen1, PhiGen1, EtaRaw1, PhiRaw1, EtGen1, EtRaw1, EtMCJ1, EtZSP1, EtJPT1, DRMAXgjet1;
  double  EtaGen2, PhiGen2, EtaRaw2, PhiRaw2, EtGen2, EtRaw2, EtMCJ2, EtZSP2, EtJPT2, DRMAXgjet2;
  
  // "Old skool"
  void analyze( const edm::Event&, 
		const std::vector<math::XYZTLorentzVector>&,
		const std::vector<math::XYZTLorentzVector>& );

  // "Old skool"
  void analyze( const edm::Event&, 
		const std::vector<LorentzVectorPair>& );
  
  /// Tags used to identify objects
  ObjectTags tags_;

  /// Verbosity flag
  bool verbose_;
  
};

inline void ObjectMatcherBase::verbose( bool v ) { verbose_ = v; }

#endif // JetMETCorrections_JetPlusTrack_ObjectMatcherBase_H
