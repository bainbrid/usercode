#ifndef bainbrid_Test_JPTPatCorrector_h
#define bainbrid_Test_JPTPatCorrector_h

#include "bainbrid/Test/test/JPTCorrector.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

/**
   \brief Jet energy correction algorithm using PAT collections
*/
class JPTPatCorrector : public JPTCorrector {

  // ---------- Public interface ----------
  
 public: 

  /// Constructor
  JPTPatCorrector( const edm::ParameterSet& );
  
  /// Destructor
  virtual ~JPTPatCorrector();
  
  // ---------- Protected interface ----------

 protected: 

  // Some useful typedefs
  typedef edm::View<pat::Muon> PatMuons;
  typedef edm::View<pat::Electron> PatElectrons;
  typedef JPTCorrector::TrackRefs TrackRefs;
  
  /// Associates tracks to jets "on-the-fly"
  bool jtaOnTheFly( const reco::Jet&, 
		    const edm::Event&, 
		    const edm::EventSetup&,
		    jpt::JetTracks& ) const;
  
  /// Categories tracks according to particle type
  void matchTracks( const jpt::JetTracks&,
		    const edm::Event&, 
		    jpt::MatchedTracks& pions, 
		    jpt::MatchedTracks& muons, 
		    jpt::MatchedTracks& electrons ) const;

  /// Get PAT muons
  bool getMuons( const edm::Event&, edm::Handle<PatMuons>& ) const;

  /// Get PAT electrons
  bool getElectrons( const edm::Event&, edm::Handle<PatElectrons>& ) const;

  /// Matches tracks to PAT muons
  bool matchMuons( reco::TrackRefVector::const_iterator,
		   const edm::Handle<PatMuons>& ) const;
  
  /// Matches tracks to PAT electrons
  bool matchElectrons( reco::TrackRefVector::const_iterator,
		       const edm::Handle<PatElectrons>& ) const;
  
  /// Private default constructor
  JPTPatCorrector() {;}

  // ---------- Protected member data ----------
  
 protected:
  
  // Some general configuration
  bool usePat_;
  bool allowOnTheFly_;
  
  // "On-the-fly" jet-track association
  edm::InputTag tracks_;
  std::string propagator_;
  double coneSize_;
  
};

// ---------- Inline methods ----------

#endif // bainbrid_Test_JPTPatCorrector_h
