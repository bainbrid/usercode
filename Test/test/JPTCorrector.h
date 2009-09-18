#ifndef bainbrid_Test_JPTCorrector_h
#define bainbrid_Test_JPTCorrector_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "boost/range/iterator_range.hpp"


// --------------------------------------------------------
// -------------------- Helper classes --------------------
// --------------------------------------------------------


namespace jpt {
  
  /// Container class for response & efficiency maps
  class Map {

  public:

    Map( std::string, bool verbose = false );
    Map();
    ~Map();
    
    uint32_t nEtaBins() const;
    uint32_t nPtBins() const;
    
    double eta( uint32_t ) const;
    double pt( uint32_t ) const;
    double value( uint32_t eta_bin, uint32_t pt_bin ) const;
    
    void clear();

  private:

    class Element {
    public:
    Element() : ieta_(0), ipt_(0), eta_(0.), pt_(0.), val_(0.) {;} 
      uint32_t ieta_;
      uint32_t ipt_;
      double eta_;
      double pt_;
      double val_;
    };
    
    typedef std::vector<double> VDouble;
    typedef std::vector<VDouble> VVDouble;
    
    std::vector<double> eta_;
    std::vector<double> pt_;
    VVDouble data_;
    
  };

  inline uint32_t Map::nEtaBins() const { return eta_.size(); }
  inline uint32_t Map::nPtBins() const { return pt_.size(); }
  
  /// Generic container class 
  class Response {

  public:

    typedef std::pair<uint16_t,double> Pair;

    uint16_t nTrks( uint32_t eta_bin, uint32_t pt_bin ) const;
    double sumE( uint32_t eta_bin, uint32_t pt_bin ) const;
    double meanE( uint32_t eta_bin, uint32_t pt_bin ) const;
    void addE( uint32_t eta_bin, uint32_t pt_bin, double energy );

    void clear();
    void resize( uint32_t eta_bins, uint32_t pt_bins, Pair value = Pair(0,0.) );

  private:

    bool check( uint32_t eta_bin, uint32_t pt_bin ) const;

    typedef std::vector<Pair> VPair;
    typedef std::vector<VPair> VVPair;
    VVPair data_;

  };
  
  /// Tracks associated to jets that are in-cone at Vertex and CaloFace
  class JetTracks {
  public:
    JetTracks();
    ~JetTracks();
    void clear();
    reco::TrackRefVector vertex_;
    reco::TrackRefVector caloFace_;
  };

  /// Particles matched to tracks that are in/in, in/out, out/in at Vertex and CaloFace
  class MatchedTracks {
  public:
    MatchedTracks();
    ~MatchedTracks();
    void clear();
    reco::TrackRefVector inVertexInCalo_;
    reco::TrackRefVector outOfVertexInCalo_;
    reco::TrackRefVector inVertexOutOfCalo_; 
  };

}


// -------------------------------------------------------
// -------------------- JPT algorithm --------------------
// -------------------------------------------------------


/**
   \brief Jet energy correction algorithm using tracks
*/
class JPTCorrector : public JetCorrector {

  // ---------- Public interface ----------
  
 public: 

  /// Constructor
  JPTCorrector( const edm::ParameterSet& );

  /// Destructor
  virtual ~JPTCorrector();

  /// Correction method
  double correction( const reco::Jet&, const edm::Event&, const edm::EventSetup& ) const;

  /// Correction method (not used)
  double correction( const reco::Jet& ) const;

  /// Correction method (not used)
  double correction( const reco::Particle::LorentzVector& ) const;

  /// Returns true
  bool eventRequired() const;

  // ---------- Extended interface ----------

  /// Can jet be JPT-corrected?
  bool canCorrect( const reco::Jet& ) const;
  
  /// Matches tracks to different particle types 
  bool matchTracks( const reco::Jet&, 
		    const edm::Event&, 
		    const edm::EventSetup&,
		    jpt::MatchedTracks& pions, 
		    jpt::MatchedTracks& muons, 
		    jpt::MatchedTracks& elecs ) const;
  
  /// Calculates corrections to be applied using pions
  double pionCorrection( const jpt::MatchedTracks& pions ) const;
  
  /// Calculates correction to be applied using muons
  double muonCorrection( const jpt::MatchedTracks& muons, bool ) const;
  
  /// Calculates correction to be applied using electrons
  double elecCorrection( const jpt::MatchedTracks& elecs ) const;

  // ---------- Protected interface ----------

 protected: 

  // Some useful typedefs
  typedef reco::MuonCollection RecoMuons;
  typedef reco::GsfElectronCollection RecoElectrons;
  typedef edm::ValueMap<float> RecoElectronIds;
  typedef reco::JetTracksAssociation::Container JetTracksAssociations;
  typedef reco::TrackRefVector TrackRefs;

  /// Associates tracks to jets
  bool jetTrackAssociation( const reco::Jet&, 
			    const edm::Event&, 
			    const edm::EventSetup&,
			    jpt::JetTracks& ) const;
  
  /// Associates tracks to jets "on-the-fly"
  virtual bool jtaOnTheFly( const reco::Jet&, 
			    const edm::Event&, 
			    const edm::EventSetup&,
			    jpt::JetTracks& ) const;
  
  /// Matches tracks to different particle types 
  virtual void matchTracks( const jpt::JetTracks&,
			    const edm::Event&, 
			    jpt::MatchedTracks& pions, 
			    jpt::MatchedTracks& muons, 
			    jpt::MatchedTracks& elecs ) const;

  /// Calculates individual pion corrections
  double pionCorrection( const TrackRefs& pions, 
			 jpt::Response& response,
			 bool in_cone_at_vertex,
			 bool in_cone_at_calo_face ) const; 

  /// Calculates individual muons corrections
  double muonCorrection( const TrackRefs& muons, 
			 bool in_cone_at_vertex,
			 bool in_cone_at_calo_face ) const;

  /// Calculates individual electron corrections
  double elecCorrection( const TrackRefs& elecs, 
			 bool in_cone_at_vertex,
			 bool in_cone_at_calo_face ) const;
  
  /// Generic method to calculates correction to be applied
  double correction( const TrackRefs&, 
		     jpt::Response&,
		     bool in_cone_at_vertex,
		     bool in_cone_at_calo_face,
		     double mass = 0.14,
		     double mip = -1. ) const;
  
  /// Correction to be applied using tracking efficiency 
  double pionEfficiency( jpt::Response&,
			 bool in_cone_at_calo_face ) const;

  /// Check scale is not negative
  double checkScale( double scale ) const;
  
  /// Get RECO muons
  bool getMuons( const edm::Event&, edm::Handle<RecoMuons>& ) const;

  /// Get RECO electrons
  bool getElectrons( const edm::Event&, 
		     edm::Handle<RecoElectrons>&, 
		     edm::Handle<RecoElectronIds>& ) const;
  
  /// Matches tracks to RECO muons
  bool matchMuons( TrackRefs::const_iterator,
		   const edm::Handle<RecoMuons>& ) const;
  
  /// Matches tracks to RECO electrons
  bool matchElectrons( TrackRefs::const_iterator,
		       const edm::Handle<RecoElectrons>&, 
		       const edm::Handle<RecoElectronIds>& ) const;
  
  /// Check on track quality
  bool failTrackQuality( TrackRefs::const_iterator ) const;

  /// Find track in JetTracks collection
  bool findTrack( const jpt::JetTracks&, 
		  TrackRefs::const_iterator,
		  TrackRefs::iterator& ) const;

  /// Find track in MatchedTracks collections
  bool findTrack( const jpt::MatchedTracks& pions, 
		  const jpt::MatchedTracks& muons,
		  const jpt::MatchedTracks& electrons,
		  TrackRefs::const_iterator ) const;

  /// Determines if any tracks in cone at CaloFace
  bool tracksInCalo( const jpt::MatchedTracks& pions, 
		     const jpt::MatchedTracks& muons,
		     const jpt::MatchedTracks& elecs ) const;
  
  /// Rebuild jet-track association 
  void rebuildJta( const reco::Jet&, 
		   const JetTracksAssociations&, 
		   TrackRefs& included,
		   TrackRefs& excluded ) const;
  
  /// Exclude jet-track association 
  void excludeJta( const reco::Jet&, 
		   const JetTracksAssociations&, 
		   TrackRefs& included,
		   const TrackRefs& excluded ) const;

  // Methods to access maps
  const jpt::Map& response() const;
  const jpt::Map& efficiency() const;
  const jpt::Map& leakage() const;

  /// Default constructor
  JPTCorrector() {;}

  // ---------- Protected member data ----------

 protected:
  
  // Some general configuration
  bool verbose_;
  bool useInConeTracks_;
  bool useOutOfConeTracks_;
  bool useOutOfVertexTracks_;
  bool usePions_;
  bool useEff_;
  bool useMuons_;
  bool useElecs_;
  bool useTrackQuality_;
  
  // Jet-track association
  edm::InputTag jetTracksAtVertex_;
  edm::InputTag jetTracksAtCalo_;
  int jetSplitMerge_;

  // Muons and electrons
  edm::InputTag muons_;
  edm::InputTag electrons_; 
  edm::InputTag electronIds_;
  
  // Filter tracks by quality
  reco::TrackBase::TrackQuality trackQuality_;

  // Response and efficiency maps  
  const jpt::Map* response_;
  const jpt::Map* efficiency_;
  const jpt::Map* leakage_;
  
};

// ---------- Inline methods ----------

inline bool JPTCorrector::eventRequired() const { return true; }
inline bool JPTCorrector::canCorrect( const reco::Jet& jet ) const { return ( fabs( jet.eta() ) <= 2.1 ); }

inline const jpt::Map& JPTCorrector::response() const { return *response_; }
inline const jpt::Map& JPTCorrector::efficiency() const { return *efficiency_; }
inline const jpt::Map& JPTCorrector::leakage() const { return *leakage_; }

inline double JPTCorrector::pionCorrection( const TrackRefs& pions, 
					    jpt::Response& response,
					    bool in_cone_at_vertex,
					    bool in_cone_at_calo_face ) const {
  return correction( pions, response, in_cone_at_vertex, in_cone_at_calo_face );
}

inline double JPTCorrector::muonCorrection( const TrackRefs& muons, 
					    bool in_cone_at_vertex,
					    bool in_cone_at_calo_face ) const {
  jpt::Response response;
  return correction( muons, response, in_cone_at_vertex, in_cone_at_calo_face, 0.105, 2. );
} 

inline double JPTCorrector::elecCorrection( const TrackRefs& elecs, 
					    bool in_cone_at_vertex,
					    bool in_cone_at_calo_face ) const {
  jpt::Response response;
  return correction( elecs, response, in_cone_at_vertex, in_cone_at_calo_face, 0.000511, 0. );
} 

inline double JPTCorrector::checkScale( double scale ) const {
  if ( scale < 0. ) { return 1.; } 
  else { return scale; }
}

#endif // bainbrid_Test_JPTCorrector_h
