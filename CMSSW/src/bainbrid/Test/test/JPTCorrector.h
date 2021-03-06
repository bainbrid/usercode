#ifndef bainbrid_Test_JPTCorrector_h
#define bainbrid_Test_JPTCorrector_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "boost/range/iterator_range.hpp"
#include <sstream>
#include <string>

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

    uint32_t etaBin( double eta ) const;
    uint32_t ptBin( double pt ) const;
    
    double value( uint32_t eta_bin, uint32_t pt_bin ) const;

    double binCenterEta( uint32_t ) const;
    double binCenterPt( uint32_t ) const;
    
    void clear();
    void print( std::stringstream& ss ) const;

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
  class Efficiency {

  public:

    Efficiency( const jpt::Map& response,
		const jpt::Map& efficiency,
		const jpt::Map& leakage );
      
    Efficiency();

    typedef std::pair<uint16_t,double> Pair;
    
    uint16_t nTrks( uint32_t eta_bin, uint32_t pt_bin ) const;
    
    double inConeCorr( uint32_t eta_bin, uint32_t pt_bin ) const;
    double outOfConeCorr( uint32_t eta_bin, uint32_t pt_bin ) const;
    
    uint32_t nEtaBins() const;
    uint32_t nPtBins() const;

    uint32_t size() const;
    bool empty() const;
    
    void addE( uint32_t eta_bin, uint32_t pt_bin, double energy );
    void reset();
    
    void print() const;
    
  private:
    
    double sumE( uint32_t eta_bin, uint32_t pt_bin ) const;
    double meanE( uint32_t eta_bin, uint32_t pt_bin ) const;

    bool check( uint32_t eta_bin, uint32_t pt_bin, std::string name = "check" ) const;

    typedef std::vector<Pair> VPair;
    typedef std::vector<VPair> VVPair;
    VVPair data_;

    jpt::Map response_;
    jpt::Map efficiency_;
    jpt::Map leakage_;
    
  };
  
  inline uint32_t Efficiency::nEtaBins() const { return response_.nEtaBins(); }
  inline uint32_t Efficiency::nPtBins() const { return response_.nPtBins(); }
  inline uint32_t Efficiency::size() const { return data_.size(); }
  inline bool Efficiency::empty() const { return data_.empty(); }

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

  // Typedefs for 4-momentum
  typedef JetCorrector::LorentzVector P4;
  typedef math::PtEtaPhiELorentzVectorD PtEtaPhiE;
  typedef math::PtEtaPhiMLorentzVectorD PtEtaPhiM;
  
  /// Vectorial correction method (corrected 4-momentum passed by reference)
  double correction( const reco::Jet&, const edm::Event&, const edm::EventSetup&, P4& ) const;
  
  /// Scalar correction method
  double correction( const reco::Jet&, const edm::Event&, const edm::EventSetup& ) const;
  
  /// Correction method (not used)
  double correction( const reco::Jet& ) const;

  /// Correction method (not used)
  double correction( const P4& ) const;
  
  /// Returns true
  bool eventRequired() const;
  
  /// Returns value of configurable
  bool vectorialCorrection() const;
  
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
  P4 pionCorrection( const P4& jet, const jpt::MatchedTracks& pions ) const;
  
  /// Calculates correction to be applied using muons
  P4 muonCorrection( const P4& jet, const jpt::MatchedTracks& muons, bool ) const;
  
  /// Calculates correction to be applied using electrons
  P4 elecCorrection( const P4& jet, const jpt::MatchedTracks& elecs ) const;
  
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
  P4 pionCorrection( const P4& jet, 
		     const TrackRefs& pions, 
		     bool in_cone_at_vertex,
		     bool in_cone_at_calo_face ) const; 

  /// Calculates individual muons corrections
  P4 muonCorrection( const P4& jet, 
		     const TrackRefs& muons, 
		     bool in_cone_at_vertex,
		     bool in_cone_at_calo_face ) const;
  
  /// Calculates individual electron corrections
  P4 elecCorrection( const P4& jet, 
		     const TrackRefs& elecs, 
		     bool in_cone_at_vertex,
		     bool in_cone_at_calo_face ) const;

  /// Calculates vectorial correction using total track 3-momentum
  P4 jetDirFromTracks( const P4& jet, 
		       const jpt::MatchedTracks& pions,
		       const jpt::MatchedTracks& muons,
		       const jpt::MatchedTracks& elecs ) const;
  
  /// Generic method to calculates 4-momentum correction to be applied
  P4 calculateCorr( const P4& jet, 
		    const TrackRefs&, 
		    bool in_cone_at_vertex,
		    bool in_cone_at_calo_face,
		    double mass,
		    bool is_mip,
		    double mip ) const;
  
  /// Correction to be applied using tracking efficiency 
  P4 pionEfficiency( const P4& jet, bool in_cone_at_calo_face ) const;
  
  /// Check corrected 4-momentum does not give negative scale
  double checkScale( const P4& jet, P4& corrected ) const;
  
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

  /// Default constructor
  JPTCorrector() {;}

  // ---------- Protected member data ----------

 protected:
  
  // Some general configuration
  bool verbose_;
  bool vectorial_;
  bool vecTracks_;
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
  const jpt::Map response_;
  const jpt::Map efficiency_;
  const jpt::Map leakage_;
  mutable jpt::Efficiency eff_;

  // Mass    
  double pionMass_;
  double muonMass_;
  double elecMass_;
  
};

// ---------- Inline methods ----------

inline double JPTCorrector::correction( const reco::Jet& fJet,
					const edm::Event& event,
					const edm::EventSetup& setup ) const {
  P4 not_used_for_scalar_correction;
  return correction( fJet, event, setup, not_used_for_scalar_correction );
}

inline bool JPTCorrector::eventRequired() const { return true; }
inline bool JPTCorrector::vectorialCorrection() const { return vectorial_; }
inline bool JPTCorrector::canCorrect( const reco::Jet& jet ) const { return ( fabs( jet.eta() ) <= 2.1 ); }

inline JPTCorrector::P4 JPTCorrector::pionCorrection( const P4& jet, 
						      const TrackRefs& pions, 
						      bool in_cone_at_vertex,
						      bool in_cone_at_calo_face ) const {
  return calculateCorr( jet, pions, in_cone_at_vertex, in_cone_at_calo_face, pionMass_, false, -1. );
}

inline JPTCorrector::P4 JPTCorrector::muonCorrection( const P4& jet, 
						      const TrackRefs& muons, 
						      bool in_cone_at_vertex,
						      bool in_cone_at_calo_face ) const {
  return calculateCorr( jet, muons, in_cone_at_vertex, in_cone_at_calo_face, muonMass_, true, 2. );
} 

inline JPTCorrector::P4 JPTCorrector::elecCorrection( const P4& jet, 
						      const TrackRefs& elecs, 
						      bool in_cone_at_vertex,
						      bool in_cone_at_calo_face ) const {
  return calculateCorr( jet, elecs, in_cone_at_vertex, in_cone_at_calo_face, elecMass_, true, 0. ); 
} 

inline double JPTCorrector::checkScale( const P4& jet, P4& corrected ) const {
  if ( jet.energy() > 0. && ( corrected.energy() / jet.energy() ) < 0. ) { 
    corrected = jet; 
  }
  return corrected.energy() / jet.energy();
}

#endif // bainbrid_Test_JPTCorrector_h
