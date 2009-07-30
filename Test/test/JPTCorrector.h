#ifndef bainbrid_Test_JPTCorrector_h
#define bainbrid_Test_JPTCorrector_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "boost/range/iterator_range.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

/**
   \brief Jet energy correction algorithm using tracks
*/
class JPTCorrector : public JetCorrector {

 public: 
  
  JPTCorrector( const edm::ParameterSet& );
  virtual ~JPTCorrector();
  
  double correction( const reco::Jet&, const edm::Event&, const edm::EventSetup& ) const;
  double correction( const reco::Jet& ) const;
  double correction( const reco::Particle::LorentzVector& ) const;
  bool eventRequired() const;
  
 private: 
  
  /// Generic container class 
  class TrackResponse : public boost::numeric::ublas::matrix< std::pair<uint16_t,double> > {
  public:
    typedef std::pair<uint16_t,double> Pair;
    uint16_t nTrks( uint32_t eta_bin, uint32_t pt_bin ) const;
    double sumE( uint32_t eta_bin, uint32_t pt_bin ) const;
    double meanE( uint32_t eta_bin, uint32_t pt_bin ) const;
    void addE( uint32_t eta_bin, uint32_t pt_bin, double energy );
  private:
    bool check( uint32_t eta_bin, uint32_t pt_bin ) const;
  };
  
  /// Generic container class for response & efficiency maps
  class Map {

  public:

    Map( std::string );
    Map();
    ~Map();
    
    uint32_t nEtaBins() const;
    uint32_t nPtBins() const;
    
    double eta( uint32_t ) const;
    double pt( uint32_t ) const;
    double value( uint32_t eta_bin, uint32_t pt_bin ) const;
    
    void clear();

  public:
    
    std::vector<double> eta_;
    std::vector<double> pt_;
    typedef boost::numeric::ublas::matrix<double> matrix;
    matrix data_;
    int nEtaBins_;
    int nPtBins_;
    std::vector<double> value_;
    
  };
  
  /// Generic container of track collections
  class JetTracks {
  public:
    JetTracks();
    ~JetTracks();
    void clear();
    reco::TrackRefVector atVertex_;
    reco::TrackRefVector atCaloFace_;
  };

  /// Generic container of in-cone & out-of-cone tracks
  class Tracks {
  public:
    Tracks();
    ~Tracks();
    void clear();
    reco::TrackRefVector inVertexInCalo_;
    reco::TrackRefVector outOfVertexInCalo_;
    reco::TrackRefVector inVertexOutOfCalo_; 
  };

 private:
  
  bool associateTracksToJets( const reco::Jet&, 
			      const edm::Event&, 
			      const edm::EventSetup&,
			      JetTracks& ) const;

  typedef reco::MuonCollection Muons;
  bool matching( reco::TrackRefVector::const_iterator,
		 const edm::Handle<Muons>& ) const;

  typedef reco::GsfElectronCollection Electrons;
  typedef edm::ValueMap<float> ElectronIDs;
  bool matching( reco::TrackRefVector::const_iterator,
		 const edm::Handle<Electrons>&, 
		 const edm::Handle<ElectronIDs>& ) const;
  
  void particleId( const JetTracks&,
		   const edm::Event&, 
		   const edm::EventSetup&,
		   Tracks& pions, 
		   Tracks& muons, 
		   Tracks& electrons ) const;
  
  double correction( const reco::TrackRefVector&, 
		     TrackResponse&,
		     bool subtract_response,
		     bool add_momentum,
		     double mass = 0.14,
		     double mip = -1. ) const;
  
  double correction( TrackResponse&,
		     bool subtract_response ) const;
  
 private:
  
  // Jet-track association
  edm::InputTag jetTracksAtVertex_;
  edm::InputTag jetTracksAtCalo_;
  
  // "On-the-fly" jet-track association
  edm::InputTag tracks_;
  std::string propagator_;
  double coneSize_;

  // Muons and electrons
  edm::InputTag muons_;
  edm::InputTag electrons_; 
  edm::InputTag electronIds_;
  
  // Handle out-of-cone tracks
  bool addOutOfConeTracks_;
  
  // Filter tracks by quality
  bool useTrackQuality_;
  reco::TrackBase::TrackQuality trackQuality_;

  // Response and efficiency maps  
  Map response_;
  Map efficiency_;
  Map leakage_;
  
};

inline bool JPTCorrector::eventRequired() const { return true; }
inline uint32_t JPTCorrector::Map::nEtaBins() const { return eta_.size(); }
inline uint32_t JPTCorrector::Map::nPtBins() const { return pt_.size(); }

#endif // bainbrid_Test_JPTCorrector_h
