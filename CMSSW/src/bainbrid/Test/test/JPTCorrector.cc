#include "bainbrid/Test/test/JPTCorrector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <fstream>
#include <vector>

using namespace std;
using namespace jpt;

// -----------------------------------------------------------------------------
//
JPTCorrector::JPTCorrector( const edm::ParameterSet& pset ) 
  : verbose_( pset.getParameter<bool>("Verbose") ),
    vectorial_( pset.getParameter<bool>("VectorialCorrection") ),
    vecTracks_( pset.getParameter<bool>("JetDirFromTracks") ),
    useInConeTracks_( pset.getParameter<bool>("UseInConeTracks") ),
    useOutOfConeTracks_( pset.getParameter<bool>("UseOutOfConeTracks") ),
    useOutOfVertexTracks_( pset.getParameter<bool>("UseOutOfVertexTracks") ),
    usePions_( pset.getParameter<bool>("UsePions") ),
    useEff_( pset.getParameter<bool>("UseEfficiency") ),
    useMuons_( pset.getParameter<bool>("UseMuons") ),
    useElecs_( pset.getParameter<bool>("UseElectrons") ),
    useTrackQuality_( pset.getParameter<bool>("UseTrackQuality") ),
    jetTracksAtVertex_( pset.getParameter<edm::InputTag>("JetTracksAssociationAtVertex") ),
    jetTracksAtCalo_( pset.getParameter<edm::InputTag>("JetTracksAssociationAtCaloFace") ),
    jetSplitMerge_( pset.getParameter<int>("JetSplitMerge") ),
    muons_( pset.getParameter<edm::InputTag>("Muons") ),
    electrons_( pset.getParameter<edm::InputTag>("Electrons") ),
    electronIds_( pset.getParameter<edm::InputTag>("ElectronIds") ),
    trackQuality_( reco::TrackBase::qualityByName( pset.getParameter<std::string>("TrackQuality") ) ),
    response_( Map( pset.getParameter<std::string>("ResponseMap"), verbose_ ) ),
    efficiency_( Map( pset.getParameter<std::string>("EfficiencyMap"), verbose_ ) ),
    leakage_( Map( pset.getParameter<std::string>("LeakageMap"), verbose_ ) ),
    eff_( response_, efficiency_, leakage_ ),
    pionMass_(0.140),
    muonMass_(0.105),
    elecMass_(0.000511)
{
  
  if ( !useInConeTracks_ || 
       !useOutOfConeTracks_ ||
       !useOutOfVertexTracks_ ) {
    std::stringstream ss;
    ss << "[JPTCorrector::" << __func__ << "]"
       << " You are using JPT algorithm in a non-standard way!" << std::endl
       << " UseInConeTracks      : " << ( useInConeTracks_ ? "true" : "false" ) << std::endl
       << " UseOutOfConeTracks   : " << ( useOutOfConeTracks_ ? "true" : "false" ) << std::endl
       << " UseOutOfVertexTracks : " << ( useOutOfVertexTracks_ ? "true" : "false" );
    edm::LogWarning("JPTCorrector") << ss.str();
  }

  LogDebug("TESTTEST") << ( vectorial_ ? "true" : "false" );

}

// -----------------------------------------------------------------------------
//
JPTCorrector::~JPTCorrector() {;}

// -----------------------------------------------------------------------------
//
double JPTCorrector::correction( const reco::Jet& fJet,
				 const edm::Event& event,
				 const edm::EventSetup& setup,
				 P4& corrected ) const 
{
  
  // Corrected 4-momentum for jet
  corrected = fJet.p4();
  
  // Check if jet can be JPT-corrected
  if ( !canCorrect(fJet) ) { return 1.; }
  
  // Match tracks to different particle types
  MatchedTracks pions;
  MatchedTracks muons;
  MatchedTracks elecs;
  bool ok = matchTracks( fJet, event, setup, pions, muons, elecs );
  if ( !ok ) { return 1.; }
  
  // Debug
  if ( verbose_ ) {
    edm::LogInfo("JPTCorrector")
      << "[JPTCorrector::" << __func__ << "]"
      << " Applying JPT corrections...";
  }

  // Pion corrections (both scalar and vectorial)
  if ( usePions_ ) { corrected += pionCorrection( fJet.p4(), pions ); }
  
  // Muon corrections (both scalar and vectorial)
  if ( useMuons_ ) { corrected += muonCorrection( fJet.p4(), muons, !pions.inVertexOutOfCalo_.empty() ); }
  
  // Electrons corrections (both scalar and vectorial)
  if ( useElecs_ ) { corrected += elecCorrection( fJet.p4(), elecs ); }

  // Define jet direction using total 3-momentum of tracks (overrides above)
  if ( vecTracks_ ) { corrected = jetDirFromTracks( corrected, pions, muons, elecs ); }
  
  // Check if corrected 4-momentum gives negative scale
  double scale = checkScale( fJet.p4(), corrected );
  
  // Debug
  if ( verbose_ ) {
    std::stringstream ss;
    ss << "Total correction:" << std::endl
       << " Uncorrected energy : " << fJet.energy() << std::endl
       << " Corrected energy   : " << corrected.energy() << std::endl
       << " Scalar correction  : " << scale;
    edm::LogVerbatim("JPTCorrector") << ss.str();
  }
  
//   LogTrace("test") << " mScale= " << scale
// 		   << " NewResponse " << corrected.energy() 
// 		   << " Jet energy " << fJet.energy()
// 		   << " event " << event.id().event();
  
  // Return energy correction
  return scale;
  
}

// -----------------------------------------------------------------------------
//
double JPTCorrector::correction( const reco::Jet& jet ) const {
  edm::LogError("JPTCorrector")
    << "JPTCorrector can be run on entire event only";
  return 1.;
}

// -----------------------------------------------------------------------------
//
double JPTCorrector::correction( const reco::Particle::LorentzVector& jet ) const {
  edm::LogError("JPTCorrector")
    << "JPTCorrector can be run on entire event only";
  return 1.;
}

// -----------------------------------------------------------------------------
//
bool JPTCorrector::matchTracks( const reco::Jet& fJet, 
				const edm::Event& event, 
				const edm::EventSetup& setup, //@@ required by method in derived class
				jpt::MatchedTracks& pions, 
				jpt::MatchedTracks& muons, 
				jpt::MatchedTracks& elecs ) const {
  
  // Associate tracks to jet at both the Vertex and CaloFace
  JetTracks jet_tracks;
  bool ok = jetTrackAssociation( fJet, event, setup, jet_tracks ); 
  if ( !ok ) { return false; }
  
  // Track collections propagated to Vertex and CaloFace for "pions", muons and electrons
  matchTracks( jet_tracks, event, pions, muons, elecs );

  // Debug
  if ( verbose_ ) {
    std::stringstream ss;
    ss << "Number of tracks:" << std::endl
       << " In-cone at Vertex   : " << jet_tracks.vertex_.size() << std::endl
       << " In-cone at CaloFace : " << jet_tracks.caloFace_.size();
    edm::LogVerbatim("JPTCorrector") << ss.str();
  }
  
  return true;
  
}

// -----------------------------------------------------------------------------
//
bool JPTCorrector::jetTrackAssociation( const reco::Jet& fJet,
					const edm::Event& event, 
					const edm::EventSetup& setup,
					JetTracks& trks ) const {
  
  // Some init
  trks.clear();
  
  // Check whether to retrieve JTA object from Event or construct "on-the-fly"
  if ( !jetTracksAtVertex_.label().empty() && 
       !jetTracksAtCalo_.label().empty() ) { 
    
    // Get Jet-track association at Vertex
    edm::Handle<reco::JetTracksAssociation::Container> jetTracksAtVertex;
    event.getByLabel( jetTracksAtVertex_, jetTracksAtVertex ); 
    if ( !jetTracksAtVertex.isValid() || jetTracksAtVertex.failedToGet() ) {
      if ( verbose_ && edm::isDebugEnabled() ) {
	edm::LogWarning("JPTCorrector")
	  << "[JPTCorrector::" << __func__ << "]"
	  << " Invalid handle to reco::JetTracksAssociation::Container (for Vertex)"
	  << " with InputTag (label:instance:process) \"" 
	  << jetTracksAtVertex_.label() << ":"
	  << jetTracksAtVertex_.instance() << ":"
	  << jetTracksAtVertex_.process() << "\"" << std::endl
	  << " Attempting to use JTA \"on-the-fly\" mode...";
      }
      return jtaOnTheFly( fJet, event, setup, trks );
    }
    
    // Retrieve jet-tracks association for given jet
    const reco::JetTracksAssociation::Container jtV = *( jetTracksAtVertex.product() );
    TrackRefs excluded; 
    if ( jetSplitMerge_ < 0 ) { trks.vertex_ = reco::JetTracksAssociation::getValue( jtV, fJet ); }
    else { rebuildJta( fJet, jtV, trks.vertex_, excluded ); }
    
    // Check if any tracks are associated to jet at vertex
    if ( trks.vertex_.empty() ) { return false; }

    // Get Jet-track association at Calo
    edm::Handle<reco::JetTracksAssociation::Container> jetTracksAtCalo;
    event.getByLabel( jetTracksAtCalo_, jetTracksAtCalo ); 
    if ( !jetTracksAtCalo.isValid() || jetTracksAtCalo.failedToGet() ) {
      if ( verbose_ && edm::isDebugEnabled() ) {
	edm::LogWarning("JPTCorrector")
	  << "[JPTCorrector::" << __func__ << "]"
	  << " Invalid handle to reco::JetTracksAssociation::Container (for CaloFace)"
	  << " with InputTag (label:instance:process) \"" 
	  << jetTracksAtCalo_.label() << ":"
	  << jetTracksAtCalo_.instance() << ":"
	  << jetTracksAtCalo_.process() << "\"" << std::endl
	  << " Attempting to use JTA \"on-the-fly\" mode...";
      }
      return jtaOnTheFly( fJet, event, setup, trks );
    }
    
    // Retrieve jet-tracks association for given jet
    const reco::JetTracksAssociation::Container jtC = *( jetTracksAtCalo.product() );
    if ( jetSplitMerge_ < 0 ) { trks.caloFace_ = reco::JetTracksAssociation::getValue( jtC, fJet ); }
    else { excludeJta( fJet, jtC, trks.caloFace_, excluded ); }
    
    // Successful
    return true;
    
  } else { return jtaOnTheFly( fJet, event, setup, trks ); }
  
}

// -----------------------------------------------------------------------------
//
bool JPTCorrector::jtaOnTheFly( const reco::Jet& fJet,
				const edm::Event& event, 
				const edm::EventSetup& setup,
				JetTracks& trks ) const {
  edm::LogWarning("JPTCorrector") 
    << "[JPTCorrector::" << __func__ << "]"
    << " \"On-the-fly\" mode not available in this version of JPT!";
  return false;
}

// -----------------------------------------------------------------------------
//
bool JPTCorrector::getMuons( const edm::Event& event, edm::Handle<RecoMuons>& reco_muons ) const {
  event.getByLabel( muons_, reco_muons ); 
  if ( !reco_muons.isValid() || reco_muons.failedToGet() ) {
    edm::LogError("JPTCorrector")
      << "[JPTCorrector::" << __func__ << "]"
      << " Invalid handle to reco::Muon collection"
      << " with InputTag (label:instance:process) \"" 
      << muons_.label() << ":"
      << muons_.instance() << ":"
      << muons_.process() << "\"";
    return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
//
void JPTCorrector::matchTracks( const JetTracks& jet_tracks, 
				const edm::Event& event, 
				MatchedTracks& pions, 
				MatchedTracks& muons,
				MatchedTracks& elecs ) const { 
  
  // Some init  
  pions.clear(); 
  muons.clear(); 
  elecs.clear(); 

  // Get RECO muons
  edm::Handle<RecoMuons> reco_muons;
  bool found_reco_muons = true;
  if ( useMuons_ ) { getMuons( event, reco_muons ); }
  
  // Get RECO electrons and their ids
  edm::Handle<RecoElectrons> reco_elecs;
  edm::Handle<RecoElectronIds> reco_elec_ids;
  bool found_reco_elecs = true;
  if ( useElecs_ ) { getElectrons( event, reco_elecs, reco_elec_ids ); }

  // Check RECO products found
  if ( !found_reco_muons || !found_reco_elecs ) {
    edm::LogError("JPTCorrector")
      << "[JPTCorrector::" << __func__ << "]"
      << " Unable to access RECO collections for muons and electrons";
    return;
  }
  
  // Identify pions/muons/electrons that are "in/in" and "in/out"
  {
    TrackRefs::const_iterator itrk = jet_tracks.vertex_.begin();
    TrackRefs::const_iterator jtrk = jet_tracks.vertex_.end();
    for ( ; itrk != jtrk; ++itrk ) {

      if ( failTrackQuality(itrk) ) { continue; }
      
      TrackRefs::iterator it = jet_tracks.caloFace_.end();
      bool found = findTrack( jet_tracks, itrk, it );

      bool is_muon = useMuons_ && matchMuons( itrk, reco_muons );
      bool is_ele  = useElecs_ && matchElectrons( itrk, reco_elecs, reco_elec_ids );

      if ( found ) { 
	if ( is_muon )     { muons.inVertexInCalo_.push_back(*it); }
	else if ( is_ele ) { elecs.inVertexInCalo_.push_back(*it); } 
	else               { pions.inVertexInCalo_.push_back(*it); } 
      } else { 
	if ( is_muon )     { muons.inVertexOutOfCalo_.push_back(*itrk); }
	//else if ( is_ele ) { elecs.inVertexOutOfCalo_.push_back(*itrk); } //@@ bug?  
	else               { pions.inVertexOutOfCalo_.push_back(*itrk); }
      } 
    } 
  }
  
  // Identify pions/muons/electrons that are "out/in"
  {
    TrackRefs::iterator itrk = jet_tracks.caloFace_.begin(); 
    TrackRefs::iterator jtrk = jet_tracks.caloFace_.end(); 
    for ( ; itrk != jtrk; ++itrk ) {
      
      if ( failTrackQuality(itrk) ) { continue; }
      
      if ( !tracksInCalo( pions, muons, elecs ) ) { continue; }

      bool found = findTrack( pions, muons, elecs, itrk );
      
      if ( !found ) {
	
	bool is_muon = useMuons_ && matchMuons( itrk, reco_muons );
	bool is_ele  = false; //@@ bug? useElecs_ && matchElectrons( itrk, reco_elecs, reco_elec_ids );
	
	if ( is_muon )     { muons.outOfVertexInCalo_.push_back(*itrk); } 
	else if ( is_ele ) { elecs.outOfVertexInCalo_.push_back(*itrk); } //@@ bug?
	else               { pions.outOfVertexInCalo_.push_back(*itrk); }
	
      }
    } 
  }
  
  if ( verbose_ && edm::isDebugEnabled() ) {
    std::stringstream ss;
    ss << "[JPTCorrector::" << __func__ << "] Number of tracks:" << std::endl 
       << " In-cone at Vertex and in-cone at CaloFace:" << std::endl  
       << "  Pions      : " << pions.inVertexInCalo_.size() << std::endl
       << "  Muons      : " << muons.inVertexInCalo_.size() << std::endl
       << "  Electrons  : " << elecs.inVertexInCalo_.size() << std::endl
       << " In-cone at Vertex and out-of-cone at CaloFace:" << std::endl  
       << "  Pions      : " << pions.inVertexOutOfCalo_.size() << std::endl
       << "  Muons      : " << muons.inVertexOutOfCalo_.size() << std::endl
       << "  Electrons  : " << elecs.inVertexOutOfCalo_.size() << std::endl
       << " Out-of-cone at Vertex and in-cone at CaloFace:" << std::endl  
       << "  Pions      : " << pions.outOfVertexInCalo_.size() << std::endl
       << "  Muons      : " << muons.outOfVertexInCalo_.size() << std::endl
       << "  Electrons  : " << elecs.outOfVertexInCalo_.size() << std::endl;
    LogTrace("JPTCorrector") << ss.str();
  }
  
}

// -----------------------------------------------------------------------------
//
bool JPTCorrector::getElectrons( const edm::Event& event, 
				 edm::Handle<RecoElectrons>& reco_elecs,
				 edm::Handle<RecoElectronIds>& reco_elec_ids ) const {
  event.getByLabel( electrons_, reco_elecs ); 
  if ( !reco_elecs.isValid() || reco_elecs.failedToGet() ) {
    edm::LogError("JPTCorrector")
      << "[JPTCorrector::" << __func__ << "]"
      << " Invalid handle to reco::GsfElectron collection"
      << " with InputTag (label:instance:process) \"" 
      << electrons_.label() << ":"
      << electrons_.instance() << ":"
      << electrons_.process() << "\"";
    return false;
  }
  event.getByLabel( electronIds_, reco_elec_ids ); 
  if ( !reco_elec_ids.isValid() || reco_elec_ids.failedToGet() ) {
    edm::LogError("JPTCorrector")
      << "[JPTCorrector::" << __func__ << "]"
      << " Invalid handle to reco::GsfElectron collection"
      << " with InputTag (label:instance:process) \"" 
      << electronIds_.label() << ":"
      << electronIds_.instance() << ":"
      << electronIds_.process() << "\"";
    return false;
  }
  return true;
} 

// -----------------------------------------------------------------------------
//
bool JPTCorrector::failTrackQuality( TrackRefs::const_iterator itrk ) const { 
  if ( useTrackQuality_ && !(*itrk)->quality(trackQuality_) ) { return true; }
  else { return false; }
}

// -----------------------------------------------------------------------------
//
bool JPTCorrector::findTrack( const JetTracks& jet_tracks,
			      TrackRefs::const_iterator itrk,
			      TrackRefs::iterator& it ) const { 
  it = find( jet_tracks.caloFace_.begin(),
	     jet_tracks.caloFace_.end(),
	     *itrk );
  if ( it != jet_tracks.caloFace_.end() ) { return true; }
  else { return false; }
}

// -----------------------------------------------------------------------------
//
bool JPTCorrector::findTrack( const MatchedTracks& pions, 
			      const MatchedTracks& muons,
			      const MatchedTracks& elecs,
			      TrackRefs::const_iterator itrk ) const { 
  TrackRefs::iterator ip = find( pions.inVertexInCalo_.begin(),
				 pions.inVertexInCalo_.end(),
				 *itrk );
  TrackRefs::iterator im = find( muons.inVertexInCalo_.begin(),
				 muons.inVertexInCalo_.end(),
				 *itrk );
  TrackRefs::iterator ie = find( elecs.inVertexInCalo_.begin(),
				 elecs.inVertexInCalo_.end(),
				 *itrk );
  if ( ip == pions.inVertexInCalo_.end() &&
       im == muons.inVertexInCalo_.end() /* && */ ) { 
    /* ie == elecs.inVertexInCalo_.end() ) { */ return false; } //@@ bug?
  else { return true; }
}

// -----------------------------------------------------------------------------
//
bool JPTCorrector::tracksInCalo( const MatchedTracks& pions, 
				 const MatchedTracks& muons,
				 const MatchedTracks& elecs ) const { 
  if ( !pions.inVertexInCalo_.empty() /* || */ ) {
    /* !muons.inVertexInCalo_.empty() || */
    /* !elecs.inVertexInCalo_.empty() ) { */ return true; } //@@ bug?
  else { return false; }
}

// -----------------------------------------------------------------------------
//
JPTCorrector::P4 JPTCorrector::pionCorrection( const P4& jet,
					       const MatchedTracks& pions ) const {
  
  P4 corr_pions;
  
  P4 corr_pions_in_cone;
  P4 corr_pions_out_of_cone;
  P4 corr_pions_out_of_vertex;

  P4 corr_pions_eff_in_cone;
  P4 corr_pions_eff_out_of_cone;

  // Corrections with pions 

  if ( useInConeTracks_ ) { 
    corr_pions_in_cone = pionCorrection( jet, pions.inVertexInCalo_, true, true ); 
    corr_pions += corr_pions_in_cone;
    if ( useEff_ ) {
      corr_pions_eff_in_cone = pionEfficiency( jet, true );
      corr_pions += corr_pions_eff_in_cone;
    }
  }
  
  if ( useOutOfConeTracks_ ) {
    corr_pions_out_of_cone = pionCorrection( jet, pions.inVertexOutOfCalo_, true, false );
    corr_pions += corr_pions_out_of_cone;
    if ( useEff_ ) {
      corr_pions_eff_out_of_cone = pionEfficiency( jet, false );
      corr_pions += corr_pions_eff_out_of_cone;
    }
  }
  
  if ( useOutOfVertexTracks_ ) {
    corr_pions_out_of_vertex = pionCorrection( jet, pions.outOfVertexInCalo_, false, true );
    corr_pions += corr_pions_out_of_vertex;
  }
    
  if ( verbose_ ) {
    std::stringstream ss;
    ss << " Pion corrections:" << std::endl  
       << "  In/In      : " << "(" << pions.inVertexInCalo_.size() << ") " << corr_pions_in_cone.energy() << std::endl  
       << "  In/Out     : " << "(" << pions.inVertexOutOfCalo_.size() << ") " << corr_pions_out_of_cone.energy() << std::endl  
       << "  Out/In     : " << "(" << pions.outOfVertexInCalo_.size() << ") " << corr_pions_out_of_vertex.energy() << std::endl;
    if ( useEff_ ) {
      ss << " Pion efficiency corrections:" << std::endl  
	 << "  In/In      : " << "(" << pions.inVertexInCalo_.size() << ") " << corr_pions_eff_in_cone.energy() << std::endl  
	 << "  In/Out     : " << "(" << pions.inVertexOutOfCalo_.size() << ") " << corr_pions_eff_out_of_cone.energy();
    }
    edm::LogVerbatim("JPTCorrector") << ss.str();
  }

  return corr_pions;

}

// -----------------------------------------------------------------------------
//
JPTCorrector::P4 JPTCorrector::muonCorrection( const P4& jet,
					       const MatchedTracks& muons,
					       bool size ) const {
  
  P4 corr_muons;
  
  P4 corr_muons_in_cone;
  P4 corr_muons_out_of_cone;
  P4 corr_muons_out_of_vertex;
  
  if ( useInConeTracks_ ) { 
    corr_muons_in_cone = muonCorrection( jet, muons.inVertexInCalo_, true, true );
    corr_muons += corr_muons_in_cone;
  }  
  
  if ( useOutOfConeTracks_ ) {
    if ( size ) { //@@ bug?
      corr_muons_out_of_cone = muonCorrection( jet, muons.inVertexOutOfCalo_, true, false );
      corr_muons += corr_muons_out_of_cone;
    }
  }    
  
  if ( useOutOfVertexTracks_ ) {
    corr_muons_out_of_vertex = muonCorrection( jet, muons.outOfVertexInCalo_, false, true );
    corr_muons += corr_muons_out_of_vertex;
  }

  if ( verbose_ ) {
    std::stringstream ss;
    ss << " Muon corrections:" << std::endl  
       << "  In/In      : " << "(" << muons.inVertexInCalo_.size() << ") " << corr_muons_in_cone.energy() << std::endl  
       << "  In/Out     : " << "(" << muons.inVertexOutOfCalo_.size() << ") " << corr_muons_out_of_cone.energy() << std::endl  
       << "  Out/In     : " << "(" << muons.outOfVertexInCalo_.size() << ") " << corr_muons_out_of_vertex.energy();
    edm::LogVerbatim("JPTCorrector") << ss.str();
  }

  return corr_muons;

}

// -----------------------------------------------------------------------------
//
JPTCorrector::P4 JPTCorrector::elecCorrection( const P4& jet,
					       const MatchedTracks& elecs ) const {

  P4 null; //@@ null 4-momentum
  
  if ( verbose_ ) {
    std::stringstream ss;
    ss << " Electron corrections:" << std::endl  
       << "  In/In      : " << "(" << elecs.inVertexInCalo_.size() << ") " << null.energy() << std::endl  
       << "  In/Out     : " << "(" << elecs.inVertexOutOfCalo_.size() << ") " << null.energy() << std::endl  
       << "  Out/In     : " << "(" << elecs.outOfVertexInCalo_.size() << ") " << null.energy();
    edm::LogVerbatim("JPTCorrector") << ss.str();
  }

  return null; //@@ to be implemented

}

// -----------------------------------------------------------------------------
//
JPTCorrector::P4 JPTCorrector::jetDirFromTracks( const P4& corrected,
						 const MatchedTracks& pions,
						 const MatchedTracks& muons,
						 const MatchedTracks& elecs ) const {
  
  // Correction to be applied to jet 4-momentum
  P4 corr;
  
  bool tracks_present = false;
  
  // Correct using pions in-cone at vertex

  if ( !pions.inVertexInCalo_.empty() ) {
    TrackRefs::iterator itrk = pions.inVertexInCalo_.begin();
    TrackRefs::iterator jtrk = pions.inVertexInCalo_.end();
    for ( ; itrk != jtrk; ++itrk ) {
      if ( (*itrk)->pt() >= 50. ) { continue; }
      corr += PtEtaPhiM( (*itrk)->pt(), (*itrk)->eta(), (*itrk)->phi(), 0. );
      tracks_present = true;
    }
  }

  if ( !pions.inVertexOutOfCalo_.empty() ) {
    TrackRefs::iterator itrk = pions.inVertexOutOfCalo_.begin();
    TrackRefs::iterator jtrk = pions.inVertexOutOfCalo_.end();
    for ( ; itrk != jtrk; ++itrk ) {
      if ( (*itrk)->pt() >= 50. ) { continue; }
      corr += PtEtaPhiM( (*itrk)->pt(), (*itrk)->eta(), (*itrk)->phi(), 0. );
      tracks_present = true;
    }
  }

  // Correct using muons in-cone at vertex

  if ( !muons.inVertexInCalo_.empty() ) {
    TrackRefs::iterator itrk = muons.inVertexInCalo_.begin();
    TrackRefs::iterator jtrk = muons.inVertexInCalo_.end();
    for ( ; itrk != jtrk; ++itrk ) {
      if ( (*itrk)->pt() >= 50. ) { continue; }
      corr += PtEtaPhiM( (*itrk)->pt(), (*itrk)->eta(), (*itrk)->phi(), 0. );
      tracks_present = true;
    }
  }

  if ( !muons.inVertexOutOfCalo_.empty() ) {
    TrackRefs::iterator itrk = muons.inVertexOutOfCalo_.begin();
    TrackRefs::iterator jtrk = muons.inVertexOutOfCalo_.end();
    for ( ; itrk != jtrk; ++itrk ) {
      if ( (*itrk)->pt() >= 50. ) { continue; }
      corr += PtEtaPhiM( (*itrk)->pt(), (*itrk)->eta(), (*itrk)->phi(), 0. );
      tracks_present = true;
    }
  }

  // Correct using electrons in-cone at vertex

  if ( !elecs.inVertexInCalo_.empty() ) {
    TrackRefs::iterator itrk = elecs.inVertexInCalo_.begin();
    TrackRefs::iterator jtrk = elecs.inVertexInCalo_.end();
    for ( ; itrk != jtrk; ++itrk ) {
      if ( (*itrk)->pt() >= 50. ) { continue; }
      corr += PtEtaPhiM( (*itrk)->pt(), (*itrk)->eta(), (*itrk)->phi(), 0. );
      tracks_present = true;
    }
  }
  
  if ( !elecs.inVertexOutOfCalo_.empty() ) {
    TrackRefs::iterator itrk = elecs.inVertexOutOfCalo_.begin();
    TrackRefs::iterator jtrk = elecs.inVertexOutOfCalo_.end();
    for ( ; itrk != jtrk; ++itrk ) {
      if ( (*itrk)->pt() >= 50. ) { continue; }
      corr += PtEtaPhiM( (*itrk)->pt(), (*itrk)->eta(), (*itrk)->phi(), 0. );
      tracks_present = true;
    }
  }
  
  //@@ Scale to corrected jet energy
  if ( !tracks_present ) { corr = corrected; }
  else { corr *= ( corr.energy() > 0. ? corrected.energy() / corr.energy() : 1. ); }
  
  return corr;
  
}

// -----------------------------------------------------------------------------
//
JPTCorrector::P4 JPTCorrector::calculateCorr( const P4& jet,
					      const TrackRefs& tracks, 
					      bool in_cone_at_vertex,
					      bool in_cone_at_calo_face,
					      double mass, 
					      bool is_mip,
					      double mip ) const { 

  // Correction to be applied to jet 4-momentum
  P4 correction;

  // Reset efficiency container
  eff_.reset();
  
  // Iterate through tracks
  if ( !tracks.empty() ) {
    TrackRefs::iterator itrk = tracks.begin();
    TrackRefs::iterator jtrk = tracks.end();
    for ( ; itrk != jtrk; ++itrk ) {

      // Ignore high-pt tracks (only when in-cone and not a mip)
      if ( in_cone_at_calo_face && !is_mip && (*itrk)->pt() >= 50. ) { continue; }
      
      // Inner track 4-momentum
      P4 inner;
      if ( vectorialCorrection() ) {
	inner = PtEtaPhiM( (*itrk)->pt(), (*itrk)->eta(), (*itrk)->phi(), mass );
      } else { 
	double energy = sqrt( (*itrk)->px() * (*itrk)->px() + 
			      (*itrk)->py() * (*itrk)->py() + 
			      (*itrk)->pz() * (*itrk)->pz() + 
			      mass * mass );
	inner = ( jet.energy() > 0. ? energy / jet.energy() : 1. ) * jet;
      }      
      
      // Add track momentum (if in-cone at vertex)
      if ( in_cone_at_vertex ) { correction += inner; }
      
      // Find appropriate eta/pt bin for given track
      double eta = fabs( (*itrk)->eta() );
      double pt = fabs( (*itrk)->pt() );
      uint32_t ieta = response_.etaBin( eta );
      uint32_t ipt = response_.ptBin( pt );

      // Check bins (not for mips)
      if ( !is_mip && ( ieta == response_.nEtaBins() || ipt == response_.nPtBins() ) ) { continue; }
      
      // Outer track 4-momentum 
      P4 outer;
      if ( in_cone_at_calo_face ) { 
	if ( vectorialCorrection() ) {
	  // Build 4-momentum from outer track (SHOULD USE IMPACT POINT?!)
	  double outer_pt  = (*itrk)->pt();
	  double outer_eta = (*itrk)->eta();
	  double outer_phi = (*itrk)->phi();
	  if ( (*itrk)->extra().isNonnull() ) {
	    outer_pt  = (*itrk)->pt();
	    outer_eta = (*itrk)->outerPosition().eta(); //@@ outerMomentum().eta()
	    outer_phi = (*itrk)->outerPosition().phi(); //@@ outerMomentum().phi()
	  }
	  outer = PtEtaPhiM( outer_pt, outer_eta, outer_phi, mass );
	  // Check if mip or not
	  if ( is_mip ) { outer *= ( outer.energy() > 0. ? mip / outer.energy() : 1. ); } //@@ Scale to mip energy
	  else { outer *= ( outer.energy() > 0. ? inner.energy() / outer.energy() : 1. ); } //@@ Scale to inner track energy
	} else {
	  // Check if mip or not
	  if ( is_mip ) { outer = ( jet.energy() > 0. ? mip / jet.energy() : 1. ) * jet; } //@@ Jet 4-mom scaled by mip energy
	  else { outer = inner; } //@@ Set to inner track 4-momentum
	}      
	if ( !is_mip ) { outer *= response_.value(ieta,ipt); } //@@ Scale by pion response
	correction -= outer; //@@ Subtract 
      }
      
      // Record inner track energy for pion efficiency correction
      if ( !is_mip ) { eff_.addE( ieta, ipt, inner.energy() ); }
	
      // Debug
      if ( verbose_ && edm::isDebugEnabled() ) {
	std::stringstream temp; 
	temp << " Response[" << ieta << "," << ipt << "]";
	std::stringstream ss;
	ss << "[JPTCorrector::" << __func__ << "]" << std::endl
	   << " Track eta / pt    : " << eta << " / " << pt << std::endl
	   << temp.str() << std::setw(21-temp.str().size()) << " : " 
	   << response_.value(ieta,ipt) << std::endl
	   << " Track momentum added : " << inner.energy() << std::endl
	   << " Response subtracted  : " << outer.energy() << std::endl
	   << " Energy correction    : " << correction.energy();
	LogDebug("JPTCorrector") << ss.str();
      }
	
    } // loop through tracks
  } // ntracks != 0

  return correction;

}

// -----------------------------------------------------------------------------
//
JPTCorrector::P4 JPTCorrector::pionEfficiency( const P4& jet,
					       bool in_cone_at_calo_face ) const { 
  
  // Total correction to be applied
  P4 correction;
  
  // Iterate through eta/pt bins
  for ( uint32_t ieta = 0; ieta < response_.nEtaBins()-1; ++ieta ) {
    for ( uint32_t ipt = 0; ipt < response_.nPtBins()-1; ++ipt ) {

      // Check tracks are found in this eta/pt bin
      if ( !eff_.nTrks(ieta,ipt) ) { continue; }

      for ( uint16_t ii = 0; ii < 2; ++ii  ) {
	
	// Check which correction should be applied
	double corr = 0.;
	if ( ii == 0 )                              { corr = eff_.outOfConeCorr( ieta, ipt ); }
	else if ( ii == 1 && in_cone_at_calo_face ) { corr = eff_.inConeCorr( ieta, ipt ); }
	else                                        { continue; }

	// Calculate correction to be applied	
	P4 corr_p4;
	if ( vectorialCorrection() ) {
	  double corr_eta = response_.binCenterEta(ieta);
	  double corr_phi = jet.phi(); //@@ jet phi!
	  double corr_pt  = response_.binCenterPt(ipt);
	  corr_p4 = PtEtaPhiM( corr_pt, corr_eta, corr_phi, pionMass_ ); //@@ E^2 = p^2 + m_pion^2, |p| calc'ed from pt bin
	  corr_p4 *= ( corr_p4.energy() > 0. ? corr / corr_p4.energy() : 1. ); //@@ p4 scaled up by mean energy for bin
	} else { 
	  corr_p4 = ( jet.energy() > 0. ? corr / jet.energy() : 1. ) * jet;
	}      

	// Apply correction
	if ( ii == 0 )      { correction += corr_p4; } //@@ Add out-of-cone
	else if ( ii == 1 ) { correction -= corr_p4; } //@@ Subtract in-cone

      }	

    }
  }

  return correction;

}

// -----------------------------------------------------------------------------
//
bool JPTCorrector::matchMuons( TrackRefs::const_iterator itrk, 
			       const edm::Handle<RecoMuons>& muons ) const {
  
  if ( muons->empty() ) { return false; }

  RecoMuons::const_iterator imuon = muons->begin(); 
  RecoMuons::const_iterator jmuon = muons->end(); 
  for ( ; imuon != jmuon; ++imuon ) {
    
    if ( imuon->innerTrack().isNull() ||
	 !muon::isGoodMuon(*imuon,muon::TMLastStationTight) ||
	 imuon->innerTrack()->pt() < 3.0 ) { continue; }
    
    if ( itrk->id() != imuon->innerTrack().id() ) {
      edm::LogError("JPTCorrector")
	<< "[JPTCorrector::" << __func__ << "]"
	<< "Product id of the tracks associated to the jet " << itrk->id() 
	<<" is different from the product id of the inner track used for muons " << imuon->innerTrack().id()
	<< "!" << std::endl
	<< "Cannot compare tracks from different collection. Configuration Error!";
      return false;
    }
    
    if ( *itrk == imuon->innerTrack() ) { return true; }
  }
  
  return false;
  
}

// -----------------------------------------------------------------------------
//
bool JPTCorrector::matchElectrons( TrackRefs::const_iterator itrk, 
				   const edm::Handle<RecoElectrons>& elecs,
				   const edm::Handle<RecoElectronIds>& elec_ids ) const {
  
  if ( elecs->empty() ) { return false; }
  
  double deltaR = 999.;
  double deltaRMIN = 999.;
	
  uint32_t electron_index = 0;
  RecoElectrons::const_iterator ielec = elecs->begin(); 
  RecoElectrons::const_iterator jelec = elecs->end(); 
  for ( ; ielec != jelec; ++ielec ) {
    
    edm::Ref<RecoElectrons> electron_ref( elecs, electron_index );
    electron_index++;
    
    if ( (*elec_ids)[electron_ref] < 1.e-6 ) { continue; } //@@ Check for null value 
    
    // DR matching b/w electron and track
    double deltaphi = fabs( ielec->phi() - (*itrk)->momentum().phi() );
    if ( deltaphi > 6.283185308 ) deltaphi -= 6.283185308;
    if ( deltaphi > 3.141592654 ) deltaphi = 6.283185308 - deltaphi;
    deltaR = abs( sqrt( pow( (ielec->eta() - (*itrk)->momentum().eta()), 2 ) + 
			pow( deltaphi , 2 ) ) ); 
    if ( deltaR < deltaRMIN ) { deltaRMIN = deltaR; }
    
  }
  
  if ( deltaR < 0.02 ) return true;
  else return false;
  
}

// -----------------------------------------------------------------------------
//
void JPTCorrector::rebuildJta( const reco::Jet& fJet, 
			       const JetTracksAssociations& jtV0, 
			       TrackRefs& tracksthis,
			       TrackRefs& Excl ) const {
  
  //std::cout<<" NEW1 Merge/Split schema "<<jetSplitMerge_<<std::endl;

  tracksthis = reco::JetTracksAssociation::getValue(jtV0,fJet);

  if(jetSplitMerge_<0) return;

  typedef std::vector<reco::JetBaseRef>::iterator JetBaseRefIterator;
  std::vector<reco::JetBaseRef> theJets = reco::JetTracksAssociation::allJets(jtV0);

  TrackRefs tracks = tracksthis;
  tracksthis.clear();

  //std::cout<<" Size of initial vector "<<tracks.size()<<" "<<fJet.et()<<" "<<fJet.eta()<<" "<<fJet.phi()<<std::endl;

  int tr=0;

  for(TrackRefs::iterator it = tracks.begin(); it != tracks.end(); it++ )
    {

      double dR2this = deltaR2( fJet.eta(), fJet.phi(), (**it).eta(), (**it).phi() );
      //       double dfi = fabs(fJet.phi()-(**it).phi());
      //       if(dfi>4.*atan(1.))dfi = 8.*atan(1.)-dfi;
      //       double deta = fJet.eta() - (**it).eta();
      //       double dR2check = sqrt(dfi*dfi+deta*deta);
      
      double scalethis = dR2this;
      if(jetSplitMerge_ == 0) scalethis = 1./fJet.et();
      if(jetSplitMerge_ == 2) scalethis = dR2this/fJet.et();
      tr++;
      int flag = 1;
      for(JetBaseRefIterator ii = theJets.begin(); ii != theJets.end(); ii++)
	{
	  if(&(**ii) == &fJet ) {continue;}
	  double dR2 = deltaR2( (*ii)->eta(), (*ii)->phi(), (**it).eta(), (**it).phi() );
	  double scale = dR2;
	  if(jetSplitMerge_ == 0) scale = 1./fJet.et();
	  if(jetSplitMerge_ == 2) scale = dR2/fJet.et();
	  if(scale < scalethis) flag = 0;

	  if(flag == 0) {
	    //std::cout<<" Track belong to another jet also "<<dR2<<" "<<
	    //(*ii)->et()<<" "<<(*ii)->eta()<<" "<< (*ii)->phi()<<" Track "<<(**it).eta()<<" "<<(**it).phi()<<" "<<scalethis<<" "<<scale<<" "<<flag<<std::endl;
	    break;
	  }
	}

      //std::cout<<" Track "<<tr<<" "<<flag<<" "<<dR2this<<" "<<dR2check<<" Jet "<<fJet.eta()<<" "<< fJet.phi()<<" Track "<<(**it).eta()<<" "<<(**it).phi()<<std::endl;
      if(flag == 1) {tracksthis.push_back (*it);}else{Excl.push_back (*it);}
    }

  //std::cout<<" The new size of tracks "<<tracksthis.size()<<" Excludede "<<Excl.size()<<std::endl;
  return;
  
}

// -----------------------------------------------------------------------------
//
void JPTCorrector::excludeJta( const reco::Jet& fJet, 
			       const JetTracksAssociations& jtV0, 
			       TrackRefs& tracksthis,
			       const TrackRefs& Excl ) const {
  
  //std::cout<<" NEW2" << std::endl;

  tracksthis = reco::JetTracksAssociation::getValue(jtV0,fJet);
  if(Excl.size() == 0) return;
  if(jetSplitMerge_<0) return;

  TrackRefs tracks = tracksthis;
  tracksthis.clear();
  
  //std::cout<<" Size of initial vector "<<tracks.size()<<" "<<fJet.et()<<" "<<fJet.eta()<<" "<<fJet.phi()<<std::endl;

  for(TrackRefs::iterator it = tracks.begin(); it != tracks.end(); it++ )
    {

      //std::cout<<" Track at calo surface "
      //<<" Track "<<(**it).eta()<<" "<<(**it).phi()<<std::endl;
      TrackRefs::iterator itold = find(Excl.begin(),Excl.end(),(*it));
      if(itold == Excl.end()) {
	tracksthis.push_back (*it);
      } 
      //else { std::cout<<"Exclude "<<(**it).eta()<<" "<<(**it).phi()<<std::endl; }

    }

  //std::cout<<" Size of calo tracks "<<tracksthis.size()<<std::endl;

  return;

}

// -----------------------------------------------------------------------------
//
Map::Map( std::string input, bool verbose )
  : eta_(),
    pt_(),
    data_()
{ 

  // Some init
  clear();
  std::vector<Element> data;

  // Parse file
  std::string file = edm::FileInPath(input).fullPath();
  std::ifstream in( file.c_str() );
  string line;
  uint32_t ieta_old = 0; 
  while ( std::getline( in, line ) ) {
    if ( !line.size() || line[0]=='#' ) { continue; }
    std::istringstream ss(line);
    Element temp;
    ss >> temp.ieta_ >> temp.ipt_ >> temp.eta_ >> temp.pt_ >> temp.val_;
    data.push_back(temp);
    if ( !ieta_old || temp.ieta_ != ieta_old ) { 
      if ( eta_.size() < temp.ieta_+1 ) { eta_.resize(temp.ieta_+1,0.); }
      eta_[temp.ieta_] = temp.eta_;
      ieta_old = temp.ieta_;
    }
    if ( pt_.size() < temp.ipt_+1 ) { pt_.resize(temp.ipt_+1,0.); }
    pt_[temp.ipt_] = temp.pt_;
  }
  
  // Populate container
  data_.resize( eta_.size(), VDouble( pt_.size(), 0. ) );
  std::vector<Element>::const_iterator idata = data.begin();
  std::vector<Element>::const_iterator jdata = data.end();
  for ( ; idata != jdata; ++idata ) { data_[idata->ieta_][idata->ipt_] = idata->val_; }

  // Check
  if ( data_.empty() || data_[0].empty() ) {
    std::stringstream ss;
    ss << "[jpt::Map::" << __func__ << "]"
       << " Problem parsing map in location \"" 
       << file << "\"! ";
    edm::LogError("JPTCorrector") << ss.str();
  }

  // Check
  if ( eta_.size() != data_.size() || 
       pt_.size() != ( data_.empty() ? 0 : data_[0].size() ) ) {
    std::stringstream ss;
    ss << "[jpt::Map::" << __func__ << "]"
       << " Discrepancy b/w number of bins!";
    edm::LogError("JPTCorrector") << ss.str();
  }

  // Debug
  if ( verbose && edm::isDebugEnabled() ) { 
    std::stringstream ss;
    ss << "[jpt::Map::" << __func__ << "]"
       << " Parsed contents of map at location:" << std::endl
       << "\"" << file << "\"" << std::endl;
    print(ss); 
    LogTrace("JPTCorrector") << ss.str();
  } 

}

// -----------------------------------------------------------------------------
//
Map::Map() 
  : eta_(),
    pt_(),
    data_()
{ 
  clear();
}

// -----------------------------------------------------------------------------
//
Map::~Map() {
  clear();
}

// -----------------------------------------------------------------------------
//
void Map::clear() {
  eta_.clear();
  pt_.clear();
  data_.clear();
}
// -----------------------------------------------------------------------------
//
double Map::eta( uint32_t eta_bin ) const {
  if ( !eta_.empty() && eta_bin < eta_.size() ) { return eta_[eta_bin]; }
  else { 
    edm::LogWarning("JPTCorrector") 
      << "[jpt::Map::" << __func__ << "]"
      << " Trying to access element " << eta_bin
      << " of a vector with size " << eta_.size()
      << "!";
    return eta_[eta_.size()-1]; 
  }
}

// -----------------------------------------------------------------------------
//
double Map::pt( uint32_t pt_bin ) const {
  if ( !pt_.empty() && pt_bin < pt_.size() ) { return pt_[pt_bin]; }
  else { 
    edm::LogWarning("JPTCorrector") 
      << "[jpt::Map::" << __func__ << "]"
      << " Trying to access element " << pt_bin
      << " of a vector with size " << pt_.size()
      << "!";
    return pt_[pt_.size()-1]; 
  }
}
 
// -----------------------------------------------------------------------------
//
double Map::binCenterEta( uint32_t eta_bin ) const {
  if ( !eta_.empty() && eta_bin+1 < eta_.size() ) { 
    return eta_[eta_bin] + ( eta_[eta_bin+1] - eta_[eta_bin] ) / 2.; 
  } else { 
    edm::LogWarning("JPTCorrector") 
      << "[jpt::Map::" << __func__ << "]"
      << " Trying to access element " << eta_bin+1
      << " of a vector with size " << eta_.size()
      << "!";
    return eta_[eta_.size()-1]; 
  }
}
 
// -----------------------------------------------------------------------------
//
double Map::binCenterPt( uint32_t pt_bin ) const {
  if ( !pt_.empty() && pt_bin+1 < pt_.size() ) { 
    return pt_[pt_bin] + ( pt_[pt_bin+1] - pt_[pt_bin] ) / 2.; 
  } else { 
    edm::LogWarning("JPTCorrector") 
      << "[jpt::Map::" << __func__ << "]"
      << " Trying to access element " << pt_bin+1
      << " of a vector with size " << pt_.size()
      << "!";
    return pt_[pt_.size()-1]; 
  }
}

// -----------------------------------------------------------------------------
//
uint32_t Map::etaBin( double val ) const {
  val = fabs( val );
  for ( uint32_t ieta = 0; ieta < nEtaBins()-1; ++ieta ) { //@@ "-1" is bug?
    if ( val > eta(ieta) && ( ieta+1 == nEtaBins() || val < eta(ieta+1) ) ) { return ieta; }
  }
  return nEtaBins();
}

// -----------------------------------------------------------------------------
//
uint32_t Map::ptBin( double val ) const {
  val = fabs( val );
  for ( uint32_t ipt = 0; ipt < nPtBins()-1; ++ipt ) { //@@ "-1" is bug?
    if ( val > pt(ipt) && ( (ipt+1) == nPtBins() || val < pt(ipt+1) ) ) { return ipt; }
  }
  return nPtBins();
}

// -----------------------------------------------------------------------------
//
double Map::value( uint32_t eta_bin, uint32_t pt_bin ) const {
  if ( eta_bin < data_.size() && 
       pt_bin < ( data_.empty() ? 0 : data_[0].size() ) ) { return data_[eta_bin][pt_bin]; }
  else { 
    edm::LogWarning("JPTCorrector") 
      << "[jpt::Map::" << __func__ << "]"
      << " Trying to access element (" << eta_bin << "," << pt_bin << ")"
      << " of a vector with size (" << data_.size() << "," << ( data_.empty() ? 0 : data_[0].size() ) << ")"
      << "!";
    return 1.; 
  }
}

// -----------------------------------------------------------------------------
//
void Map::print( std::stringstream& ss ) const {
  ss << " Number of bins in eta : " << data_.size() << std::endl 
     << " Number of bins in pt  : " << ( data_.empty() ? 0 : data_[0].size() ) << std::endl;
  VVDouble::const_iterator ieta = data_.begin();
  VVDouble::const_iterator jeta = data_.end();
  for ( ; ieta != jeta; ++ieta ) {
    VDouble::const_iterator ipt = ieta->begin();
    VDouble::const_iterator jpt = ieta->end();
    for ( ; ipt != jpt; ++ipt ) {
      uint32_t eta_bin = static_cast<uint32_t>( ieta - data_.begin() );
      uint32_t pt_bin  = static_cast<uint32_t>( ipt - ieta->begin() );
      ss << " EtaBinNumber: " << eta_bin 
	 << " PtBinNumber: " << pt_bin 
	 << " EtaValue: " << eta_[ eta_bin ]
	 << " PtValue: " << pt_[ pt_bin ]
	 << " Value: " << data_[eta_bin][pt_bin]
	 << std::endl;
    }
  }
}

// -----------------------------------------------------------------------------
//
MatchedTracks::MatchedTracks() 
  : inVertexInCalo_(),
    outOfVertexInCalo_(),
    inVertexOutOfCalo_()
{ 
  clear();
}

// -----------------------------------------------------------------------------
//
MatchedTracks::~MatchedTracks() {
  clear();
}

// -----------------------------------------------------------------------------
//
void MatchedTracks::clear() {
  inVertexInCalo_.clear();
  outOfVertexInCalo_.clear();
  inVertexOutOfCalo_.clear();
}
 
// -----------------------------------------------------------------------------
//
JetTracks::JetTracks() 
  : vertex_(),
    caloFace_()
{ 
  clear();
}

// -----------------------------------------------------------------------------
//
JetTracks::~JetTracks() {
  clear();
}

// -----------------------------------------------------------------------------
//
void JetTracks::clear() {
  vertex_.clear();
  caloFace_.clear();
}

// -----------------------------------------------------------------------------
//
Efficiency::Efficiency( const jpt::Map& response,
			const jpt::Map& efficiency,
			const jpt::Map& leakage ) 
  : response_(response),
    efficiency_(efficiency),
    leakage_(leakage)
{
  reset();
}

// -----------------------------------------------------------------------------
//
Efficiency::Efficiency() 
  : response_(),
    efficiency_(),
    leakage_()
{;}

// -----------------------------------------------------------------------------
//
double Efficiency::inConeCorr( uint32_t eta_bin, uint32_t pt_bin ) const {
  if ( check(eta_bin,pt_bin,__func__) ) { 
    return ( outOfConeCorr( eta_bin, pt_bin ) * 
	     leakage_.value( eta_bin, pt_bin ) * 
	     response_.value( eta_bin, pt_bin ) ); 
  } else { return 0.; }
}

// -----------------------------------------------------------------------------
//
double Efficiency::outOfConeCorr( uint32_t eta_bin, uint32_t pt_bin ) const {
  if ( check(eta_bin,pt_bin,__func__) ) { 
    uint16_t ntrks = nTrks( eta_bin, pt_bin );
    double mean    = meanE( eta_bin, pt_bin );
    double eff     = ( 1. - efficiency_.value( eta_bin, pt_bin ) ) / efficiency_.value( eta_bin, pt_bin );
    if ( !ntrks ) { return 0.; }
    return ( ntrks * eff * mean ); 
  } else { return 0.; }
}

// -----------------------------------------------------------------------------
//
uint16_t Efficiency::nTrks( uint32_t eta_bin, uint32_t pt_bin ) const {
  if ( check(eta_bin,pt_bin,__func__) ) { 
    return data_[eta_bin][pt_bin].first; 
  } else { return 0; }
}

// -----------------------------------------------------------------------------
//
double Efficiency::sumE( uint32_t eta_bin, uint32_t pt_bin ) const {
  if ( check(eta_bin,pt_bin,__func__) ) { 
    return data_[eta_bin][pt_bin].second; 
  } else { return 0.; }
}

// -----------------------------------------------------------------------------
//
double Efficiency::meanE( uint32_t eta_bin, uint32_t pt_bin ) const {
  if ( check(eta_bin,pt_bin,__func__) ) { 
    Pair tmp = data_[eta_bin][pt_bin]; 
    if ( tmp.first ) { return tmp.second / tmp.first; }
    else { return 0.; }
  } else { return 0.; }
}

// -----------------------------------------------------------------------------
//
void Efficiency::addE( uint32_t eta_bin, uint32_t pt_bin, double energy ) {
  if ( check(eta_bin,pt_bin,__func__) ) { 
    data_[eta_bin][pt_bin].first++; 
    data_[eta_bin][pt_bin].second += energy;
  } 
}

// -----------------------------------------------------------------------------
//
bool Efficiency::check( uint32_t eta_bin, uint32_t pt_bin, std::string method ) const {
  if ( eta_bin < data_.size() && pt_bin < ( data_.empty() ? 0 : data_[0].size() ) ) { return true; }
  else { 
    edm::LogWarning("JPTCorrector") 
      << "[jpt::Efficiency::" << method << "]"
      << " Trying to access element (" << eta_bin << "," << pt_bin << ")"
      << " of a vector with size (" << data_.size() << "," << ( data_.empty() ? 0 : data_[0].size() ) << ")"
      << "!";
    return false; 
  }
}

// -----------------------------------------------------------------------------
//
void Efficiency::reset() { 
  data_.clear();
  data_.resize( response_.nEtaBins(), VPair( response_.nPtBins(), Pair(0,0.) ) );
}

// -----------------------------------------------------------------------------
//
void Efficiency::print() const { 
  if ( edm::isDebugEnabled() ) { 
    std::stringstream ss;
    ss << "[jpt::Efficiency::" << __func__ << "]"
       << " Contents of maps:" << std::endl;
    ss << "Response map: " << std::endl;
    response_.print(ss);
    ss << "Efficiency map: " << std::endl;
    efficiency_.print(ss);
    ss << "Leakage map: " << std::endl;
    leakage_.print(ss);
    LogTrace("JPTCorrector") << ss.str();
  } 
}
