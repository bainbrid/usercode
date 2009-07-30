#include "bainbrid/Test/test/JPTCorrector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 
#include "RecoJets/JetAssociationAlgorithms/interface/JetTracksAssociationDRCalo.h"
#include "RecoJets/JetAssociationAlgorithms/interface/JetTracksAssociationDRVertex.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

// -----------------------------------------------------------------------------
//
JPTCorrector::JPTCorrector( const edm::ParameterSet& pset ) 
  : jetTracksAtVertex_( pset.getParameter<edm::InputTag>("JetTracksAssociationAtVertex") ),
    jetTracksAtCalo_( pset.getParameter<edm::InputTag>("JetTracksAssociationAtCaloFace") ),
    tracks_( pset.getParameter<edm::InputTag>("Tracks") ),
    propagator_( pset.getParameter<std::string>("Propagator") ),
    coneSize_( pset.getParameter<double>("ConeSize") ),
    muons_( pset.getParameter<edm::InputTag>("Muons") ),
    electrons_( pset.getParameter<edm::InputTag>("Electrons") ),
    electronIds_( pset.getParameter<edm::InputTag>("ElectronIds") ),
    addOutOfConeTracks_( pset.getParameter<bool>("AddOutOfConeTracks") ),
    useTrackQuality_( pset.getParameter<bool>("UseTrackQuality") ),
    trackQuality_( reco::TrackBase::qualityByName( pset.getParameter<std::string>("TrackQuality") ) ),
    response_( pset.getParameter<std::string>("ResponseMap") ),
    efficiency_( pset.getParameter<std::string>("EfficiencyMap") ),
    leakage_( pset.getParameter<std::string>("LeakageMap") )
{  
  
//   std::string file_resp = pset.getParameter<std::string>("ResponseMap");
//   std::string file_eff  = pset.getParameter<std::string>("EfficiencyMap");
//   std::string file_leak = pset.getParameter<std::string>("LeakageMap");
  
//   std::string path   = "JetMETCorrections/Configuration/data/";
//   std::string suffix = ".txt";
  
//   file_resp = path + file_resp.substr( 0, file_resp.find_first_of(suffix) ) + suffix;
//   file_eff  = path + file_eff.substr( 0, file_eff.find_first_of(suffix) ) + suffix;
//   file_leak = path + file_leak.substr( 0, file_leak.find_first_of(suffix) ) + suffix;
  
//   std::stringstream ss;
  
//   if ( edm::isDebugEnabled() ) { 
//     ss << "[JPTCorrector::" << __func__ << "]"
//        << " Parsing response map   : " << file_resp << std::endl;
//   }
//   response_.clear();
//   response_ = Map( edm::FileInPath(file_resp).fullPath() );
  
//   if ( edm::isDebugEnabled() ) {
//     ss << "[JPTCorrector::" << __func__ << "]"
//        << " Parsing efficiency map : " << file_resp << std::endl;
//   }
//   efficiency_.clear();
//   efficiency_ = Map( edm::FileInPath(file_eff).fullPath() );
  
//   if ( edm::isDebugEnabled() ) {
//     ss << "[JPTCorrector::" << __func__ << "]"
//        << " Parsing leakage map    : " << file_resp << std::endl;
//   }
//   leakage_.clear();
//   leakage_ = Map( edm::FileInPath(file_leak).fullPath() );
  
//   if ( edm::isDebugEnabled() ) { LogTrace("JPTCorrector") << ss.str(); }
  
}

// -----------------------------------------------------------------------------
//
JPTCorrector::~JPTCorrector() {;}

// -----------------------------------------------------------------------------
//
double JPTCorrector::correction( const reco::Jet& fJet,
				 const edm::Event& event,
				 const edm::EventSetup& setup ) const 
{
  
  // Jet energy to correct
  double jet_energy = fJet.energy();
  
  // Check that jet falls within |eta| < 2.1
  if ( fabs( fJet.eta() ) > 2.1 ) { return jet_energy / fJet.energy(); }
  
  // Associate tracks to jet at both the Vertex and CaloFace
  JetTracks jet_tracks;
  bool ok = associateTracksToJets( fJet, event, setup, jet_tracks ); 
  if ( !ok ) { return ( jet_energy / fJet.energy() ); }
  
  // Track collections propagated to Vertex and CaloFace for "pions" and muons
  Tracks pions;
  Tracks muons;
  Tracks electrons;
  particleId( jet_tracks, event, setup, pions, muons, electrons );
  
  // -------------------- In-cone at Vertex and in-cone at CaloFace --------------------
  
  // Pions: subtract expected response and add track momentum
  TrackResponse in_cone;
  double corr_pions_in_cone = correction( pions.inVertexInCalo_, in_cone, true, true );
  jet_energy += corr_pions_in_cone;
  
  // Muons: subtract expected response and add track momentum
  TrackResponse not_used3;
  double corr_muons_in_cone = correction( muons.inVertexInCalo_, not_used3, true, true, 0.105, 2. );
  jet_energy += corr_muons_in_cone;

  // Pions: correct for tracking inefficiencies
  double corr_pion_eff_in_cone = correction( in_cone, true );
  jet_energy += corr_pion_eff_in_cone;

  // -------------------- In-cone at Vertex and out-of-cone at CaloFace --------------------

  // Pions: add track momentum 
  TrackResponse out_of_cone;
  double corr_pions_out_of_cone = correction( pions.inVertexOutOfCalo_, out_of_cone, false, true );
  jet_energy += corr_pions_out_of_cone;
  
  // Pions: add track momentum 
  TrackResponse not_used4;
  double corr_muons_out_of_cone = correction( muons.inVertexOutOfCalo_, out_of_cone, false, true , 0.105, 2. );
  jet_energy += corr_muons_out_of_cone;

  // Pions: correct for tracking inefficiencies
  double corr_pion_eff_out_of_cone = correction( out_of_cone, false );
  jet_energy += corr_pion_eff_out_of_cone;

  // -------------------- Out-of-cone at Vertex and in-cone at CaloFace --------------------

  // Pions: subtract expected response
  TrackResponse not_used1;
  double pion_out_of_vertex = correction( pions.outOfVertexInCalo_, not_used1, true, false );
  jet_energy += pion_out_of_vertex;

  // Muons: subtract expected response
  TrackResponse not_used2;
  double muons_out_of_vertex = correction( muons.outOfVertexInCalo_, not_used2, true, false, 0.105, 2. );
  jet_energy += muons_out_of_vertex;
   
  // -------------------- Return corrected energy -------------------- 
  
  return jet_energy / fJet.energy();

}

// -----------------------------------------------------------------------------
//
double JPTCorrector::correction( const reco::Jet& jet ) const {
  edm::LogError("Invalid correction use")
    << "JPTCorrector can be run on entire event only";
  return 1.;
}

// -----------------------------------------------------------------------------
//
double JPTCorrector::correction( const reco::Particle::LorentzVector& jet ) const {
  edm::LogError("Invalid correction use")
    << "JPTCorrector can be run on entire event only";
  return 1.;
}

// -----------------------------------------------------------------------------
//
bool JPTCorrector::associateTracksToJets( const reco::Jet& fJet,
					  const edm::Event& event, 
					  const edm::EventSetup& setup,
					  JetTracks& trks ) const {
  
  // Some init
  trks.clear();
  
  // Check whether to retrieve JTA object from Event or construct "on-the-fly"
  if ( !jetTracksAtVertex_.label().empty() && !jetTracksAtCalo_.label().empty() ) { 
    
    // Get Jet-track association at Vertex
    edm::Handle<reco::JetTracksAssociation::Container> jetTracksAtVertex;
    event.getByLabel( jetTracksAtVertex_, jetTracksAtVertex );
    
    // Check if handle is valid 
    if ( !jetTracksAtVertex.isValid() ) { return false; }
    
    // Retrieve jet-tracks association for given jet
    const reco::JetTracksAssociation::Container jtV = *( jetTracksAtVertex.product() );
    trks.atVertex_ = reco::JetTracksAssociation::getValue( jtV, fJet );
    
    // Check if any tracks are associated to jet at vertex
    if ( trks.atVertex_.empty() ) { return false; }
  
    // Get Jet-track association at Calo
    edm::Handle<reco::JetTracksAssociation::Container> jetTracksAtCalo;
    event.getByLabel( jetTracksAtCalo_, jetTracksAtCalo );
  
    // Check if handle is valid 
    if ( jetTracksAtCalo.isValid() ) { 
      // Retrieve jet-tracks association for given jet
      const reco::JetTracksAssociation::Container jtC = *( jetTracksAtCalo.product() );
      trks.atCaloFace_ = reco::JetTracksAssociation::getValue( jtC, fJet );
    }

  } else {

    // Construct objects that perform association 
    static JetTracksAssociationDRVertex vrtx(coneSize_);
    static JetTracksAssociationDRCalo   calo(coneSize_);

    // Container for propagated tracks
    static JetTracksAssociationDR::TrackRefs propagated;
    
    // Perform once per event
    static uint32_t last_event = 0;
    if ( event.id().event() != last_event ) {
      last_event = event.id().event();

      // Retrieve magnetic field and track propagator
      edm::ESHandle<MagneticField> field;
      setup.get<IdealMagneticFieldRecord>().get( field );
      edm::ESHandle<Propagator> propagator;
      setup.get<TrackingComponentsRecord>().get( propagator_, propagator );

      // Retrieve global tracks 
      edm::Handle<reco::TrackCollection> tracks;
      event.getByLabel( tracks_, tracks );
      if ( !tracks.isValid() || tracks.failedToGet() ) {
	edm::LogError("JEC")
	  << "[JPTCorrector::" << __func__ << "]"
	  << " Invalid handle to \"reco::TrackCollection\""
	  << " with InputTag (label:instance:process) \"" 
	  << tracks_.label() << ":"
	  << tracks_.instance() << ":"
	  << tracks_.process() << "\"";
	return false;
      }

      // Propagate tracks for to calo face 
      JetTracksAssociationDR::createTrackRefs( propagated, tracks, trackQuality_ );
      vrtx.propagateTracks( propagated ); //@@ needed?
      calo.propagateTracks( propagated, *field, *propagator );
      
    } 

    // Associate tracks to jets at both vertex and calo face
    vrtx.associateTracksToJet( trks.atVertex_, fJet, propagated );
    calo.associateTracksToJet( trks.atCaloFace_, fJet, propagated );
    
    // Check if any tracks are associated to jet at vertex
    if ( trks.atVertex_.empty() ) { return false; }
    
  }

  return true;
  
}

// -----------------------------------------------------------------------------
//
void JPTCorrector::particleId( const JetTracks& jet_tracks, 
			       const edm::Event& event, 
			       const edm::EventSetup& setup,
			       Tracks& pion_trks, 
			       Tracks& muon_trks,
			       Tracks& elec_trks ) const { 
  
  // Some init  
  pion_trks.clear(); 
  muon_trks.clear(); 
  elec_trks.clear(); 
  
  // Get muons
  edm::Handle<Muons> muon_collection;
  event.getByLabel( muons_, muon_collection );
  
  // Get electrons
  edm::Handle<Electrons> electron_collection;
  event.getByLabel( electrons_, electron_collection );
  
  // Get electron IDs
  edm::Handle<ElectronIDs> electron_ids;
  event.getByLabel( electronIds_, electron_ids );
  
  // Loop through tracks at "Vertex"
  {
    reco::TrackRefVector::const_iterator itrk = jet_tracks.atVertex_.begin();
    reco::TrackRefVector::const_iterator jtrk = jet_tracks.atVertex_.end();
    for ( ; itrk != jtrk; ++itrk ) {

      if ( useTrackQuality_ && !(*itrk)->quality(trackQuality_) ) { continue; }
      
      reco::TrackRefVector::iterator it = find( jet_tracks.atCaloFace_.begin(),
						jet_tracks.atCaloFace_.end(),
						*itrk );

      bool is_muon = matching( itrk, muon_collection );
      bool is_elec = false; //matching( itrk, electron_collection, electron_ids );
      
      if ( it != jet_tracks.atCaloFace_.end() ) { 
	if ( is_muon )      { muon_trks.inVertexInCalo_.push_back(*it); }
	else if ( is_elec ) { elec_trks.inVertexInCalo_.push_back(*it); } 
	else                { pion_trks.inVertexInCalo_.push_back(*it); } 
      } else { 
	if ( is_muon )      { muon_trks.inVertexOutOfCalo_.push_back(*itrk); }
	else                { pion_trks.inVertexOutOfCalo_.push_back(*itrk); }
      } 
    } 
    
  }
  
  // Loop through tracks at "CaloFace"
  {
    reco::TrackRefVector::iterator itrk = jet_tracks.atCaloFace_.begin(); 
    reco::TrackRefVector::iterator jtrk = jet_tracks.atCaloFace_.end(); 

    for ( ; itrk != jtrk; ++itrk ) {
      
      if ( useTrackQuality_ && !(*itrk)->quality(trackQuality_) ) { continue; }
      
      if( !pion_trks.inVertexInCalo_.empty() ) {
	
	reco::TrackRefVector::iterator it = find( pion_trks.inVertexInCalo_.begin(),
						  pion_trks.inVertexInCalo_.end(),
						  *itrk );
	
	reco::TrackRefVector::iterator im = find( muon_trks.inVertexInCalo_.begin(),
						  muon_trks.inVertexInCalo_.end(),
						  *itrk );
	
	if ( it == pion_trks.inVertexInCalo_.end() && 
	     im == muon_trks.inVertexInCalo_.end() ) {
	  
	  bool is_muon = matching( itrk, muon_collection );
	  bool is_elec = false; //@@ This is JW's implementation - needs checking!
	  
	  if ( is_muon )      { muon_trks.outOfVertexInCalo_.push_back(*itrk); } 
	  else if ( is_elec ) { elec_trks.outOfVertexInCalo_.push_back(*itrk); } 
	  else                { pion_trks.outOfVertexInCalo_.push_back(*itrk); }
	  
	}
      }
    } 
    
  }
  
  if ( edm::isDebugEnabled() ) {
    LogTrace("JPTCorrector")
      << " Size of theRecoTracksInConeInVertex " << pion_trks.inVertexInCalo_.size() << std::endl  
      << " Size of theRecoTracksOutOfConeInVertex " << pion_trks.inVertexOutOfCalo_.size() << std::endl
      << " Size of theRecoTracksInConeOutOfVertex " << pion_trks.outOfVertexInCalo_.size() << std::endl 
      << " Size of theRecoMuonTracksInConeInVertex " << muon_trks.inVertexInCalo_.size() << std::endl  
      << " Size of theRecoMuonTracksOutOfConeInVertex " << muon_trks.inVertexOutOfCalo_.size() << std::endl
      << " Size of theRecoMuonTracksInConeOutOfVertex " << muon_trks.outOfVertexInCalo_.size() << std::endl
      << " muon collection size = " << muon_collection->size() << std::endl;
  }

}

// -----------------------------------------------------------------------------
//
bool JPTCorrector::matching( reco::TrackRefVector::const_iterator itrk, 
			     const edm::Handle<Muons>& muons ) const {
  
  if ( muons->empty() ) { return false; }

  Muons::const_iterator imuon = muons->begin(); 
  Muons::const_iterator jmuon = muons->end(); 
  for ( ; imuon != jmuon; ++imuon ) {
    
    if ( !imuon->isGood(reco::Muon::TMLastStationTight) ) { continue; }
    if ( imuon->innerTrack()->pt() < 3.0 ) { continue; }
    
    if ( itrk->id() != imuon->innerTrack().id() ) {
      edm::LogError("JPTCorrector")
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
bool JPTCorrector::matching( reco::TrackRefVector::const_iterator itrk, 
			     const edm::Handle<Electrons>& electrons,
			     const edm::Handle<ElectronIDs>& electron_id ) const {

  if ( electrons->empty() ) { return false; }
  
  double deltaR = 999.;
  double deltaRMIN = 999.;
	
  uint32_t electron_index = 0;
  Electrons::const_iterator ielec = electrons->begin(); 
  Electrons::const_iterator jelec = electrons->end(); 
  for ( ; ielec != jelec; ++ielec ) {
    
    edm::Ref<Electrons> electron_ref( electrons, electron_index );

    if ( (*electron_id)[electron_ref] < 1.e-6 ) { continue; } //@@ Check for null value 

    // DR matching b/w electron and track
    double deltaphi = fabs( ielec->phi() - (*itrk)->momentum().phi() );
    if ( deltaphi > 6.283185308 ) deltaphi -= 6.283185308;
    if ( deltaphi > 3.141592654 ) deltaphi = 6.283185308 - deltaphi;
    deltaR = abs( sqrt( pow( (ielec->eta() - (*itrk)->momentum().eta()), 2 ) + 
			pow( deltaphi , 2 ) ) ); 
    if ( deltaR < deltaRMIN ) { deltaRMIN = deltaR; }

    electron_index++;
  }
  
  if ( deltaR < 0.02 ) return true;
  else return false;
  
}

// -----------------------------------------------------------------------------
//
double JPTCorrector::correction( const reco::TrackRefVector& tracks, 
				 TrackResponse& response,
				 bool subtract_response,
				 bool add_momentum,
				 double mass, 
				 double mip ) const { 
  
  // Correction to be applied
  double correction = 0.;
  
  // Clear and reset the response matrix
  response.clear();
  response.resize( response_.nEtaBins(), response_.nPtBins() );

  // Iterate through tracks
  if ( !tracks.empty() ) {
    reco::TrackRefVector::iterator itrk = tracks.begin();
    reco::TrackRefVector::iterator jtrk = tracks.end();
    for ( ; itrk != jtrk; ++itrk ) {

      // Track momentum
      double momentum = sqrt( (*itrk)->px() * (*itrk)->px() + 
			      (*itrk)->py() * (*itrk)->py() + 
			      (*itrk)->pz() * (*itrk)->pz() + 
			      mass * mass );
      
      // Add track momentum (if appropriate)
      if ( add_momentum ) { correction += momentum; }

      // Check if particle is mip or not
      if ( mip > 0. ) { 
	if ( subtract_response ) { correction -= mip; }
      } else { 
	// Find appropriate eta/pt bin for given track
	for ( uint32_t ieta = 0; ieta < response_.nEtaBins(); ++ieta ) {
	  for ( uint32_t ipt = 0; ipt < response_.nPtBins(); ++ipt ) {
	    double eta = fabs( (*itrk)->eta() );
	    if ( ( eta > response_.eta(ieta) && ieta+1 == response_.nEtaBins() ) ||
		 ( eta > response_.eta(ieta) && eta < response_.eta(ieta+1) ) ) {
	      double pt = fabs( (*itrk)->pt() );
	      if ( ( pt > response_.pt(ipt) && ipt+1 == response_.nPtBins() ) ||
		   ( pt > response_.pt(ipt) && pt < response_.pt(ipt+1) ) ) {
		
		// Subtract expected response (if appropriate)
		if ( subtract_response ) { correction -= ( momentum * response_.value( ieta, ipt ) ); } 
		
		// Record track momentum for efficiency correction
		response.addE( ieta, ipt, momentum );
		
	      }
	    }
	  }
	}
      }

    }
  }

  return correction;

}

// -----------------------------------------------------------------------------
//
double JPTCorrector::correction( TrackResponse& response,
				 bool subtract_response ) const { 
  
  // Correction to be applied
  double correction = 0.;
  
  // Iterate through eta/pt bins
  for ( uint32_t ieta = 0; ieta < response_.nEtaBins(); ++ieta ) {
    for ( uint32_t ipt = 0; ipt < response_.nPtBins(); ++ipt ) {
      double ntrks = response.nTrks( ieta, ipt );
      double mean  = response.meanE( ieta, ipt );
      double eff   = ( 1. - efficiency_.value(ieta,ipt) ) / efficiency_.value(ieta,ipt);
      correction  += ntrks * eff * mean;
      if ( subtract_response ) { 
	correction -= ntrks * eff * mean * leakage_.value(ieta,ipt) * response_.value(ieta,ipt);
      }
    }
  }

  return correction;

}


// -----------------------------------------------------------------------------
//
JPTCorrector::Map::Map( std::string file ) { 

  // Clear map
  clear();

  // Some debug
  if ( edm::isDebugEnabled() ) {
    std::stringstream ss;
    ss << "[JPTCorrector::" << __func__ << "] Parsing map from file \""
       << file 
       << "\"...";
    edm::LogVerbatim("JPTCorrector") << ss.str();
  }

  // Some init
  std::stringstream sss;
  std::ifstream in( file.c_str() );
  int ieta_prev = -1; 
  string line;

  // Parse file
  while ( std::getline( in, line ) ) {

    if ( !line.size() || line[0]=='#' ) { continue; }
    
    istringstream linestream(line);
    int ieta, ipt;
    double eta, pt, val;
    linestream >> ieta >> ipt >> eta >> pt >> val;
    
    if ( edm::isDebugEnabled() ) {
      sss << " EtaBinNumber: " << ieta
	  << " PtBinNumber: " << ipt
	  << " EtaValue: " << eta
	  << " PtValue: " << pt
	  << " Value: " << val
	  << std::endl;
    }      
    
    if ( ieta != ieta_prev ) {
      data_.resize( ieta+1, data_.size2() );
      eta_.push_back(eta);
      nEtaBins_ = ieta+1;
      ieta_prev = ieta;
    }
    
    if ( ieta_prev == 0 ) {
      data_.resize( data_.size1(), ipt+1 );
      pt_.push_back(pt); 
      nPtBins_ = ipt+1;
    }

    value_.push_back(val);
    data_(ieta,ipt) = val;

  }
  
  if ( edm::isDebugEnabled() ) {
    LogTrace("JPTCorrector") << sss.str();
  }
  
  if ( edm::isDebugEnabled() ) {
    std::stringstream ss;
    ss << "[JPTCorrector::" << __func__ << "]"
       << " nEtaBins: " << nEtaBins_ 
       << " nPtBins: " << nPtBins_ 
       << " nValues: " << value_.size();
    edm::LogVerbatim("JPTCorrector") << ss.str();
  }
  
  if ( edm::isDebugEnabled() ) {
    std::stringstream sss;
    std::vector<double>::const_iterator ieta = eta_.begin();
    std::vector<double>::const_iterator jeta = eta_.end();
    for ( ; ieta != jeta; ++ieta ) {
      std::vector<double>::const_iterator ipt = pt_.begin();
      std::vector<double>::const_iterator jpt = pt_.end();
      for ( ; ipt != jpt; ++ipt ) {
	uint16_t iter = ( ieta - eta_.begin() ) * eta_.size() + ( ipt - pt_.begin() );
	sss << " EtaBinNumber: " << uint16_t( ieta - eta_.begin() )
	    << " PtBinNumber: " << uint16_t( ipt - pt_.begin() )
	    << " EtaValue: " << *ieta
	    << " PtValue: " << *ipt
	    << " Value: " << value_[iter]
	    << std::endl;
      }
    }
    edm::LogVerbatim("JPTCorrector") << sss.str();
  }

  if ( edm::isDebugEnabled() ) {
    std::stringstream sss;
    matrix::const_iterator1 ieta = data_.begin1();
    matrix::const_iterator1 jeta = data_.end1();
    for ( ; ieta != jeta; ++ieta ) {
      matrix::const_iterator2 ipt = data_.begin2();
      matrix::const_iterator2 jpt = data_.end2();
      for ( ; ipt != jpt; ++ipt ) {
	sss << " EtaBinNumber: " << uint16_t( ieta - data_.begin1() )
	    << " PtBinNumber: " << uint16_t( ipt - data_.begin2() )
	    << " EtaValue: " << eta_[ uint16_t( ieta - data_.begin1() ) ]
	    << " PtValue: " << pt_[ uint16_t( ipt - data_.begin2() ) ]
	    << " Value: " << data_( uint16_t( ieta - data_.begin1() ),
				    uint16_t( ipt - data_.begin2() ) )
	    << " test: " << *ieta << " " << *ipt
	    << std::endl;
      }
    }
    edm::LogVerbatim("JPTCorrector") << sss.str();
  }

}

// -----------------------------------------------------------------------------
//
JPTCorrector::Map::Map() 
  : eta_(),
    pt_(),
    data_(),
    nEtaBins_(0),
    nPtBins_(0),
    value_()
{ 
  clear();
}

// -----------------------------------------------------------------------------
//
JPTCorrector::Map::~Map() {
  clear();
}

// -----------------------------------------------------------------------------
//
void JPTCorrector::Map::clear() {
  eta_.clear();
  pt_.clear();
  data_.clear();
  nEtaBins_ = 0;
  nPtBins_ = 0;
  value_.clear();
}
 
// -----------------------------------------------------------------------------
//
double JPTCorrector::Map::eta( uint32_t eta_bin ) const {
  if ( eta_bin < eta_.size() ) { return eta_[eta_bin]; }
  else { 
    edm::LogWarning("JPTCorrector") 
      << "[JPTCorrector::Map::eta]"
      << " Trying to access element " << eta_bin
      << " of a vector with size " << eta_.size()
      << "!";
    return -1.; 
  }
}

// -----------------------------------------------------------------------------
//
double JPTCorrector::Map::pt( uint32_t pt_bin ) const {
  if ( pt_bin < pt_.size() ) { return pt_[pt_bin]; }
  else { 
    edm::LogWarning("JPTCorrector") 
      << "[JPTCorrector::Map::eta]"
      << " Trying to access element " << pt_bin
      << " of a vector with size " << pt_.size()
      << "!";
    return -1.; 
  }
}

// -----------------------------------------------------------------------------
//
double JPTCorrector::Map::value( uint32_t eta_bin, uint32_t pt_bin ) const {
  if ( eta_bin < data_.size1() && pt_bin < data_.size2() ) { return data_(eta_bin,pt_bin); }
  else { 
    edm::LogWarning("JPTCorrector") 
      << "[JPTCorrector::Map::eta]"
      << " Trying to access element (" << eta_bin << "," << pt_bin << ")"
      << " of a vector with size (" << data_.size1() << "," << data_.size2() << ")"
      << "!";
    return 1.; 
  }
}
 
// -----------------------------------------------------------------------------
//
JPTCorrector::Tracks::Tracks() 
  : inVertexInCalo_(),
    outOfVertexInCalo_(),
    inVertexOutOfCalo_()
{ 
  clear();
}

// -----------------------------------------------------------------------------
//
JPTCorrector::Tracks::~Tracks() {
  clear();
}

// -----------------------------------------------------------------------------
//
void JPTCorrector::Tracks::clear() {
  inVertexInCalo_.clear();
  outOfVertexInCalo_.clear();
  inVertexOutOfCalo_.clear();
}
 
// -----------------------------------------------------------------------------
//
JPTCorrector::JetTracks::JetTracks() 
  : atVertex_(),
    atCaloFace_()
{ 
  clear();
}

// -----------------------------------------------------------------------------
//
JPTCorrector::JetTracks::~JetTracks() {
  clear();
}

// -----------------------------------------------------------------------------
//
void JPTCorrector::JetTracks::clear() {
  atVertex_.clear();
  atCaloFace_.clear();
}

// -----------------------------------------------------------------------------
//
uint16_t JPTCorrector::TrackResponse::nTrks( uint32_t eta_bin, uint32_t pt_bin ) const {
  if ( check(eta_bin,pt_bin) ) { 
    return ( operator() (eta_bin,pt_bin) ).first; 
  } else { return 0; }
}

// -----------------------------------------------------------------------------
//
double JPTCorrector::TrackResponse::sumE( uint32_t eta_bin, uint32_t pt_bin ) const {
  if ( check(eta_bin,pt_bin) ) { 
    return ( operator() (eta_bin,pt_bin) ).second; 
  } else { return 0.; }
}

// -----------------------------------------------------------------------------
//
double JPTCorrector::TrackResponse::meanE( uint32_t eta_bin, uint32_t pt_bin ) const {
  if ( check(eta_bin,pt_bin) ) { 
    Pair tmp = operator() (eta_bin,pt_bin); 
    if ( tmp.first ) { return tmp.second / tmp.first; }
    else { return 0.; }
  } else { return 0.; }
}

// -----------------------------------------------------------------------------
//
void JPTCorrector::TrackResponse::addE( uint32_t eta_bin, uint32_t pt_bin, double energy ) {
  if ( check(eta_bin,pt_bin) ) { 
    Pair tmp = operator() (eta_bin,pt_bin); 
    tmp.first  += 1;
    tmp.second += energy;
  } 
}

// -----------------------------------------------------------------------------
//
bool JPTCorrector::TrackResponse::check( uint32_t eta_bin, uint32_t pt_bin ) const {
  if ( eta_bin < size1() && pt_bin < size2() ) { return true; }
  else { 
    edm::LogWarning("JPTCorrector") 
      << "[JPTCorrector::TrackResponse::check]"
      << " Trying to access element (" << eta_bin << "," << pt_bin << ")"
      << " of a vector with size (" << size1() << "," << size2() << ")"
      << "!";
    return false; 
  }
}










//   // Tracks that are out-of-cone at Vertex and in-cone at CaloFace
//   // Subtract expected response
//   if ( !propagated_pions.outOfVertexInCalo_.empty() ) {
//     reco::TrackRefVector::iterator itrk = propagated_pions.outOfVertexInCalo_.begin();
//     reco::TrackRefVector::iterator jtrk = propagated_pions.outOfVertexInCalo_.end();
//     for ( ; itrk != jtrk; ++itrk ) {
//       double momentum = sqrt( (*itrk)->px() * (*itrk)->px() + 
// 			      (*itrk)->py() * (*itrk)->py() + 
// 			      (*itrk)->pz() * (*itrk)->pz() + 
// 			      0.14 * 0.14 );
//       std::vector<double>::const_iterator ieta = response_.eta_.begin();
//       std::vector<double>::const_iterator jeta = response_.eta_.end();
//       for ( ; ieta != jeta; ++ieta ) {
// 	std::vector<double>::const_iterator ipt = response_.pt_.begin();
// 	std::vector<double>::const_iterator jpt = response_.pt_.end();
// 	for ( ; ipt != jpt; ++ipt ) {
// 	  double eta = fabs( (*itrk)->eta() );
// 	  if( ( eta > *ieta && ieta+1 == jeta ) ||
// 	      ( eta > *ieta && eta < *(ieta+1) ) ) {
// 	    double pt = fabs( (*itrk)->pt() );
// 	    if( ( pt > *ipt && ipt+1 == jpt ) ||
// 		( pt > *ipt && pt < *(ipt+1) ) ) {
// 	      uint16_t iter = ( jeta - ieta ) * response_.eta_.size() + ( jpt - ipt );
// 	      jet_energy -= ( response_.value_[iter] * momentum ); //@@ subtract expected response
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
  
//   std::vector<double> emean_incone( response_.nEtaBins_ * response_.nPtBins_, 0. );
//   std::vector<double> ntracks_incone( response_.nEtaBins_ * response_.nPtBins_, 0. );

//   // Tracks that are in-cone at Vertex and in-cone at CaloFace
//   // Subtract expected response and add track momentum 
//   if ( !propagated_pions.inVertexInCalo_.empty() ) {
//     reco::TrackRefVector::iterator itrk = propagated_pions.inVertexInCalo_.begin();
//     reco::TrackRefVector::iterator jtrk = propagated_pions.inVertexInCalo_.end();
//     for ( ; itrk != jtrk; ++itrk ) {
//       if ( (*itrk)->pt() >= 50. ) { continue; }
//       double momentum = sqrt( (*itrk)->px() * (*itrk)->px() + 
// 			      (*itrk)->py() * (*itrk)->py() + 
// 			      (*itrk)->pz() * (*itrk)->pz() + 
// 			      0.14 * 0.14 );
//       jet_energy += momentum; //@@ add track momentum
//       std::vector<double>::const_iterator ieta = response_.eta_.begin();
//       std::vector<double>::const_iterator jeta = response_.eta_.end();
//       for ( ; ieta != jeta; ++ieta ) {
// 	std::vector<double>::const_iterator ipt = response_.pt_.begin();
// 	std::vector<double>::const_iterator jpt = response_.pt_.end();
// 	for ( ; ipt != jpt; ++ipt ) {
// 	  double eta = fabs( (*itrk)->eta() );
// 	  if( ( eta > *ieta && ieta+1 == jeta ) ||
// 	      ( eta > *ieta && eta < *(ieta+1) ) ) {
// 	    double pt = fabs( (*itrk)->pt() );
// 	    if( ( pt > *ipt && ipt+1 == jpt ) ||
// 		( pt > *ipt && pt < *(ipt+1) ) ) {
// 	      uint16_t iter = ( jeta - ieta ) * response_.eta_.size() + ( jpt - ipt );
// 	      jet_energy -= response_.value_[iter] * momentum; //@@ subtract expected response
// 	      ntracks_incone[iter]++; 
// 	      emean_incone[iter] += momentum;
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
  
//   // Tracks that are in-cone at Vertex and in-cone at CaloFace
//   // Correct for tracking inefficiencies and energy leakage 
//   if ( !propagated_pions.inVertexInCalo_.empty() ) {
//     double correction = 0.;
//     std::vector<double>::const_iterator ieta = response_.eta_.begin();
//     std::vector<double>::const_iterator jeta = response_.eta_.end();
//     for ( ; ieta != jeta; ++ieta ) {
//       std::vector<double>::const_iterator ipt = response_.pt_.begin();
//       std::vector<double>::const_iterator jpt = response_.pt_.end();
//       for ( ; ipt != jpt; ++ipt ) {
// 	uint16_t iter = ( jeta - ieta ) * response_.eta_.size() + ( jpt - ipt );
// 	if ( iter < ntracks_incone.size() && ntracks_incone[iter] ) { 
// 	  emean_incone[iter] /= ntracks_incone[iter];
// 	  correction += ( ( ntracks_incone[iter] ) * 
// 			  ( ( 1. - efficiency_.value_[iter] ) / efficiency_.value_[iter] ) * emean_incone[iter] * 
// 			  ( ( 1. - leakage_.value_[iter] * response_.value_[iter] ) )
// 			  );
// 	}
//       }
//     }
//     jet_energy += correction;
//   }
  
//   std::vector<double> emean_outcone( response_.nEtaBins_ * response_.nPtBins_, 0. );
//   std::vector<double> ntracks_outcone( response_.nEtaBins_ * response_.nPtBins_, 0. );
  
//   if ( addOutOfConeTracks_ ) {

//     // Tracks that are in-cone at Vertex and out-of-cone at CaloFace
//     // Add track momentum 
//     if ( !propagated_pions.inVertexOutOfCalo_.empty() ) {
//       reco::TrackRefVector::iterator itrk = propagated_pions.inVertexOutOfCalo_.begin();
//       reco::TrackRefVector::iterator jtrk = propagated_pions.inVertexOutOfCalo_.end();
//       for ( ; itrk != jtrk; ++itrk ) {
// 	if ( (*itrk)->pt() >= 50. ) { continue; }
// 	double momentum = sqrt( (*itrk)->px() * (*itrk)->px() + 
// 				(*itrk)->py() * (*itrk)->py() + 
// 				(*itrk)->pz() * (*itrk)->pz() + 
// 				0.14 * 0.14 );
// 	jet_energy += momentum; //@@ add track momentum
// 	std::vector<double>::const_iterator ieta = response_.eta_.begin();
// 	std::vector<double>::const_iterator jeta = response_.eta_.end();
// 	for ( ; ieta != jeta; ++ieta ) {
// 	  std::vector<double>::const_iterator ipt = response_.pt_.begin();
// 	  std::vector<double>::const_iterator jpt = response_.pt_.end();
// 	  for ( ; ipt != jpt; ++ipt ) {
// 	    double eta = fabs( (*itrk)->eta() );
// 	    if( ( eta > *ieta && ieta+1 == jeta ) ||
// 		( eta > *ieta && eta < *(ieta+1) ) ) {
// 	      double pt = fabs( (*itrk)->pt() );
// 	      if( ( pt > *ipt && ipt+1 == jpt ) ||
// 		  ( pt > *ipt && pt < *(ipt+1) ) ) {
// 		uint16_t iter = ( jeta - ieta ) * response_.eta_.size() + ( jpt - ipt );
// 		if ( iter < ntracks_outcone.size() ) { 
// 		  ntracks_outcone[iter]++;
// 		  emean_outcone[iter] += momentum;
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//       }
//     }

//     // Tracks that are in-cone at Vertex and out-of-cone at CaloFace
//     // Correct for tracking inefficiencies
//     if ( !propagated_pions.inVertexOutOfCalo_.empty() ) {
//       double correction = 0.;
//       std::vector<double>::const_iterator ieta = response_.eta_.begin();
//       std::vector<double>::const_iterator jeta = response_.eta_.end();
//       for ( ; ieta != jeta; ++ieta ) {
// 	std::vector<double>::const_iterator ipt = response_.pt_.begin();
// 	std::vector<double>::const_iterator jpt = response_.pt_.end();
// 	for ( ; ipt != jpt; ++ipt ) {
// 	  uint16_t iter = ( jeta - ieta ) * response_.eta_.size() + ( jpt - ipt );
// 	  if ( iter < ntracks_outcone.size() && ntracks_outcone[iter] ) { 
// 	    emean_outcone[iter] /= ntracks_outcone[iter];
// 	    correction += ( ( ntracks_outcone[iter] ) * 
// 			    ( ( 1. - efficiency_.value_[iter] ) / efficiency_.value_[iter] ) * emean_outcone[iter]
// 			    );
// 	  }
// 	}
//       }
//       jet_energy += correction;
//     }
  
//   }


