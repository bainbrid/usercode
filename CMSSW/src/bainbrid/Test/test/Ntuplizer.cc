#include "Ntuplizer.h"
#include "CommonTools/Utils/interface/EtComparator.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cmath>

// -----------------------------------------------------------------------------
//
const uint32_t Ntuplizer::SIZE = 50;

// -----------------------------------------------------------------------------
//
Ntuplizer::Ntuplizer( const edm::ParameterSet& pset ) 
  : tree_(0),
    tags_(),
    names_(),
    n_(),
    e_(), 
    et_(),
    p_(),
    pt_(),
    px_(),
    py_(),
    pz_(),
    eta_(),
    phi_(),
    mass_(),
    charge_(),
    pdgId_(),
    vx_(),
    vy_(),
    vz_(),
    k_( pset.getUntrackedParameter<bool>("Kinematics",true) ),
    t_( pset.getUntrackedParameter<bool>("Transverse",true) ),
    d_( pset.getUntrackedParameter<bool>("Direction",true) ),
    i_( pset.getUntrackedParameter<bool>("ParticleId",true) ),
    v_( pset.getUntrackedParameter<bool>("Vertex",true) )
{
  
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("Data","data");
  tree_->SetAutoSave(1);
  
  std::vector<edm::InputTag> tags = pset.getParameter< std::vector<edm::InputTag> >("src");
  std::vector<edm::InputTag>::const_iterator ii = tags.begin();
  std::vector<edm::InputTag>::const_iterator jj = tags.end();
  for ( ; ii != jj; ++ii ) {
    if ( ii->label().empty() || ii->instance().empty() ) {
      edm::LogError("TEST") << "Null label and/or instance for: " << *ii; 
      continue;
    }
    std::vector<edm::InputTag>::iterator iii = tags_.begin();
    std::vector<edm::InputTag>::iterator jjj = tags_.end();
    for ( ; iii != jjj; ++iii ) { if ( ii->label() == iii->label() &&
				       ii->instance() == iii->instance() ) { break; } }
    if ( iii != tags_.end() ) {
      edm::LogError("TEST") << "InputTag " << *ii << " already exists in map";
      continue;
    }
    tags_.push_back( *ii );
  }
  
  init( tags_.size() );

}

// -----------------------------------------------------------------------------
//
void Ntuplizer::analyze( const edm::Event& event, 
			 const edm::EventSetup& setup ) {

  reset( tags_.size() );
  
  uint32_t index = 0;
  std::vector<edm::InputTag>::iterator ii = tags_.begin();
  std::vector<edm::InputTag>::iterator jj = tags_.end();
  for ( ; ii != jj; ++ii ) {
    
    edm::Handle< edm::View<reco::LeafCandidate> > handle;
    event.getByLabel(*ii,handle);
    if ( handle.isValid() && !handle.failedToGet() ) {
      edm::View<reco::LeafCandidate>::const_iterator iii = handle->begin();
      edm::View<reco::LeafCandidate>::const_iterator jjj = handle->end();
      for ( ; iii != jjj; ++iii ) { add( index, *iii ); } 
    } else {
      edm::LogError("TEST") << "Error retrieving LeafCandidate collection for: " << *ii; 
    }
    index++;
    
  }

  if ( tree_ ) { tree_->Fill(); }

}

// -----------------------------------------------------------------------------
//
void Ntuplizer::init( uint32_t size ) {

  // Clear containers
  n_.clear();
  e_.clear();
  et_.clear();
  p_.clear();
  pt_.clear();
  px_.clear();
  py_.clear();
  pz_.clear();
  eta_.clear();
  phi_.clear();
  mass_.clear();
  charge_.clear();
  pdgId_.clear();
  vx_.clear();
  vy_.clear();
  vz_.clear();

  // Resize containers
  n_.resize( size, 0 );
  e_.resize( size, VDouble(SIZE,0.) );
  et_.resize( size, VDouble(SIZE,0.) );
  p_.resize( size, VDouble(SIZE,0.) );
  pt_.resize( size, VDouble(SIZE,0.) );
  px_.resize( size, VDouble(SIZE,0.) );
  py_.resize( size, VDouble(SIZE,0.) );
  pz_.resize( size, VDouble(SIZE,0.) );
  eta_.resize( size, VDouble(SIZE,0.) );
  phi_.resize( size, VDouble(SIZE,0.) );
  mass_.resize( size, VDouble(SIZE,0.) );
  charge_.resize( size, VInt(SIZE,0) );
  pdgId_.resize( size, VInt(SIZE,0) );
  vx_.resize( size, VDouble(SIZE,0.) );
  vy_.resize( size, VDouble(SIZE,0.) );
  vz_.resize( size, VDouble(SIZE,0.) );

  // Define branches
  if ( tree_ ) {
    uint32_t index = 0;
    std::vector<edm::InputTag>::iterator ii = tags_.begin();
    std::vector<edm::InputTag>::iterator jj = tags_.end();
    for ( ; ii != jj; ++ii ) { 
      tree_->Branch( name(*ii,"Entries"), &n_[index], name(*ii,"Entries","i") );
      if ( k_ ) { tree_->Branch( name(*ii,"Energy"), &e_[index].front(), name(*ii,"Energy","D","Entries") ); }
      if ( k_ ) { tree_->Branch( name(*ii,"Momentum"), &p_[index].front(), name(*ii,"Momentum","D","Entries") ); }
      if ( t_ ) { tree_->Branch( name(*ii,"Et"), &et_[index].front(), name(*ii,"Et","D","Entries") ); }
      if ( t_ ) { tree_->Branch( name(*ii,"Pt"), &pt_[index].front(), name(*ii,"Pt","D","Entries") ); }
      if ( k_ ) { tree_->Branch( name(*ii,"Px"), &px_[index].front(), name(*ii,"Px","D","Entries") ); }
      if ( k_ ) { tree_->Branch( name(*ii,"Py"), &py_[index].front(), name(*ii,"Py","D","Entries") ); }
      if ( k_ ) { tree_->Branch( name(*ii,"Pz"), &pz_[index].front(), name(*ii,"Pz","D","Entries") ); }
      if ( d_ ) { tree_->Branch( name(*ii,"Eta"), &eta_[index].front(), name(*ii,"Eta","D","Entries") ); }
      if ( d_ ) { tree_->Branch( name(*ii,"Phi"), &phi_[index].front(), name(*ii,"Phi","D","Entries") ); }
      if ( i_ ) { tree_->Branch( name(*ii,"Mass"), &mass_[index].front(), name(*ii,"Mass","D","Entries") ); }
      if ( i_ ) { tree_->Branch( name(*ii,"Charge"), &charge_[index].front(), name(*ii,"Charge","i","Entries") ); }
      if ( i_ ) { tree_->Branch( name(*ii,"PdgId"), &pdgId_[index].front(), name(*ii,"PdgId","i","Entries") ); }
      if ( v_ ) { tree_->Branch( name(*ii,"Vx"), &vx_[index].front(), name(*ii,"Vx","D","Entries") ); }
      if ( v_ ) { tree_->Branch( name(*ii,"Vy"), &vy_[index].front(), name(*ii,"Vy","D","Entries") ); }
      if ( v_ ) { tree_->Branch( name(*ii,"Vz"), &vz_[index].front(), name(*ii,"Vz","D","Entries") ); }
      index++;
    }
  } else {
    edm::LogError("TEST") << "Null pointer to TTree";
  }

}

// -----------------------------------------------------------------------------
//
void Ntuplizer::reset( uint32_t size ) {
  std::memset( &n_.front(), 0, size*sizeof(Int) );
  for ( uint32_t index = 0; index < size; ++index ) {  
    std::memset( &et_[index].front(), 0, SIZE*sizeof(Double) );
    std::memset( &p_[index].front(), 0, SIZE*sizeof(Double) );
    std::memset( &pt_[index].front(), 0, SIZE*sizeof(Double) );
    std::memset( &px_[index].front(), 0, SIZE*sizeof(Double) );
    std::memset( &py_[index].front(), 0, SIZE*sizeof(Double) );
    std::memset( &pz_[index].front(), 0, SIZE*sizeof(Double) );
    std::memset( &eta_[index].front(), 0, SIZE*sizeof(Double) );
    std::memset( &phi_[index].front(), 0, SIZE*sizeof(Double) );
    std::memset( &mass_[index].front(), 0, SIZE*sizeof(Double) );
    std::memset( &charge_[index].front(), 0, SIZE*sizeof(Int) );
    std::memset( &pdgId_[index].front(), 0, SIZE*sizeof(Int) );
    std::memset( &vx_[index].front(), 0, SIZE*sizeof(Double) );
    std::memset( &vy_[index].front(), 0, SIZE*sizeof(Double) );
    std::memset( &vz_[index].front(), 0, SIZE*sizeof(Double) );
  }
}

// -----------------------------------------------------------------------------
//
void Ntuplizer::add( uint32_t index, const reco::LeafCandidate& cand ) {
  if ( index < n_.size() && n_[index] < SIZE ) {
    e_[index][n_[index]] = cand.energy();
    et_[index][n_[index]] = cand.et();
    p_[index][n_[index]] = cand.p();
    pt_[index][n_[index]] = cand.pt();
    px_[index][n_[index]] = cand.px();
    py_[index][n_[index]] = cand.py();
    pz_[index][n_[index]] = cand.pz();
    eta_[index][n_[index]] = cand.eta();
    phi_[index][n_[index]] = cand.phi();
    mass_[index][n_[index]] = cand.mass();
    charge_[index][n_[index]] = cand.charge();
    pdgId_[index][n_[index]] = cand.pdgId();
    vx_[index][n_[index]] = cand.vx();
    vy_[index][n_[index]] = cand.vy();
    vz_[index][n_[index]] = cand.vz();
    std::stringstream ss;
    if ( edm::isDebugEnabled() ) {
      std::string tag = "UNKNOWN";
      if ( index < tags_.size() ) { 
	tag = (tags_.begin()+index)->label(); 
	tag += ":";
	tag = (tags_.begin()+index)->instance(); 
      }
      ss << " Object #" << n_[index]
	 << " from \"" << tag << "\": "
	 << " E/p=" 
	 << e_[index][n_[index]] << "/"
	 << p_[index][n_[index]] 
	 << ", Et/Pt="
	 << et_[index][n_[index]] << "/"
	 << pt_[index][n_[index]] 
	 << ", Px/Py/Pz="
	 << px_[index][n_[index]] << "/"
	 << py_[index][n_[index]] << "/"
	 << pz_[index][n_[index]] 
	 << ", Eta/Phi="
	 << eta_[index][n_[index]] << "/"
	 << phi_[index][n_[index]];
      LogTrace("TEST") << ss.str();
    }
    n_[index]++;
  }
}

// -----------------------------------------------------------------------------
//
void Ntuplizer::beginJob( const edm::EventSetup& ) {;}

// -----------------------------------------------------------------------------
//
const char* Ntuplizer::name( const edm::InputTag& tag,
			     const std::string& name,
			     const std::string& type,
			     const std::string& index ) {
  if ( type.empty() ) {
    names_.push_back( std::string( tag.label() + "_" + 
				   tag.instance() + "_" + 
				   name ) );
  } else {
    if ( index.empty() ) {
      names_.push_back( std::string( tag.label() + "_" + 
				     tag.instance() + "_" + 
				     name + "/" + 
				     type ) );
    } else {
      names_.push_back( std::string( tag.label() + "_" + 
				     tag.instance() + "_" + 
				     name + "[" + 
				     tag.label() + "_" + 
				     tag.instance() + "_" + 
				     index + "]" + "/" + 
				     type ) );
    }
  }
  return names_.back().c_str();
}

// -----------------------------------------------------------------------------
// 
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Ntuplizer);
