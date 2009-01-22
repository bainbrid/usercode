#include "JetMETCorrections/JetPlusTrack/interface/ObjectMatcher.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <sstream>

// -----------------------------------------------------------------------------
// 
template<class GEN, class RECO>
ObjectMatcher<GEN,RECO>::ObjectMatcher( const edm::ParameterSet& pset ) :
  ObjectMatcherBase(pset)
{
  LogTrace("EnergyScale")
    << "[ObjectMatcher" << id() << "::" << __func__ << "]"
    << " Constructing...";  
}

// -----------------------------------------------------------------------------
// 
template<class GEN, class RECO>
ObjectMatcher<GEN,RECO>::ObjectMatcher() :
  ObjectMatcherBase()
{;}

// -----------------------------------------------------------------------------
// 
template<class GEN, class RECO>
ObjectMatcher<GEN,RECO>::~ObjectMatcher() {
  LogTrace("EnergyScale")
    << "[ObjectMatcher" << id() << "::" << __func__ << "]"
    << " Destructing...";  
}

// -----------------------------------------------------------------------------
//
template<class GEN, class RECO>
void ObjectMatcher<GEN,RECO>::gen( const edm::Event& event, 
				   const edm::EventSetup& setup, 
				   std::vector<HepLorentzVector>& objects ) {
  
  // Retrieve collection
  edm::Handle< edm::View<GEN> > collection;
  event.getByLabel( tags_.gen().tag(), collection );
  
  // Check if collection is found
  if ( !collection.isValid() ) {
    edm::LogError("EnergyScale")
      << "[ObjectMatcher" << id() << "::" << __func__ << "]"
      << " Invalid handle to \"" << collection.provenance()->className()
      << "\" with label \"" << collection.provenance()->moduleLabel() << "\"!";
    return;
  }
  
  // Check 
  if ( collection.failedToGet() ) {
    edm::LogError("EnergyScale")
      << "[ObjectMatcher" << id() << "::" << __func__ << "]"
      << " Invalid handle to object!" << std::endl
      << " Type                    : \"" << tags_.gen().type_ << "\"" << std::endl
      << " InputTag: Label         : \"" << tags_.gen().label_ << "\""
      << "           Instance      : \"" << tags_.gen().instance_ << "\""
      << "           Process       : \"" << tags_.gen().process_ << "\""
      << " Provenance: ClassName   : \"" << collection.provenance()->className() << "\""
      << "             ModuleLabel : \"" << collection.provenance()->moduleLabel() << "\"";
    return;
  }

  // Two "leptons" in z-dir (why?!)
  HepLorentzVector l1(0.,0.,1.,1.);
  HepLorentzVector l2(0.,0.,1.,1.); //@@ should be (0.,0.,-1.,1.) ??
  
  // pT threshold
  float pt_threshold = 20.;
  
  // Iterate through collection
  typename edm::View<GEN>::const_iterator ii = collection->begin(); 
  typename edm::View<GEN>::const_iterator jj = collection->end(); 
  for ( ; ii != jj; ++ii ) {
    
    // Only consider objects above energy threshold
    if ( ii->pt() < pt_threshold ) { continue; } //@@ configurable ??
    
    // Build Lorentz vector for each object
    HepLorentzVector vec( ii->px(), 
			  ii->py(), 
			  ii->pz(), 
			  ii->energy() );

    double dr1 = l1.deltaR(vec);
    double dr2 = l2.deltaR(vec);

    // Check jet not co-linear with either of "leptons"
    if ( dr1 < 1. || dr2 < 1. ) { continue; }
    
    // Push back first 'N' objects
    if ( objects.size() < 2 ) { objects.push_back(vec); } //@@ configurable ??
    
    if ( edm::isDebugEnabled() ) {
      std::stringstream ss;
      ss << "[ObjectMatcher" << id() << "::" << __func__ << "]"
	 << " Object"
	 << " #: " << static_cast<uint32_t>( ii - collection->begin() ) 
	 << std::endl
	 << "  e: " << ii->energy() << std::endl
	 << "  pt: " << ii->pt() << std::endl
	 << "  px: " << ii->px() << std::endl
	 << "  py: " << ii->py() << std::endl
	 << "  pz: " << ii->pz() << std::endl
	 << "  eta: " << ii->eta() << std::endl
	 << "  phi: " << ii->phi();
      LogTrace("EnergyScale") << ss.str();
    }
  }

  if ( edm::isDebugEnabled() ) {
    LogTrace("EnergyScale")
      << "[ObjectMatcher" << id() << "::" << __func__ << "]"
      << " Found " << objects.size() 
      << " objects with E > " << pt_threshold
      << " GeV!";
  }

}

// -----------------------------------------------------------------------------
//
template<class GEN, class RECO>
void ObjectMatcher<GEN,RECO>::reco( const edm::Event& event, 
				    const edm::EventSetup& setup,
				    std::vector<HepLorentzVector>& objects ) {
  
  // Get collection
  edm::Handle< edm::View<RECO> > collection;
  event.getByLabel( tags_.reco().tag(), collection );
  
  // Check if collection is found
  if ( !collection.isValid() ) {
    edm::LogError("EnergyScale")
      << "[ObjectMatcher" << id() << "::" << __func__ << "]"
      << " Invalid handle to \"" << collection.provenance()->className()
      << "\" with label \"" << collection.provenance()->moduleLabel() << "\"!";
    return;
  }

  // Check 
  if ( collection.failedToGet() ) {
    edm::LogError("EnergyScale")
      << "[ObjectMatcher" << id() << "::" << __func__ << "]"
      << " Invalid handle to object!" << std::endl
      << " Type                    : \"" << tags_.reco().type_ << "\"" << std::endl
      << " InputTag: Label         : \"" << tags_.reco().label_ << "\""
      << "           Instance      : \"" << tags_.reco().instance_ << "\""
      << "           Process       : \"" << tags_.reco().process_ << "\""
      << " Provenance: ClassName   : \"" << collection.provenance()->className() << "\""
      << "             ModuleLabel : \"" << collection.provenance()->moduleLabel() << "\"";
    return;
  }
  
  // pT threshold
  float pt_threshold = 0.;

  // Iterate through collection
  typename edm::View<RECO>::const_iterator ii = collection->begin(); 
  typename edm::View<RECO>::const_iterator jj = collection->end(); 
  for ( ; ii != jj; ++ii ) {
    
    // Only consider objects above energy threshold
    if ( ii->pt() < pt_threshold ) { continue; } //@@ configurable ??
    
    // Build Lorentz vector for each object
    HepLorentzVector vec( ii->px(), 
			  ii->py(), 
			  ii->pz(), 
			  ii->energy() );
    
    // Push back first 'N' objects
    if ( objects.size() < 1000. ) { objects.push_back(vec); } //@@ configurable ??
    
    if ( edm::isDebugEnabled() ) {
      std::stringstream ss;
      ss << "[ObjectMatcher" << id() << "::" << __func__ << "]"
	 << " Object"
	 << " #: " << static_cast<uint32_t>( ii - collection->begin() ) 
	 << std::endl
	 << "  e: " << ii->energy() << std::endl
	 << "  pt: " << ii->pt() << std::endl
	 << "  px: " << ii->px() << std::endl
	 << "  py: " << ii->py() << std::endl
	 << "  pz: " << ii->pz() << std::endl
	 << "  eta: " << ii->eta() << std::endl
	 << "  phi: " << ii->phi();
      LogTrace("EnergyScale") << ss.str();
    }
  }

  if ( edm::isDebugEnabled() ) {
    LogTrace("EnergyScale")
      << "[ObjectMatcher" << id() << "::" << __func__ << "]"
      << " Found " << objects.size() 
      << " objects with E > " << pt_threshold
      << " GeV!";
  }

}

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
template class ObjectMatcher<reco::GenJet,reco::CaloJet>;
template class ObjectMatcher<reco::GenJet,pat::Jet>;
