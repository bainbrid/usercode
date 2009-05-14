#include "bainbrid/Test/test/ObjectMatcher.h"
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
  if ( !collection.isValid() || collection.failedToGet() ) {
    std::stringstream ss;
    ss << "[ObjectMatcher" << id() << "::" << __func__ << "]"
       << " Unable to retrieve collection from event!" << std::endl
       << " Type     : \"" << tags_.gen().type_ << "\"" << std::endl
       << " Label    : \"" << tags_.gen().label_ << "\"" << std::endl
       << " Instance : \"" << tags_.gen().instance_ << "\"" << std::endl
       << " Process  : \"" << tags_.gen().process_ << "\"";
    edm::LogError("EnergyScale") << ss.str();
    return;
  } else {
    if ( edm::isDebugEnabled() && verbose_ ) {
      std::stringstream ss;
      ss << "[ObjectMatcher" << id() << "::" << __func__ << "]"
	 << " Retrieved collection from event:" << std::endl
	 << " ClassName           : \"" << collection.provenance()->className() << "\"" << std::endl
	 << " FriendlyClassName   : \"" << collection.provenance()->friendlyClassName() << "\"" << std::endl
	 << " ModuleName          : \"" << collection.provenance()->moduleName() << "\"" << std::endl
	 << " ModuleLabel         : \"" << collection.provenance()->moduleLabel() << "\"" << std::endl
	 << " ProductInstanceName : \"" << collection.provenance()->productInstanceName() << "\"" << std::endl
	 << " ProcessName         : \"" << collection.provenance()->processName() << "\"" << std::endl
	 << " CollectionSize      : "   << collection->size();
      LogTrace("EnergyScale") << ss.str();
    }      
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

//     std::cout << "ALL GenJets:"
// 	      << " Event: " << event.id().event()
// 	      << " GenJet# " << int( ii - collection->begin() )
// 	      << " e= " << ii->energy()
// 	      << " et= " << ii->et()
// 	      << " pt= " << ii->pt()
// 	      << " px= " << ii->px()
// 	      << " py= " << ii->py()
// 	      << " pz= " << ii->pz()
// 	      << std::endl;
    
//     if ( edm::isDebugEnabled() && verbose_ ) {
//       std::stringstream ss;
//       ss << "[ObjectMatcher" << id() << "::" << __func__ << "]"
// 	 << " Event: " << event.id().event()
// 	 << " Object# " << static_cast<uint32_t>( ii - collection->begin() ) 
// 	 << std::endl
// 	 << " e: " << ii->energy() << std::endl
// 	 << " et: " << ii->et()
// 	 << " pt: " << ii->pt()
// 	 << " px: " << ii->px()
// 	 << " py: " << ii->py()
// 	 << " pz: " << ii->pz()
// 	 << " eta: " << ii->eta()
// 	 << " phi: " << ii->phi()
// 	 << std::endl 
// 	 << ii->print();
//       LogTrace("EnergyScale") << ss.str();
//     }
    
//     // Only consider objects above energy threshold
//     if ( ii->pt() < pt_threshold ) { continue; } //@@ configurable ??
    
    // Build Lorentz vector for each object
    HepLorentzVector vec( ii->px(), 
			  ii->py(), 
			  ii->pz(), 
			  ii->energy() );
    
    // Check jet not co-linear with either of "leptons"
    if ( l1.deltaR(vec) < 1. || 
	 l2.deltaR(vec) < 1. ) { continue; }
    
    // Push back first 'N' objects
    if ( objects.size() < 2 ) { objects.push_back(vec); } //@@ configurable ??
    
  }

  if ( edm::isDebugEnabled() && verbose_ ) {
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
  if ( !collection.isValid() || collection.failedToGet() ) {
    std::stringstream ss;
    ss << "[ObjectMatcher" << id() << "::" << __func__ << "]"
       << " Unable to retrieve collection from event!" << std::endl
       << " Type     : \"" << tags_.reco().type_ << "\"" << std::endl
       << " Label    : \"" << tags_.reco().label_ << "\"" << std::endl
       << " Instance : \"" << tags_.reco().instance_ << "\"" << std::endl
       << " Process  : \"" << tags_.reco().process_ << "\"";
    edm::LogError("EnergyScale") << ss.str();
    return;
  } else {
    if ( edm::isDebugEnabled() && verbose_ ) {
      std::stringstream ss;
      ss << "[ObjectMatcher" << id() << "::" << __func__ << "]"
	 << " Retrieved collection from event:" << std::endl
	 << " ClassName           : \"" << collection.provenance()->className() << "\"" << std::endl
	 << " FriendlyClassName   : \"" << collection.provenance()->friendlyClassName() << "\"" << std::endl
	 << " ModuleName          : \"" << collection.provenance()->moduleName() << "\"" << std::endl
	 << " ModuleLabel         : \"" << collection.provenance()->moduleLabel() << "\"" << std::endl
	 << " ProductInstanceName : \"" << collection.provenance()->productInstanceName() << "\"" << std::endl
	 << " ProcessName         : \"" << collection.provenance()->processName() << "\"" << std::endl
	 << " CollectionSize      : "   << collection->size();
      LogTrace("EnergyScale") << ss.str();
    }      
  }
  
  // pT threshold
  float pt_threshold = 0.;

  // Iterate through collection
  typename edm::View<RECO>::const_iterator ii = collection->begin(); 
  typename edm::View<RECO>::const_iterator jj = collection->end(); 
  for ( ; ii != jj; ++ii ) {

//     std::cout << "ALL RecoJets:"
// 	      << " Event: " << event.id().event()
// 	      << " RecoJet# " << int( ii - collection->begin() )
// 	      << " e= " << ii->energy()
// 	      << " et= " << ii->et()
// 	      << " pt= " << ii->pt()
// 	      << " px= " << ii->px()
// 	      << " py= " << ii->py()
// 	      << " pz= " << ii->pz()
// 	      << std::endl;
    
    if ( edm::isDebugEnabled() && verbose_ ) {
      std::stringstream ss;
      ss << "[ObjectMatcher" << id() << "::" << __func__ << "]"
	 << " Event: " << event.id().event()
	 << " Object# " << static_cast<uint32_t>( ii - collection->begin() ) 
	 << std::endl
// 	 << " e: " << ii->energy() << std::endl
// 	 << " et: " << ii->et()
// 	 << " pt: " << ii->pt()
// 	 << " px: " << ii->px()
// 	 << " py: " << ii->py()
// 	 << " pz: " << ii->pz()
// 	 << " eta: " << ii->eta()
// 	 << " phi: " << ii->phi()
// 	 << std::endl 
	 << ii->print();
      LogTrace("EnergyScale") << ss.str();
    }

    // Only consider objects above energy threshold
    if ( ii->pt() < pt_threshold ) { continue; } //@@ configurable ??
    
    // Build Lorentz vector for each object
    HepLorentzVector vec( ii->px(),
			  ii->py(),
			  ii->pz(),
			  ii->energy() );
    
    // Push back first 'N' objects
    if ( objects.size() < 100 ) { objects.push_back(vec); } //@@ configurable ??
    
  }

  if ( edm::isDebugEnabled() && verbose_ ) {
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
