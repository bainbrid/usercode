#include "JetMETCorrections/JetPlusTrack/interface/ObjectMatcherBase.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "JetMETCorrections/JetPlusTrack/interface/EnergyScaleHistogrammer.h"
#include "JetMETCorrections/JetPlusTrack/interface/EnergyScaleHistograms.h"
#include "JetMETCorrections/JetPlusTrack/interface/ObjectTags.h"
#include "JetMETCorrections/JetPlusTrack/src/ObjectMatcher.cc"
#include <sstream>
#include <string>

// -----------------------------------------------------------------------------
// 
ObjectMatcherBase::ObjectMatcherBase( const edm::ParameterSet& pset ) 
  : tags_( ObjectTags(pset) )
{
  LogTrace("EnergyScale")
    << "[ObjectMatcherBase::"<<__func__<<"]"
    << " Constructing...";  
}

// -----------------------------------------------------------------------------
// 
ObjectMatcherBase::ObjectMatcherBase() 
  : tags_()
{;}

// -----------------------------------------------------------------------------
// 
ObjectMatcherBase::~ObjectMatcherBase() {
  LogTrace("EnergyScale")
    << "[ObjectMatcherBase::"<<__func__<<"]"
    << " Destructing...";  
}


// -----------------------------------------------------------------------------
// 
ObjectMatcherBase* ObjectMatcherBase::Instance( const edm::ParameterSet& pset ) {
  LogTrace("EnergyScale")
    << "[ObjectMatcherBase::"<<__func__<<"]"
    << " Creating instance...";  
  
  std::string g = pset.getParameter<std::string>( "GenObjectType" );
  std::string r = pset.getParameter<std::string>( "RecoObjectType" );
  
  if ( g == "GenJet" && 
       r == "CaloJet" ) { 
    return new ObjectMatcher<reco::GenJet,reco::CaloJet>(pset);
  } else if ( g == "GenJet" && 
	      r == "PatJet" ) { 
    return new ObjectMatcher<reco::GenJet,pat::Jet>(pset);
  } else {
    edm::LogError("EnergyScale")
      << "[ObjectMatcherBase::"<<__func__<<"]"
      << " Unexpected string value for ObjectType! "
      << " GenObjectType: \"" << g << "\""
      << " RecoObjectType: \"" << r << "\"";
    return 0;
  } 
  
}

// -----------------------------------------------------------------------------
//
void ObjectMatcherBase::analyze( const edm::Event& event, 
				 const edm::EventSetup& setup ) {
  
  // Create 4-vector pairs
  std::vector<LorentzVectorPair> pairs;
  
  // Build vector of suitable gen objects:
  std::vector<HepLorentzVector> gen_objects;
  gen( event, setup, gen_objects );
  
  // Check if gen objects were found
  if ( gen_objects.empty() ) { 
    edm::LogWarning("EnergyScale")
      << "[ObjectMatcherBase::" << __func__ << "]"
      << " No suitable GenObjects found to generate 4-vectors!";
    return; 
  } else {
    LogTrace("EnergyScale")
      << "[ObjectMatcherBase::" << __func__ << "]"
      << " Found " << gen_objects.size() << " GenObjects!";
  }
  
  // Build vector of suitable reco objects:
  std::vector<HepLorentzVector> reco_objects;
  reco( event, setup, reco_objects );
  
  // Check if reco objects were found
  if ( reco_objects.empty() ) { 
    edm::LogWarning("EnergyScale")
      << "[ObjectMatcherBase::" << __func__ << "]"
      << " No suitable RecoObjects found to generate 4-vectors!";
    return; 
  } else {
    LogTrace("EnergyScale")
      << "[ObjectMatcherBase::" << __func__ << "]"
      << " Found " << reco_objects.size() << " RecoObjects!";
  }
  
  // Record of matched reco objects
  std::vector< std::vector<HepLorentzVector>::const_iterator > matched;
  
  // Iterate through 4-vectors of gen objects
  std::vector<HepLorentzVector>::const_iterator ig = gen_objects.begin();
  std::vector<HepLorentzVector>::const_iterator jg = gen_objects.end();
  for ( ; ig != jg; ++ig ) {
    
    // Create 4-vector pair
    LorentzVectorPair match( *ig );
    
    // Iterate through 4-vectors of reco objects and find match
    double dr = 1.e6;
    std::vector<HepLorentzVector>::const_iterator iter = reco_objects.end(); 
    std::vector<HepLorentzVector>::const_iterator ir = reco_objects.begin(); 
    std::vector<HepLorentzVector>::const_iterator jr = reco_objects.end(); 
    for ( ; ir != jr; ++ir ) {
      //       if ( std::find( matched.begin(), 
      // 		      matched.end(), 
      // 		      ir ) == matched.end() ) { // Check not already matched
	double idr = ig->deltaR(*ir);
	if ( idr < dr ) { // Check if closer in dR
	  iter = ir;
	  dr = idr;
	}
	//       }
    }
    
    // Check if match is found 
    if ( iter != jr ) { 
      
      matched.push_back(iter);
      match.reco(*iter);
      
      // Debug
      if ( edm::isDebugEnabled() ) {
	std::stringstream ss;
	ss << "[ObjectMatcherBase::" << __func__ << "]"
	   << " Matched GenObject #" << static_cast<uint32_t>( ig - gen_objects.begin() )
	   << " to RecoObject #" << static_cast<uint32_t>( iter - reco_objects.begin() ) << std::endl
	   << "  GenObject:  pt/eta/phi: " 
	   << ig->perp() << "/" 
	   << ig->eta() << "/"
	   << ig->phi() << std::endl
	   << "  RecoObject: pt/eta/phi: " 
	   << iter->perp() << "/" 
	   << iter->eta() << "/"
	   << iter->phi() << std::endl
	   << "  DR: " << dr << std::endl;
	LogTrace("EnergyScale") << ss.str();
      }

    }
    
    pairs.push_back(match);
    
  }
  
  // Access to histos object
  EnergyScaleHistograms* histos = EnergyScaleHistogrammer::Histograms( tags_ ); 
  if ( histos ) { histos->analyze( pairs ); }
  else {
    edm::LogError("EnergyScale")
      << "[ObjectMatcherBase::" << __func__ << "]"
      << " NULL pointer to EnergyScaleHistograms object!";
    return;
  }
  
}

// -----------------------------------------------------------------------------
//
std::string ObjectMatcherBase::id() {
  return std::string( "<" + tags_.gen().type_ + "," + tags_.reco().type_ + ">" );
}
