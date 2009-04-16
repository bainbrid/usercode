#include "bainbrid/Test/test/PATPhotonIDProducer.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <memory>
#include <vector>

// -----------------------------------------------------------------------------
//
pat::PATPhotonIDProducer::PATPhotonIDProducer( const edm::ParameterSet& pset ) 
  : algo_(),
    phoLabel_( pset.getParameter<edm::InputTag>("PhotonTag") )
{
  algo_.setup( pset );
  produces< std::vector<pat::Photon> >();
}

// -----------------------------------------------------------------------------
//
void pat::PATPhotonIDProducer::produce( edm::Event& event, 
					const edm::EventSetup& setup ) {

  // Retrieve photons from Event
  edm::Handle< std::vector<pat::Photon> > orig_photons;
  event.getByLabel( phoLabel_, orig_photons );
  
  // Create copy of photon collection
  std::vector<pat::Photon>* new_photons = new std::vector<pat::Photon>();
  new_photons->resize( orig_photons->size() );
  std::copy( orig_photons->begin(),
	     orig_photons->end(), 
	     new_photons->begin() );
  
  // Iterate through new photon collection
  std::vector<pat::Photon>::iterator iphoton = new_photons->begin(); 
  std::vector<pat::Photon>::iterator jphoton = new_photons->end(); 
  for ( ; iphoton != jphoton; ++iphoton ) {

    // Retrieve internal ID object of new photon
    const reco::PhotonID* id = iphoton->photonID();
    if ( id ) {
      
      // Make copy of ID with modified ECAL isolation value
      reco::PhotonID modified_id( id->isLooseEM(),
				  id->isLoosePhoton(),
				  id->isTightPhoton(),
				  id->isolationSolidTrkCone(),
				  id->isolationHollowTrkCone(), // iphoton->trackIso(), //@@ take TRK isolation value from PAT
				  id->nTrkSolidCone(),
				  id->nTrkHollowCone(), 
				  iphoton->ecalIso(), //@@ take ECAL isolation value from PAT
				  iphoton->hcalIso(), //@@ take HCAL isolation value from PAT
				  id->r9(),
				  id->isEBPho(),
				  id->isEEPho(),
				  id->isEBGap(),
				  id->isEEGap(),
				  id->isEBEEGap(),
				  id->isAlsoElectron() );
      
      // Set ID according to thresholds on new ECAL isolation value
      if ( modified_id.isEBPho() ) { algo_.decideEB( modified_id, &*iphoton ); }
      else if ( modified_id.isEEPho() ) { algo_.decideEE( modified_id, &*iphoton ); }
      else { std::cout << "PROBLEM!" << std::endl; }
      
      // Overwrite internal ID with newly constructed one
      iphoton->setPhotonID( modified_id );
      
    }
    
  }    

  // Write new photons to Event
  std::auto_ptr< std::vector<pat::Photon> > photons( new_photons );
  event.put( photons );
  
}

// -----------------------------------------------------------------------------
// 
#include "FWCore/Framework/interface/MakerMacros.h"
using namespace pat;
DEFINE_FWK_MODULE(PATPhotonIDProducer);