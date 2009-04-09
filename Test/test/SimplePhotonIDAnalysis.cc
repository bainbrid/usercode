#include "bainbrid/Test/test/SimplePhotonIDAnalysis.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"
#include "TH2D.h"
#include <memory>
#include <cmath>

// -----------------------------------------------------------------------------
//
SimplePhotonIDAnalysis::SimplePhotonIDAnalysis( const edm::ParameterSet& pset ) 
  : histos_(),
    labels_(),
    photons_( pset.getUntrackedParameter<edm::InputTag>("Photons") ),
    others_( pset.getUntrackedParameter<edm::InputTag>("OtherPhotons") )
{
  labels_.push_back("SignalPhotons");
  labels_.push_back("OtherPhotons");
  labels_.push_back("FakePhotons");
}

// -----------------------------------------------------------------------------
//
void SimplePhotonIDAnalysis::analyze( const edm::Event& event, 
				      const edm::EventSetup& setup ) {
  
  // Photon counter
  std::vector<uint16_t> multiplicity( labels_.size(), 0 );
  
  // Retrieve PAT photon collection from Event
  edm::Handle<edm::View<pat::Photon> > photons;
  event.getByLabel( photons_, photons );
  
  // Retrieve PAT photonID collection from Event
  edm::Handle<edm::View<pat::Photon> > others;
  event.getByLabel( others_, others );
  
  // Iterate through PAT photons
  edm::View<pat::Photon>::const_iterator iphoton = photons->begin(); 
  edm::View<pat::Photon>::const_iterator jphoton = photons->end(); 
  for ( ; iphoton != jphoton; ++iphoton ) {
    
    // GMSB PhotonID based on MC truth
    bool photon_id = false; 
    bool gmsb_photon_id = false; 
    reco::Particle* part = const_cast<reco::Particle*>( iphoton->genPhoton() );
    reco::Candidate* cand = dynamic_cast<reco::Candidate*>( part );
    if ( cand && cand->pdgId() == 22 ) { 
      photon_id = true; 
      const reco::Candidate* mother = cand->mother();
      if ( mother && cand->pdgId() == cand->mother()->pdgId() ) { mother = mother->mother(); }
      if ( mother && mother->pdgId() == 1000022 ) { gmsb_photon_id = true; }
    }
    
    // Some selection
    //if ( iphoton->et() < 100. || fabs( iphoton->eta() ) > 2. ) { continue; }

    // Iterate though categories
    for ( uint16_t ii = 0; ii < labels_.size(); ++ii ) {
      
      if ( !( ( labels_[ii] == "SignalPhotons" && gmsb_photon_id ) || 
	      ( labels_[ii] == "OtherPhotons" && photon_id && !gmsb_photon_id ) ||
	      ( labels_[ii] == "FakePhotons" && !photon_id ) ) ) { continue; }
	   
      multiplicity[ii]++; // Multiplicity

      std::string name = labels_[ii];
      
      // Photon ID flags
      histos_[name+"_PhotonID"]->Fill( 0. );
      if ( iphoton->isLooseEM() )     { histos_[name+"_PhotonID"]->Fill( 1. ); }
      if ( iphoton->isLoosePhoton() ) { histos_[name+"_PhotonID"]->Fill( 2. ); }
      if ( iphoton->isTightPhoton() ) { histos_[name+"_PhotonID"]->Fill( 3. ); }
      if ( !iphoton->isLooseEM() &&
	   !iphoton->isLoosePhoton() &&
	   !iphoton->isTightPhoton() ) { histos_[name+"_PhotonID"]->Fill( 4. ); }
      
      // Kinematics
      histos_[name+"_Energy"]->Fill( iphoton->energy() );
      histos_[name+"_Et"]->Fill( iphoton->et() );
      histos_[name+"_Eta"]->Fill( iphoton->eta() );
      
      // Track isolation from RECO
      histos_[name+"_NTrkSolid"]->Fill( iphoton->nTrkSolidCone() );
      histos_[name+"_NTrkHollow"]->Fill( iphoton->nTrkHollowCone() );
      histos_[name+"_IsoSolidTrk"]->Fill( iphoton->isolationSolidTrkCone() );
      histos_[name+"_IsoHollowTrk"]->Fill( iphoton->isolationHollowTrkCone() );
      
      // ECAL/HCAL isolation from RECO
      histos_[name+"_IsoEcalRecHit"]->Fill( iphoton->isolationEcalRecHit() );
      histos_[name+"_IsoHcalRecHit"]->Fill( iphoton->isolationHcalRecHit() );
      histos_[name+"_HadronicOverEm"]->Fill( iphoton->hadronicOverEm() );
      if ( fabs( iphoton->isolationEcalRecHit() ) > 0. ) {
	histos_[name+"_IsoHcalRecHitOverIsoEcalRecHit"]->Fill( iphoton->isolationHcalRecHit() / 
							       iphoton->isolationEcalRecHit() );
      }
      histos_[name+"_R9"]->Fill( iphoton->r9() );
      
      // Isolation from PAT
      histos_[name+"_IsoEcal"]->Fill( iphoton->ecalIso() );
      histos_[name+"_IsoHcal"]->Fill( iphoton->hcalIso() );
      histos_[name+"_IsoTrk"]->Fill( iphoton->trackIso() );
      histos_[name+"_IsoAll"]->Fill( iphoton->caloIso() );
      histos_[name+"_HcalIsoOverEcalIso"]->Fill( iphoton->hcalIso() / iphoton->ecalIso() );

      // Iterate through OtherPhotonID
      if ( !others_.label().empty() ) {

	edm::View<pat::Photon>::const_iterator iother = others->begin(); 
	edm::View<pat::Photon>::const_iterator jother = others->end(); 
	for ( ; iother != jother; ++iother ) {

	  if ( uint32_t( iphoton - photons->begin() ) !=
	       uint32_t( iother - others->begin() ) ) { continue; }

	  histos2d_[name+"_CompareEcalIso"]->Fill( iphoton->isolationEcalRecHit(),
						   iother->isolationEcalRecHit() );
	  histos2d_[name+"_CompareHcalIso"]->Fill( iphoton->isolationHcalRecHit(),
						   iother->isolationHcalRecHit() );
	  histos2d_[name+"_CompareTrkIso"]->Fill( iphoton->isolationHollowTrkCone(),
						  iother->isolationHollowTrkCone() );
	  
	  float iii = 0.;
	  if ( iphoton->isTightPhoton() ) { iii = 3.; }
	  else if ( iphoton->isLoosePhoton() ) { iii = 2.; }
	  else if ( iphoton->isLooseEM() ) { iii = 1.; }
	  else { iii = 0.; }

	  float jjj = 0.;
	  if ( iother->isTightPhoton() ) { jjj = 3.; }
	  else if ( iother->isLoosePhoton() ) { jjj = 2.; }
	  else if ( iother->isLooseEM() ) { jjj = 1.; }
	  else { jjj = 0.; }

	  histos_[name+"_OtherPhotonID"]->Fill( jjj );
	  histos2d_[name+"_ComparePhotonID"]->Fill( iii, jjj );
	  
	}
      }

    }
    
  } // photon loop    
  
  // Photon multiplicities
  histos_["Multiplicity"]->Fill( photons->size() );
  for ( uint16_t ii = 0; ii < labels_.size(); ++ii ) { 
    std::string name = labels_[ii];
    histos_[labels_[ii]+"_Multiplicity"]->Fill( multiplicity[ii] );
  }
  
}

// -----------------------------------------------------------------------------
//
void SimplePhotonIDAnalysis::beginJob( const edm::EventSetup& ) {

  edm::Service<TFileService> fs;

  histos_["Multiplicity"] = fs->make<TH1D>("Multiplicity","Multiplicity",21,-0.5,20.5);
  
  for ( uint16_t ii = 0; ii < 3; ++ii ) {
    
    std::string name = labels_[ii]; // SignalPhotons, OtherPhotons, FakePhotons

    TFileDirectory dir = fs->mkdir(name);

    histos_[name+"_Multiplicity"] = dir.make<TH1D>("Multiplicity","Multiplicity",21,-0.5,20.5);
    
    // Kinematics
    histos_[name+"_Energy"] = dir.make<TH1D>("Energy",";E / GeV;",1000,0,1000.);
    histos_[name+"_Et"] = dir.make<TH1D>("Et",";E_{T} / GeV;",1000,0,1000.);
    histos_[name+"_Eta"] = dir.make<TH1D>("Eta",";#eta;",800,-4.,4.);
  
    // Track isolation from RECO
    histos_[name+"_NTrkSolid"] = dir.make<TH1D>("nTrkSolid",";N_{trk};",101,-0.5,100.5);
    histos_[name+"_NTrkHollow"] = dir.make<TH1D>("nTrkHollow",";N_{trk};",101,-0.5,100.5);
    histos_[name+"_IsoSolidTrk"] = dir.make<TH1D>("IsoSolidTrk",";#sum E_{T} [GeV];",300,0.,300.);
    histos_[name+"_IsoHollowTrk"] = dir.make<TH1D>("IsoHollowTrk",";#sum E_{T} [GeV];",300,0.,300.);

    // ECAL/HCAL isolation from RECO
    histos_[name+"_IsoEcalRecHit"] = dir.make<TH1D>("IsoEcalRecHit",";#sum E_{T} [GeV];",300,0.,300.);
    histos_[name+"_IsoHcalRecHit"] = dir.make<TH1D>("IsoHcalRecHit",";#sum E_{T} [GeV];",300,0.,300.);
    histos_[name+"_HadronicOverEm"] = dir.make<TH1D>("HadronicOverEm","HCAL/ECAL",100,0.,2.);
    histos_[name+"_IsoHcalRecHitOverIsoEcalRecHit"] = 
      dir.make<TH1D>("IsoHcalRecHitOverIsoEcalRecHit","HCAL/ECAL",100,0.,2.);
    histos_[name+"_R9"] = dir.make<TH1D>("R9",";R9;",100,0.,3.);
  
    // Isolation from PAT
    histos_[name+"_IsoEcal"] = dir.make<TH1D>("IsoEcal",";#sum E_{T} [GeV];",300,0.,300.);
    histos_[name+"_IsoHcal"] = dir.make<TH1D>("IsoHcal",";#sum E_{T} [GeV];",300,0.,300.);
    histos_[name+"_IsoTrk"] = dir.make<TH1D>("IsoTrk",";#sum E_{T} [GeV];",300,0.,300.);
    histos_[name+"_IsoAll"] = dir.make<TH1D>("IsoAll",";#sum E_{T} [GeV];",300,0.,300.);
    histos_[name+"_HcalIsoOverEcalIso"] = dir.make<TH1D>("HcalIsoOverEcalIso",";HCAL/ECAL;",200,0.,2.);
    
    // Photon ID flags
    histos_[name+"_PhotonID"] = dir.make<TH1D>("PhotonID","",5,-0.5,4.5);
    if ( !others_.label().empty() ) { 
      histos_[name+"_OtherPhotonID"] = dir.make<TH1D>("OtherPhotonID","",5,-0.5,4.5);
      histos2d_[name+"_ComparePhotonID"] = dir.make<TH2D>("ComparePhotonID","ComparePhotonID",4,-0.5,3.5,4,-0.5,3.5);
      histos2d_[name+"_CompareEcalIso"] = dir.make<TH2D>("CompareEcalIso","CompareEcalIso",300,0.,300.,300,0.,300.);
      histos2d_[name+"_CompareHcalIso"] = dir.make<TH2D>("CompareHcalIso","CompareHcalIso",300,0.,300.,300,0.,300.);
      histos2d_[name+"_CompareTrkIso"] = dir.make<TH2D>("CompareTrkIso","CompareTrkIso",300,0.,300.,300,0.,300.);
    }
    
  } 
  
}

// -----------------------------------------------------------------------------
// Define class to be a plug-in module
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SimplePhotonIDAnalysis);



//iphoton->ecalIsoDeposit()->begin()->dR();


//       if ( iphoton->isLooseEM() ) { histos2d_["PhotonID_Old"]->Fill( iphoton->isolationEcalRecHit(), 1. ); }
//       if ( iphoton->isLoosePhoton() ) { histos2d_["PhotonID_Old"]->Fill( iphoton->isolationEcalRecHit(), 2. ); }
//       if ( iphoton->isTightPhoton() ) { histos2d_["PhotonID_Old"]->Fill( iphoton->isolationEcalRecHit(), 3. ); }
//       if ( !iphoton->isLooseEM() &&
// 	   !iphoton->isLoosePhoton() &&
// 	   !iphoton->isTightPhoton() ) { histos2d_["PhotonID_Old"]->Fill( iphoton->isolationEcalRecHit(), 0. ); }
      
//       if ( iphotonid->isLooseEM() ) { histos2d_["PhotonID_New"]->Fill( iphotonid->isolationEcalRecHit(), 1. ); }
//       if ( iphotonid->isLoosePhoton() ) { histos2d_["PhotonID_New"]->Fill( iphotonid->isolationEcalRecHit(), 2. ); }
//       if ( iphotonid->isTightPhoton() ) { histos2d_["PhotonID_New"]->Fill( iphotonid->isolationEcalRecHit(), 3. ); }
//       if ( !iphotonid->isLooseEM() &&
// 	   !iphotonid->isLoosePhoton() &&
// 	   !iphotonid->isTightPhoton() ) { histos2d_["PhotonID_Old"]->Fill( iphotonid->isolationEcalRecHit(), 0. ); }

//     // Some selection criteria (based on RECO!)
//     if ( !( iphoton->et() > 10. && 
// 	    iphoton->eta() > -3.0 &&
// 	    iphoton->eta() < 3.0 && 
// 	    iphoton->r9() > 0. &&
// 	    iphoton->isolationHcalRecHit() / iphoton->isolationEcalRecHit() < 100. )
// 	 ) { continue; }


//     if ( iphoton->genPhoton() &&                                // check pointer
//     	 iphoton->genPhoton()->mother() &&                      // check mother pointer
//     	 iphoton->genPhoton()->mother()->pdgId() == 1000022 ) { // check if mother is neutralino
//       gmsb_photon_id = true;
//     }

//     std::cout << " TEST " 
// 	      << " SignalPhoton: " << gmsb_photon_id 
// 	      << " TruthPhoton: " << photon_id 
// 	      << " FakePhoton: " << !photon_id
// 	      << " part " << part
// 	      << " cand " << cand
// 	      << " cand->mother() " 
// 	      << ( cand ? (uint32_t)cand->mother() : 10000001 )
// 	      << " cand->mother()->pdgId() " 
// 	      << ( cand ? ( (uint32_t)cand->mother() ? cand->mother()->pdgId() : 10000002 ) : 10000003 )
// 	      << std::endl;
      
