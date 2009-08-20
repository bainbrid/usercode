#include "TestPatCrossCleaner.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1.h"
#include <algorithm>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cmath>

// -----------------------------------------------------------------------------
//
TestPatCrossCleaner::TestPatCrossCleaner( const edm::ParameterSet& pset ) 
  : //
    rawPhotons_( pset.getParameter<edm::InputTag>("RawPhotons") ),
    rawJets_( pset.getParameter<edm::InputTag>("RawJets") ),
    rawMuons_( pset.getParameter<edm::InputTag>("RawMuons") ),
    rawElectrons_( pset.getParameter<edm::InputTag>("RawElectrons") ),
    //
    ccPhotons_( pset.getParameter<edm::InputTag>("CCPhotons") ),
    ccJets_( pset.getParameter<edm::InputTag>("CCJets") ),
    ccMuons_( pset.getParameter<edm::InputTag>("CCMuons") ),
    ccElectrons_( pset.getParameter<edm::InputTag>("CCElectrons") ),
    //
    droppedPhotons_( pset.getParameter<edm::InputTag>("DroppedPhotons") ),
    droppedJets_( pset.getParameter<edm::InputTag>("DroppedJets") ),
    droppedMuons_( pset.getParameter<edm::InputTag>("DroppedMuons") ),
    droppedElectrons_( pset.getParameter<edm::InputTag>("DroppedElectrons") ),
    //
    jetMatchedCCPhotons_( pset.getParameter<edm::InputTag>("JetMatchedCCPhotons") ),
    jetMatchedCCElectrons_( pset.getParameter<edm::InputTag>("JetMatchedCCElectrons") ),
    jetMatchedCCPartons_( pset.getParameter<edm::InputTag>("JetMatchedCCPartons") ),
    //
    jetMatchedDroppedPhotons_( pset.getParameter<edm::InputTag>("JetMatchedDroppedPhotons") ),
    jetMatchedDroppedElectrons_( pset.getParameter<edm::InputTag>("JetMatchedDroppedElectrons") ),
    jetMatchedDroppedPartons_( pset.getParameter<edm::InputTag>("JetMatchedDroppedPartons") ),
    //
    photonMatchedCCPhotons_( pset.getParameter<edm::InputTag>("PhotonMatchedCCPhotons") ),
    photonMatchedDroppedPhotons_( pset.getParameter<edm::InputTag>("PhotonMatchedDroppedPhotons") ),
    //
    histos_(),
    histos2d_()
{;}

// -----------------------------------------------------------------------------
//
TH1* TestPatCrossCleaner::histo( const std::string& histogram_name ) {
  std::map<std::string,TH1D*>::const_iterator ii = histos_.find(histogram_name);
  if ( ii != histos_.end() ) { return ii->second; }
  else {
    std::map<std::string,TH2D*>::const_iterator jj = histos2d_.find(histogram_name);
    if ( jj != histos2d_.end() ) { return jj->second; }
  }
  edm::LogWarning("TEST") << "Cannot find string: " << histogram_name;
  return 0;
}

// -----------------------------------------------------------------------------
//
void TestPatCrossCleaner::analyze( const edm::Event& iEvent, 
				   const edm::EventSetup& iSetup ) {


  static const float em_fraction = 0.95;


  // -------------------- Get collections --------------------


  edm::Handle< edm::View<pat::Photon> > raw_photons;
  iEvent.getByLabel(rawPhotons_,raw_photons); 
  edm::Handle< edm::View<pat::Jet> > raw_jets;
  iEvent.getByLabel(rawJets_,raw_jets); 
  edm::Handle< edm::View<pat::Electron> > raw_electrons;
  iEvent.getByLabel(rawElectrons_,raw_electrons); 
  edm::Handle< edm::View<pat::Muon> > raw_muons;
  iEvent.getByLabel(rawMuons_,raw_muons); 

  edm::Handle< edm::View<pat::Photon> > cc_photons;
  iEvent.getByLabel(ccPhotons_,cc_photons); 
  edm::Handle< edm::View<pat::Jet> > cc_jets;
  iEvent.getByLabel(ccJets_,cc_jets); 
  edm::Handle< edm::View<pat::Electron> > cc_electrons;
  iEvent.getByLabel(ccElectrons_,cc_electrons); 
  edm::Handle< edm::View<pat::Muon> > cc_muons;
  iEvent.getByLabel(ccMuons_,cc_muons); 

  edm::Handle< edm::View<pat::Photon> > dropped_photons;
  iEvent.getByLabel(droppedPhotons_,dropped_photons); 
  edm::Handle< edm::View<pat::Jet> > dropped_jets;
  iEvent.getByLabel(droppedJets_,dropped_jets); 
  edm::Handle< edm::View<pat::Electron> > dropped_electrons;
  iEvent.getByLabel(droppedElectrons_,dropped_electrons); 
  edm::Handle< edm::View<pat::Muon> > dropped_muons;
  iEvent.getByLabel(droppedMuons_,dropped_muons); 

  edm::Handle< edm::Association< reco::GenParticleCollection > > jet_matched_cc_photons;
  iEvent.getByLabel( jetMatchedCCPhotons_, jet_matched_cc_photons ); 
  edm::Handle< edm::Association< reco::GenParticleCollection > > jet_matched_cc_electrons;
  iEvent.getByLabel( jetMatchedCCElectrons_, jet_matched_cc_electrons ); 
  edm::Handle< edm::Association< reco::GenParticleCollection > > jet_matched_cc_partons;
  iEvent.getByLabel( jetMatchedCCPartons_, jet_matched_cc_partons ); 

  edm::Handle< edm::Association< reco::GenParticleCollection > > jet_matched_dropped_photons;
  iEvent.getByLabel( jetMatchedDroppedPhotons_, jet_matched_dropped_photons ); 
  edm::Handle< edm::Association< reco::GenParticleCollection > > jet_matched_dropped_electrons;
  iEvent.getByLabel( jetMatchedDroppedElectrons_, jet_matched_dropped_electrons ); 
  edm::Handle< edm::Association< reco::GenParticleCollection > > jet_matched_dropped_partons;
  iEvent.getByLabel( jetMatchedDroppedPartons_, jet_matched_dropped_partons ); 

  edm::Handle< edm::Association< reco::GenParticleCollection > > photon_matched_cc_photons;
  iEvent.getByLabel( photonMatchedCCPhotons_, photon_matched_cc_photons ); 
  edm::Handle< edm::Association< reco::GenParticleCollection > > photon_matched_dropped_photons;
  iEvent.getByLabel( photonMatchedDroppedPhotons_, photon_matched_dropped_photons ); 
  
  if ( !raw_photons.isValid() || 
       !raw_jets.isValid() || 
       !raw_electrons.isValid() || 
       !raw_muons.isValid() || 
       !cc_photons.isValid() || 
       !cc_jets.isValid() || 
       !cc_electrons.isValid() || 
       !cc_muons.isValid() ||
       !dropped_photons.isValid() || 
       !dropped_jets.isValid() || 
       !dropped_electrons.isValid() || 
       !dropped_muons.isValid() ||
       !jet_matched_cc_photons.isValid() ||
       !jet_matched_dropped_photons.isValid() ||
       !photon_matched_cc_photons.isValid() ||
       !photon_matched_dropped_photons.isValid() ) {
    edm::LogWarning("TEST") << " Missing collection!";
    return;
  }
  
      
  // -------------------- Cross-cleaned jets --------------------


  { 
    histo("CCJets_Num")->Fill( cc_jets->size() ); 
    edm::View<pat::Jet>::const_iterator ii = cc_jets->begin();
    edm::View<pat::Jet>::const_iterator jj = cc_jets->end();
    for ( ; ii != jj; ++ii ) {
      histo("CCJets_EMfraction")->Fill( ii->emEnergyFraction() ); 
      if ( ii->emEnergyFraction() > em_fraction ) {

	// Nearest GenPhoton
	{
	  
	  float dr1 = -1.;
	  float et1 = -1.;
	  uint32_t index1 = ii - cc_jets->begin();
	  edm::RefToBase<pat::Jet> jet_ref = cc_jets->refAt(index1);
	  reco::GenParticleRef gen_jet_ref;
	  if ( jet_ref.isNonnull() ) { 
	    gen_jet_ref = (*jet_matched_cc_photons)[jet_ref];
	    if ( gen_jet_ref.isNonnull() ) { 
	      reco::GenParticle* gen_particle = const_cast<reco::GenParticle*>( &(*gen_jet_ref) );
	      reco::Candidate* gen_candidate = dynamic_cast<reco::Candidate*>( gen_particle );
	      if ( gen_particle ) {
		dr1 = reco::deltaR( ii->p4().eta(), gen_candidate->p4().eta(), 
				    ii->p4().phi(), gen_candidate->p4().phi() );
		et1 = fabs( ii->et() - gen_candidate->et() ) / gen_candidate->et();
	      }
	    } 
	  }
	  histo("CCJets_MatchedGenPhoton_DeltaR")->Fill( dr1 ); 
	  if ( dr1 > -0.5 ) { histo("CCJets_MatchedGenPhoton_DeltaEt")->Fill( et1 ); }
	  
	  // Nearest (dropped) reco::Photon
	  edm::View<pat::Photon>::const_iterator iii = dropped_photons->begin();
	  edm::View<pat::Photon>::const_iterator jjj = dropped_photons->end();
	  for ( ; iii != jjj; ++iii ) {
	    uint32_t index2 = iii - dropped_photons->begin();
	    edm::RefToBase<pat::Photon> photon_ref = dropped_photons->refAt(index2);
	    if ( photon_ref.isNonnull() ) { 
	      reco::GenParticleRef gen_photon_ref = (*photon_matched_dropped_photons)[photon_ref];
	      if ( gen_jet_ref.isNonnull() && gen_photon_ref.isNonnull() && gen_photon_ref == gen_jet_ref ) { 
		float dr2 = reco::deltaR( iii->p4().eta(), ii->p4().eta(), 
					  iii->p4().phi(), ii->p4().phi() );
		float et2 = fabs( iii->et() - ii->et() ) / ii->et();
		histo("CCJets_MatchedRecoPhoton_DeltaR")->Fill( dr2 ); 
		histo("CCJets_MatchedRecoPhoton_DeltaEt")->Fill( et2 ); 
// 		if      ( iii->isTightPhoton() ) { histo("CCJets_MatchedRecoPhoton_PhotonID")->Fill( 3. ); }
// 		else if ( iii->isLoosePhoton() ) { histo("CCJets_MatchedRecoPhoton_PhotonID")->Fill( 2. ); }
// 		else if ( iii->isLooseEM() )     { histo("CCJets_MatchedRecoPhoton_PhotonID")->Fill( 1. ); }
// 		else                             { histo("CCJets_MatchedRecoPhoton_PhotonID")->Fill( 0. ); }
	      }
	    }
	  }
	  
	}
	
	// Nearest GenElectron
	{
	  float dr = -1.;
	  float et = -1.;
	  uint32_t index = ii - cc_jets->begin();
	  edm::RefToBase<pat::Jet> jet_ref = cc_jets->refAt(index);
	  if ( jet_ref.isNonnull() ) { 
	    reco::GenParticleRef gen_ref = (*jet_matched_cc_electrons)[jet_ref];
	    if ( gen_ref.isNonnull() ) { 
	      reco::GenParticle* gen_particle = const_cast<reco::GenParticle*>( &(*gen_ref) );
	      reco::Candidate* gen_candidate = dynamic_cast<reco::Candidate*>( gen_particle );
	      dr = reco::deltaR( ii->p4().eta(), gen_candidate->p4().eta(),
				 ii->p4().phi(), gen_candidate->p4().phi() );
	      et = fabs( ii->et() - gen_candidate->et() ) / gen_candidate->et();
	    } 
	  } 
	  histo("CCJets_MatchedGenElectron_DeltaR")->Fill( dr ); 
	  if ( dr > -0.5 ) { histo("CCJets_MatchedGenElectron_DeltaEt")->Fill( et ); }
	}
	
	// Nearest GenParton
	{
	  float dr = -1.;
	  float et = -1.;
	  uint32_t index = ii - cc_jets->begin();
	  edm::RefToBase<pat::Jet> jet_ref = cc_jets->refAt(index);
	  if ( jet_ref.isNonnull() ) { 
	    reco::GenParticleRef gen_ref = (*jet_matched_cc_partons)[jet_ref];
	    if ( gen_ref.isNonnull() ) { 
	      reco::GenParticle* gen_particle = const_cast<reco::GenParticle*>( &(*gen_ref) );
	      reco::Candidate* gen_candidate = dynamic_cast<reco::Candidate*>( gen_particle );
	      dr = reco::deltaR( ii->p4().eta(), gen_candidate->p4().eta(),
				 ii->p4().phi(), gen_candidate->p4().phi() );
	      et = fabs( ii->et() - gen_candidate->et() ) / gen_candidate->et();
	    } 
	  } 
	  histo("CCJets_MatchedGenParton_DeltaR")->Fill( dr ); 
	  if ( dr > -0.5 ) { histo("CCJets_MatchedGenParton_DeltaEt")->Fill( et ); }
	}

      }
    }
  }
  

  // -------------------- Dropped jets --------------------


  { 
    histo("DroppedJets_Num")->Fill( dropped_jets->size() ); 
    edm::View<pat::Jet>::const_iterator ii = dropped_jets->begin();
    edm::View<pat::Jet>::const_iterator jj = dropped_jets->end();
    for ( ; ii != jj; ++ii ) {
      histo("DroppedJets_EMfraction")->Fill( ii->emEnergyFraction() ); 
      if ( ii->emEnergyFraction() > em_fraction ) {
	
	// Nearest GenPhoton
	{

	  float dr1 = -1.;
	  float et1 = -1.;
	  uint32_t index1 = ii - dropped_jets->begin();
	  edm::RefToBase<pat::Jet> jet_ref = dropped_jets->refAt(index1);
	  reco::GenParticleRef gen_jet_ref;
	  if ( jet_ref.isNonnull() ) { 
	    gen_jet_ref = (*jet_matched_dropped_photons)[jet_ref];
	    if ( gen_jet_ref.isNonnull() ) { 
	      reco::GenParticle* gen_particle = const_cast<reco::GenParticle*>( &(*gen_jet_ref) );
	      reco::Candidate* gen_candidate = dynamic_cast<reco::Candidate*>( gen_particle );
	      dr1 = reco::deltaR( ii->p4().eta(), gen_candidate->p4().eta(),
				  ii->p4().phi(), gen_candidate->p4().phi() );
	      et1 = fabs( ii->et() - gen_candidate->et() ) / gen_candidate->et();
	    } 
	  }
	  histo("DroppedJets_MatchedGenPhoton_DeltaR")->Fill( dr1 ); 
	  if ( dr1 > -0.5 ) { histo("DroppedJets_MatchedGenPhoton_DeltaEt")->Fill( et1 ); }
	  
	  // Nearest (cc'ed) reco::Photon
	  edm::View<pat::Photon>::const_iterator iii = cc_photons->begin();
	  edm::View<pat::Photon>::const_iterator jjj = cc_photons->end();
	  for ( ; iii != jjj; ++iii ) {
	    uint32_t index2 = iii - cc_photons->begin();
	    edm::RefToBase<pat::Photon> photon_ref = cc_photons->refAt(index2);
	    if ( photon_ref.isNonnull() ) { 
	      reco::GenParticleRef gen_photon_ref = (*photon_matched_cc_photons)[photon_ref];
	      if ( gen_photon_ref == gen_jet_ref ) { 
		float dr2 = reco::deltaR( iii->p4().eta(), ii->p4().eta(),
					  iii->p4().phi(), ii->p4().phi() );
		float et2 = fabs( iii->et() - ii->et() ) / ii->et();
		histo("DroppedJets_MatchedRecoPhoton_DeltaR")->Fill( dr2 ); 
		histo("DroppedJets_MatchedRecoPhoton_DeltaEt")->Fill( et2 ); 
// 		if      ( iii->isTightPhoton() ) { histo("DroppedJets_MatchedRecoPhoton_PhotonID")->Fill( 3. ); }
// 		else if ( iii->isLoosePhoton() ) { histo("DroppedJets_MatchedRecoPhoton_PhotonID")->Fill( 2. ); }
// 		else if ( iii->isLooseEM() )     { histo("DroppedJets_MatchedRecoPhoton_PhotonID")->Fill( 1. ); }
// 		else                             { histo("DroppedJets_MatchedRecoPhoton_PhotonID")->Fill( 0. ); }
	      }
	    }
	  }

	}

	// Nearest GenElectron
	{
	  float dr = -1.;
	  float et = -1.;
	  uint32_t index1 = ii - dropped_jets->begin();
	  edm::RefToBase<pat::Jet> jet_ref = dropped_jets->refAt(index1);
	  if ( jet_ref.isNonnull() ) { 
	    reco::GenParticleRef gen_jet_ref = (*jet_matched_dropped_electrons)[jet_ref];
	    if ( gen_jet_ref.isNonnull() ) { 
	      reco::GenParticle* gen_particle = const_cast<reco::GenParticle*>( &(*gen_jet_ref) );
	      reco::Candidate* gen_candidate = dynamic_cast<reco::Candidate*>( gen_particle );
	      dr = reco::deltaR( ii->p4().eta(), gen_candidate->p4().eta(),
				 ii->p4().phi(), gen_candidate->p4().phi() );
	      et = fabs( ii->et() - gen_candidate->et() ) / gen_candidate->et();
	    } 
	  }
	  histo("DroppedJets_MatchedGenElectron_DeltaR")->Fill( dr ); 
	  if ( dr > -0.5 ) { histo("DroppedJets_MatchedGenElectron_DeltaEt")->Fill( et ); }
	}
	
	// Nearest GenParton
	{
	  float dr = -1.;
	  float et = -1.;
	  uint32_t index1 = ii - dropped_jets->begin();
	  edm::RefToBase<pat::Jet> jet_ref = dropped_jets->refAt(index1);
	  if ( jet_ref.isNonnull() ) { 
	    reco::GenParticleRef gen_jet_ref = (*jet_matched_dropped_partons)[jet_ref];
	    if ( gen_jet_ref.isNonnull() ) { 
	      reco::GenParticle* gen_particle = const_cast<reco::GenParticle*>( &(*gen_jet_ref) );
	      reco::Candidate* gen_candidate = dynamic_cast<reco::Candidate*>( gen_particle );
	      if ( gen_candidate ) {
		dr = reco::deltaR( ii->p4().eta(), gen_candidate->p4().eta(),
				   ii->p4().phi(), gen_candidate->p4().phi() );
		et = fabs( ii->et() - gen_candidate->et() ) / gen_candidate->et();
	      }
	    } 
	  }
	  histo("DroppedJets_MatchedGenParton_DeltaR")->Fill( dr ); 
	  if ( dr > -0.5 ) { histo("DroppedJets_MatchedGenParton_DeltaEt")->Fill( et ); }
	}
	
      }
    }
  }
  



// 	// Nearest cross-cleaned photon
// 	float et = -1.;
// 	float min_dr = 1.e6;
// 	edm::View<pat::Photon>::const_iterator iii = cc_photons->begin();
// 	edm::View<pat::Photon>::const_iterator jjj = cc_photons->end();
// 	for ( ; iii != jjj; ++iii ) {
// 	  if ( fabs( iii->et() - ii->et() ) / ii->et() < 0.5 ) { // allow 50% error
	    
// 	    float dr = reco::deltaR<pat::Photon,pat::Jet>( *ii, *iii );
// 	    if ( dr < min_dr ) { 
// 	      min_dr = dr;
// 	      et = iii->et();
// 	    }
// 	  }
// 	histo("DroppedJets_NearestCCPhoton_DeltaR")->Fill( min ); 
// 	histo("DroppedJets_NearestCCPhoton_Et")->Fill( et ); 
	




// 	  if      ( ii->isTightPhoton() ) { histo("PhotonID_DroppedPhotons")->Fill( 3. ); }
// 	  else if ( ii->isLoosePhoton() ) { histo("PhotonID_DroppedPhotons")->Fill( 2. ); }
// 	  else if ( ii->isLooseEM() )     { histo("PhotonID_DroppedPhotons")->Fill( 1. ); }
// 	  else                            { histo("PhotonID_DroppedPhotons")->Fill( 0. ); }
// 	  histo("Isolation_DroppedPhotons")->Fill( ii->isolation(pat::CaloIso) );
	  
// 	  float et = -1.;
// 	  reco::Candidate* candidate = 0;
// 	  int32_t daughter_id = static_cast<int32_t>(1.e8+0.5);
// 	  std::vector<reco::GenParticleRef> gen_particles = ii->genParticleRefs();
// 	  std::vector<reco::GenParticleRef>::const_iterator igen = gen_particles.begin();
// 	  std::vector<reco::GenParticleRef>::const_iterator jgen = gen_particles.end();
// 	  for ( ; igen != jgen; ++igen ) {
// 	    if ( igen->isNonnull() ) { 
// 	      reco::GenParticle* part = const_cast<reco::GenParticle*>( &(**igen) );
// 	      reco::Candidate* cand = dynamic_cast<reco::Candidate*>( part );
// 	      if ( cand ) {
// 		candidate = cand;
// 		daughter_id = cand->pdgId();
// 		et = cand->et();
// 		break;
// 	      }
// 	    }
// 	  }
// 	  histo("PdgId_DroppedPhotons")->Fill( daughter_id*1. ); 
// 	  if ( et >= 0. ) {
// 	    histo("GenEt_DroppedPhotons")->Fill( et );
// 	    histo("Et_GenEt_DroppedPhotons")->Fill( fabs( et - ii->et() ) );
// 	  }
	  
// 	}
//       }
      

//       histo("DeltaEt_DroppedJets_RawJets")->Fill( fabs( ii->et() - iii->et() ) ); 
// 	  histo("Et_DroppedJets_Vs_RawJets")->Fill( ii->et(), iii->et() ); 
// 	  histo("EMfraction_DroppedJets")->Fill( ii->emEnergyFraction() ); 
	  
// 	  float max = -1.;
// 	  float et = -1.;
// 	  reco::Candidate* candidate = 0;
// 	  int32_t max_id = static_cast<int32_t>(1.e8+0.5);
// 	  int32_t daughter_id = static_cast<int32_t>(1.e8+0.5);
// 	  const reco::GenJet* jet = ii->genJet();
// 	  if ( jet ) {
// 	    histo("GenEt_DroppedJets")->Fill( jet->et() );
// 	    histo("Et_GenEt_DroppedJets")->Fill( fabs( jet->et() - ii->et() ) );
// 	    std::vector<const reco::GenParticle*> gen_particles = jet->getGenConstituents();
// 	    std::vector<const reco::GenParticle*>::const_iterator igen = gen_particles.begin();
// 	    std::vector<const reco::GenParticle*>::const_iterator jgen = gen_particles.end();
// 	    for ( ; igen != jgen; ++igen ) {
// 	      if ( *igen ) { 
// 		reco::GenParticle* part = const_cast<reco::GenParticle*>( *igen );
// 		reco::Candidate* cand = dynamic_cast<reco::Candidate*>( part );
// 		if ( cand ) {
// 		  candidate = cand;
// 		  daughter_id = cand->pdgId();
// 		  if ( daughter_id == 22 ) { et = cand->et(); }
// 		  if ( cand->et() > max ) { max_id = daughter_id; max = cand->et(); }
// 		  break;
// 		}
// 	      }
// 	    }
// 	    if ( max >= 0. ) {
// 	      histo("PdgId_DroppedJets")->Fill( max_id*1. ); 
// 	      histo("MaxEt_DroppedJets")->Fill( max );
// 	    }
// 	  }

// 	}
//       }
//     }
//   }




} 

// -----------------------------------------------------------------------------
//
void TestPatCrossCleaner::beginJob( const edm::EventSetup& ) {

  edm::Service<TFileService> fs;

  { // Misc
    
    TFileDirectory dir = fs->mkdir("Misc");

    histos_["CCJets_Num"] = dir.make<TH1D>("CCJets_Num","",21,-0.5,20.5);
    histos_["CCJets_EMfraction"] = dir.make<TH1D>("CCJets_EMfraction","",120,0.,1.2);
    
    histos_["CCJets_MatchedGenPhoton_DeltaR"] = dir.make<TH1D>("CCJets_MatchedGenPhoton_DeltaR","",100,0.,1.);
    histos_["CCJets_MatchedGenElectron_DeltaR"] = dir.make<TH1D>("CCJets_MatchedGenElectron_DeltaR","",100,0.,1.);
    histos_["CCJets_MatchedGenParton_DeltaR"] = dir.make<TH1D>("CCJets_MatchedGenParton_DeltaR","",100,0.,1.);

    histos_["CCJets_MatchedGenPhoton_DeltaEt"] = dir.make<TH1D>("CCJets_MatchedGenPhoton_DeltaEt","",100,0.,1.);
    histos_["CCJets_MatchedGenElectron_DeltaEt"] = dir.make<TH1D>("CCJets_MatchedGenElectron_DeltaEt","",100,0.,1.);
    histos_["CCJets_MatchedGenParton_DeltaEt"] = dir.make<TH1D>("CCJets_MatchedGenParton_DeltaEt","",100,0.,1.);

    histos_["DroppedJets_Num"] = dir.make<TH1D>("DroppedJets_Num","",21,-0.5,20.5);
    histos_["DroppedJets_EMfraction"] = dir.make<TH1D>("DroppedJets_EMfraction","",120,0.,1.2);

    histos_["DroppedJets_MatchedGenPhoton_DeltaR"] = dir.make<TH1D>("DroppedJets_MatchedGenPhoton_DeltaR","",100,0.,1.);
    histos_["DroppedJets_MatchedGenElectron_DeltaR"] = dir.make<TH1D>("DroppedJets_MatchedGenElectron_DeltaR","",100,0.,1.);
    histos_["DroppedJets_MatchedGenParton_DeltaR"] = dir.make<TH1D>("DroppedJets_MatchedGenParton_DeltaR","",100,0.,1.);

    histos_["DroppedJets_MatchedGenPhoton_DeltaEt"] = dir.make<TH1D>("DroppedJets_MatchedGenPhoton_DeltaEt","",100,0.,1.);
    histos_["DroppedJets_MatchedGenElectron_DeltaEt"] = dir.make<TH1D>("DroppedJets_MatchedGenElectron_DeltaEt","",100,0.,1.);
    histos_["DroppedJets_MatchedGenParton_DeltaEt"] = dir.make<TH1D>("DroppedJets_MatchedGenParton_DeltaEt","",100,0.,1.);

    histos_["CCJets_MatchedRecoPhoton_DeltaR"] = dir.make<TH1D>("CCJets_MatchedRecoPhoton_DeltaR","",100,0.,1.);
    histos_["CCJets_MatchedRecoPhoton_DeltaEt"] = dir.make<TH1D>("CCJets_MatchedRecoPhoton_DeltaEt","",100,0.,1.);
    histos_["CCJets_MatchedRecoPhoton_PhotonID"] = dir.make<TH1D>("CCJets_MatchedRecoPhoton_PhotonID","",4,-0.5,3.5);

    histos_["DroppedJets_MatchedRecoPhoton_DeltaR"] = dir.make<TH1D>("DroppedJets_MatchedRecoPhoton_DeltaR","",100,0.,1.);
    histos_["DroppedJets_MatchedRecoPhoton_DeltaEt"] = dir.make<TH1D>("DroppedJets_MatchedRecoPhoton_DeltaEt","",100,0.,1.);
    histos_["DroppedJets_MatchedRecoPhoton_PhotonID"] = dir.make<TH1D>("DroppedJets_MatchedRecoPhoton_PhotonID","",4,-0.5,3.5);
    
  }
  
}

// -----------------------------------------------------------------------------
// 
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TestPatCrossCleaner);













    
//     histos_["Num_RawPhotons"] = dir.make<TH1D>("Num_RawPhotons","",21,-0.5,20.5);
//     histos_["Num_RawJets"] = dir.make<TH1D>("Num_RawJets","",21,-0.5,20.5);
//     //histos_["Num_RawElectrons"] = dir.make<TH1D>("Num_RawElectrons","",21,-0.5,20.5);
    
//     histos_["Num_CCPhotons"] = dir.make<TH1D>("Num_CCPhotons","",21,-0.5,20.5);
//     histos_["Num_CCJets"] = dir.make<TH1D>("Num_CCJets","",21,-0.5,20.5);
//     //histos_["Num_CCElectrons"] = dir.make<TH1D>("Num_CCElectrons","",21,-0.5,20.5);

//     histos_["Num_DroppedJets"] = dir.make<TH1D>("Num_DroppedJets","",21,-0.5,20.5);
//     histos_["Num_DroppedPhotons"] = dir.make<TH1D>("Num_DroppedPhotons","",21,-0.5,20.5);

    
//     histos2d_["Num_CCPhotons_Vs_RawPhotons"] = dir.make<TH2D>("Num_CCPhotons_Vs_RawPhotons","",21,-0.5,20.5,21,-0.5,20.5);
//     histos2d_["Num_CCJets_Vs_RawJets"] = dir.make<TH2D>("Num_CCJets_Vs_RawJets","",21,-0.5,20.5,21,-0.5,20.5);
//     //histos2d_["Num_CCElectrons_Vs_RawElectrons"] = dir.make<TH2D>("Num_CCElectrons_Vs_RawElectrons","",21,-0.5,20.5,21,-0.5,20.5);
    
//     histos_["PhotonID_CCPhotons"] = dir.make<TH1D>("PhotonID_CCPhotons","",4,-0.5,3.5);
//     histos_["PhotonID_DroppedPhotons"] = dir.make<TH1D>("PhotonID_DroppedPhotons","",4,-0.5,3.5);

//     histos_["Isolation_DroppedPhotons"] = dir.make<TH1D>("Isolation_DroppedPhotons","",100,0.,100.);
//     histos_["PdgId_DroppedPhotons"] = dir.make<TH1D>("PdgId_DroppedPhotons","",201,-100.5,100.5);
//     histos_["GenEt_DroppedPhotons"] = dir.make<TH1D>("GenEt_DroppedPhotons","",1000,0.,1000.);
//     histos_["Et_GenEt_DroppedPhotons"] = dir.make<TH1D>("Et_GenEt_DroppedPhotons","",100,0.,100.);
    
//     histos_["GenEt_DroppedJets"] = dir.make<TH1D>("GenEt_DroppedJets","",1000,0.,1000.);
//     histos_["Et_GenEt_DroppedJets"] = dir.make<TH1D>("Et_GenEt_DroppedJets","",1000,0.,1000.);
//     histos_["PdgId_DroppedJets"] = dir.make<TH1D>("PdgId_DroppedJets","",201,-100.5,100.5);
//     histos_["MaxEt_DroppedJets"] = dir.make<TH1D>("MaxEt_DroppedJets","",1000,0.,1000.);
 
//     histos2d_["Et_CCPhotons_Vs_CCJets"] = dir.make<TH2D>("Et_CCPhotons_Vs_CCJets","",200,0.,1000.,200,0.,1000.);
//     histos2d_["Et_CCJets_Vs_RawJets"] = dir.make<TH2D>("Et_CCJets_Vs_RawJets","",200,0.,1000.,200,0.,1000.);
//     histos2d_["Et_DroppedJets_Vs_RawJets"] = dir.make<TH2D>("Et_DroppedJets_Vs_RawJets","",200,0.,1000.,200,0.,1000.);
    
//     histos_["DeltaEt_CCJets_RawJets"] = dir.make<TH1D>("DeltaEt_CCJets_RawJets","",100,100.,100.);
//     histos_["DeltaEt_DroppedJets_RawJets"] = dir.make<TH1D>("DeltaEt_DroppedJets_RawJets","",100,0.,100.);
    
//     histos_["DeltaR_RawPhotons_RawJets"] = dir.make<TH1D>("DeltaR_RawPhotons_RawJets","",100,0.,5.);
//     //histos_["DeltaR_RawPhotons_RawElectrons"] = dir.make<TH1D>("DeltaR_RawPhotons_RawElectrons","",100,0.,5.);
    
//     histos_["DeltaR_CCPhotons_CCJets"] = dir.make<TH1D>("DeltaR_CCPhotons_CCJets","",100,0.,5.);
//     //histos_["DeltaR_CCPhotons_CCElectrons"] = dir.make<TH1D>("DeltaR_CCPhotons_CCElectrons","",100,0.,5.);
    
//     histos_["EMfraction_DroppedJets"] = dir.make<TH1D>("EMfraction_DroppedJets","",120,0.,1.2);




//   // Photons
//   { 
//     histo("Num_RawPhotons")->Fill( raw_photons->size() ); 
//     histo("Num_CCPhotons")->Fill( cc_photons->size() ); 
//     histo("Num_CCPhotons_Vs_RawPhotons")->Fill( cc_photons->size(), raw_photons->size() ); 
//   }
  
//   // Jets
//   { 
//     histo("Num_RawJets")->Fill( raw_jets->size() ); 
//     histo("Num_CCJets")->Fill( cc_jets->size() ); 
//     histo("Num_CCJets_Vs_RawJets")->Fill( cc_jets->size(), raw_jets->size() ); 
    
//     edm::View<pat::Jet>::const_iterator ii = raw_jets->begin();
//     edm::View<pat::Jet>::const_iterator jj = raw_jets->end();
//     for ( ; ii != jj; ++ii ) {
//       edm::View<pat::Jet>::const_iterator iii = cc_jets->begin();
//       edm::View<pat::Jet>::const_iterator jjj = cc_jets->end();
//       for ( ; iii != jjj; ++iii ) {
// 	if ( ii->originalObjectRef() == iii->originalObjectRef() ) {
// 	  histo("DeltaEt_CCJets_RawJets")->Fill( fabs( ii->et() - iii->et() ) ); 
// 	  histo("Et_CCJets_Vs_RawJets")->Fill( iii->et(), ii->et() ); 
// 	}
//       }
//     }
//   }

//   // Electrons
//   {
//     histo("Num_RawElectrons")->Fill( raw_electrons->size() ); 
//     histo("Num_CCElectrons")->Fill( cc_electrons->size() ); 
//     histo("Num_CCElectrons_Vs_RawElectrons")->Fill( cc_electrons->size(), raw_electrons->size() ); 
//   }


//   // -------------------- Checks b/w dropped objects --------------------


  
//   // Photons
//   { 
//     histo("Num_DroppedPhotons")->Fill( dropped_photons->size() ); 
//     edm::View<pat::Photon>::const_iterator ii = dropped_photons->begin();
//     edm::View<pat::Photon>::const_iterator jj = dropped_photons->end();
//     for ( ; ii != jj; ++ii ) {
//       if      ( ii->isTightPhoton() ) { histo("PhotonID_DroppedPhotons")->Fill( 3. ); }
//       else if ( ii->isLoosePhoton() ) { histo("PhotonID_DroppedPhotons")->Fill( 2. ); }
//       else if ( ii->isLooseEM() )     { histo("PhotonID_DroppedPhotons")->Fill( 1. ); }
//       else                            { histo("PhotonID_DroppedPhotons")->Fill( 0. ); }
//       histo("Isolation_DroppedPhotons")->Fill( ii->isolation(pat::CaloIso) );
      
//       float et = -1.;
//       reco::Candidate* candidate = 0;
//       int32_t daughter_id = static_cast<int32_t>(1.e8+0.5);
//       std::vector<reco::GenParticleRef> gen_particles = ii->genParticleRefs();
//       std::vector<reco::GenParticleRef>::const_iterator igen = gen_particles.begin();
//       std::vector<reco::GenParticleRef>::const_iterator jgen = gen_particles.end();
//       for ( ; igen != jgen; ++igen ) {
// 	if ( igen->isNonnull() ) { 
// 	  reco::GenParticle* part = const_cast<reco::GenParticle*>( &(**igen) );
//  	  reco::Candidate* cand = dynamic_cast<reco::Candidate*>( part );
//  	  if ( cand ) {
// 	    candidate = cand;
// 	    daughter_id = cand->pdgId();
// 	    et = cand->et();
// 	    break;
// 	  }
// 	}
//       }
//       histo("PdgId_DroppedPhotons")->Fill( daughter_id*1. ); 
//       if ( et >= 0. ) {
// 	histo("GenEt_DroppedPhotons")->Fill( et );
// 	histo("Et_GenEt_DroppedPhotons")->Fill( fabs( et - ii->et() ) );
//       }
      
//     }
//   }
  
//   // Jets
//   { 
//     histo("Num_DroppedJets")->Fill( dropped_jets->size() ); 
//     edm::View<pat::Jet>::const_iterator ii = dropped_jets->begin();
//     edm::View<pat::Jet>::const_iterator jj = dropped_jets->end();
//     for ( ; ii != jj; ++ii ) {
//       edm::View<pat::Jet>::const_iterator iii = raw_jets->begin();
//       edm::View<pat::Jet>::const_iterator jjj = raw_jets->end();
//       for ( ; iii != jjj; ++iii ) {
// 	if ( ii->originalObjectRef() == iii->originalObjectRef() ) {

// 	  histo("DeltaEt_DroppedJets_RawJets")->Fill( fabs( ii->et() - iii->et() ) ); 
// 	  histo("Et_DroppedJets_Vs_RawJets")->Fill( ii->et(), iii->et() ); 
// 	  histo("EMfraction_DroppedJets")->Fill( ii->emEnergyFraction() ); 
	  
// 	  float max = -1.;
// 	  float et = -1.;
// 	  reco::Candidate* candidate = 0;
// 	  int32_t max_id = static_cast<int32_t>(1.e8+0.5);
// 	  int32_t daughter_id = static_cast<int32_t>(1.e8+0.5);
// 	  const reco::GenJet* jet = ii->genJet();
// 	  if ( jet ) {
// 	    histo("GenEt_DroppedJets")->Fill( jet->et() );
// 	    histo("Et_GenEt_DroppedJets")->Fill( fabs( jet->et() - ii->et() ) );
// 	    std::vector<const reco::GenParticle*> gen_particles = jet->getGenConstituents();
// 	    std::vector<const reco::GenParticle*>::const_iterator igen = gen_particles.begin();
// 	    std::vector<const reco::GenParticle*>::const_iterator jgen = gen_particles.end();
// 	    for ( ; igen != jgen; ++igen ) {
// 	      if ( *igen ) { 
// 		reco::GenParticle* part = const_cast<reco::GenParticle*>( *igen );
// 		reco::Candidate* cand = dynamic_cast<reco::Candidate*>( part );
// 		if ( cand ) {
// 		  candidate = cand;
// 		  daughter_id = cand->pdgId();
// 		  if ( daughter_id == 22 ) { et = cand->et(); }
// 		  if ( cand->et() > max ) { max_id = daughter_id; max = cand->et(); }
// 		  break;
// 		}
// 	      }
// 	    }
// 	    if ( max >= 0. ) {
// 	      histo("PdgId_DroppedJets")->Fill( max_id*1. ); 
// 	      histo("MaxEt_DroppedJets")->Fill( max );
// 	    }
// 	  }

// 	}
//       }
//     }
//   }


//   // -------------------- Checks b/w UN-cross-cleaned objects --------------------

  
//   // UNcross-cleaned photons-jets
//   {
//     edm::View<pat::Photon>::const_iterator ii = raw_photons->begin();
//     edm::View<pat::Photon>::const_iterator jj = raw_photons->end();
//     for ( ; ii != jj; ++ii ) {
//       float min = 1.e6;
//       edm::View<pat::Jet>::const_iterator iii = raw_jets->begin();
//       edm::View<pat::Jet>::const_iterator jjj = raw_jets->end();
//       for ( ; iii != jjj; ++iii ) {
// 	float delta = reco::deltaR<pat::Photon,pat::Jet>( *ii, *iii );
// 	if ( delta < min ) { min = delta; }
//       }
//       histo("DeltaR_RawPhotons_RawJets")->Fill( min ); 
//     }
//   }
  
// //   // UNcross-cleaned photons-electrons
// //   {
// //     edm::View<pat::Photon>::const_iterator ii = raw_photons->begin();
// //     edm::View<pat::Photon>::const_iterator jj = raw_photons->end();
// //     for ( ; ii != jj; ++ii ) {
// //       float min = 1.e6;
// //       edm::View<pat::Electron>::const_iterator iii = raw_electrons->begin();
// //       edm::View<pat::Electron>::const_iterator jjj = raw_electrons->end();
// //       for ( ; iii != jjj; ++iii ) {
// // 	float delta = reco::deltaR<pat::Photon,pat::Electron>( *ii, *iii );
// // 	if ( delta < min ) { min = delta; }
// //       }
// //       histo("DeltaR_RawPhotons_RawElectrons")->Fill( min ); 
// //     }
// //   }
  
  
//   // -------------------- Checks b/w cross-cleaned objects --------------------


//   // Cross-cleaned photons-jets
//   {
//     edm::View<pat::Photon>::const_iterator ii = cc_photons->begin();
//     edm::View<pat::Photon>::const_iterator jj = cc_photons->end();
//     for ( ; ii != jj; ++ii ) {
//       float et = -1.;
//       float min = 1.e6;
//       edm::View<pat::Jet>::const_iterator iii = cc_jets->begin();
//       edm::View<pat::Jet>::const_iterator jjj = cc_jets->end();
//       for ( ; iii != jjj; ++iii ) {
// 	float delta = reco::deltaR<pat::Photon,pat::Jet>( *ii, *iii );
// 	if ( delta < min ) { 
// 	  min = delta;
// 	  et = iii->et();
// 	}
//       }
//       histo("DeltaR_CCPhotons_CCJets")->Fill( min ); 
//       if ( min < 0.4 ) { 
// 	if      ( ii->isTightPhoton() ) { histo("PhotonID_CCPhotons")->Fill( 3. ); }
// 	else if ( ii->isLoosePhoton() ) { histo("PhotonID_CCPhotons")->Fill( 2. ); }
// 	else if ( ii->isLooseEM() )     { histo("PhotonID_CCPhotons")->Fill( 1. ); }
// 	else                            { histo("PhotonID_CCPhotons")->Fill( 0. ); }
// 	histo("Et_CCPhotons_Vs_CCJets")->Fill( ii->et(), et ); 
//       }
//     }
//   }

// //   // Cross-cleaned photons-electrons
// //   {
// //     edm::View<pat::Photon>::const_iterator ii = cc_photons->begin();
// //     edm::View<pat::Photon>::const_iterator jj = cc_photons->end();
// //     for ( ; ii != jj; ++ii ) {
// //       float min = 1.e6;
// //       edm::View<pat::Electron>::const_iterator iii = cc_electrons->begin();
// //       edm::View<pat::Electron>::const_iterator jjj = cc_electrons->end();
// //       for ( ; iii != jjj; ++iii ) {
// // 	float delta = reco::deltaR<pat::Photon,pat::Electron>( *ii, *iii );
// // 	if ( delta < min ) { min = delta; }
// //       }
// //       histo("DeltaR_CCPhotons_CCElectrons")->Fill( min ); 
// //     }
// //   }




  
