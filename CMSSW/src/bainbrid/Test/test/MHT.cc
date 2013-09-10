#include "MHT.h"
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
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cmath>

// -----------------------------------------------------------------------------
//
MHT::MHT( const edm::ParameterSet& pset ) : 
  // Primary Objects
  jets_( pset.getParameter<edm::InputTag>("Jets") ),
  muons_( pset.getParameter<edm::InputTag>("Muons") ),
  electrons_( pset.getParameter<edm::InputTag>("Electrons") ),
  photons_( pset.getParameter<edm::InputTag>("Photons") ),
  // MET 
  caloMet_( pset.getParameter<edm::InputTag>("CaloMET") ),
  genMet_( pset.getParameter<edm::InputTag>("GenMET") ),
  // Jet kinematics
  jetEt_( pset.getParameter<double>("JetEt") ),
  jetEta_( pset.getParameter<double>("JetEta") ),
  jetEMfrac_( pset.getParameter<double>("JetEMfraction") ),
  // Muon kinematics
  muonPt_( pset.getParameter<double>("MuonPt") ),
  muonEta_( pset.getParameter<double>("MuonEta") ),
  muonTrkIso_( pset.getParameter<double>("MuonTrkIso") ),
  // Electron kinematics
  electronPt_( pset.getParameter<double>("ElectronPt") ),
  electronEta_( pset.getParameter<double>("ElectronEta") ),
  electronTrkIso_( pset.getParameter<double>("ElectronTrkIso") ),
  // Photon kinematics
  photonEt_( pset.getParameter<double>("PhotonEt") ),
  photonEta_( pset.getParameter<double>("PhotonEta") ),
  // Event selection
  totalHt_( pset.getParameter<double>("TotalEt") ),
  minObjects_( pset.getParameter<int>("MinObjects") ),
  minJets_( pset.getParameter<int>("MinJets") ),
  minMuons_( pset.getParameter<int>("MinMuons") ),
  minElectrons_( pset.getParameter<int>("MinElectrons") ),
  minPhotons_( pset.getParameter<int>("MinPhotons") )
{
  produces<Candidates>("PrimaryObjects");
  produces<Candidates>("Jets");
  produces<Candidates>("Muons");
  produces<Candidates>("Electrons");
  produces<Candidates>("Photons");
  produces<Candidates>("MHTs");
  produces<Candidates>("CaloMETs");
  produces<Candidates>("GenMETs");
}

// -----------------------------------------------------------------------------
//
void MHT::produce( edm::Event& event, 
		   const edm::EventSetup& setup ) {

  // -------------------- Retrieve objects --------------------

  // Get photons
  std::auto_ptr<Candidates> photons( new Candidates );
  if ( getPhotons( event, *photons ) ) { return; }
  
  // Get jets
  std::auto_ptr<Candidates> jets( new Candidates );
  std::auto_ptr<Candidates> gen_jets( new Candidates );
  if ( getJets( event, *jets, *gen_jets ) ) { return; }
    
  // Get muons
  std::auto_ptr<Candidates> muons( new Candidates );
  if ( getMuons( event, *muons ) ) { return; }

  // Get electrons
  std::auto_ptr<Candidates> electrons( new Candidates );
  if ( getElectrons( event, *electrons ) ) { return; }
  
  // -------------------- Common objects --------------------
    
  std::auto_ptr<Candidates> objects( new Candidates );

  if ( !jets->empty() ) { objects->insert( objects->end(), jets->begin(), jets->end() ); }
  if ( !muons->empty() ) { objects->insert( objects->end(), muons->begin(), muons->end() ); }
  if ( !electrons->empty() ) { objects->insert( objects->end(), electrons->begin(), electrons->end() ); }
  if ( !photons->empty() ) { objects->insert( objects->end(), photons->begin(), photons->end() ); }

  std::sort( objects->begin(), objects->end(), GreaterByEt<Candidate>() );

  // -------------------- Event selection --------------------
  
  if ( objects->size() < minObjects_ ) { return; }
  if ( jets->size() < minJets_ ) { return; }
  if ( muons->size() < minMuons_ ) { return; }
  if ( electrons->size() < minElectrons_ ) { return; }
  if ( photons->size() < minPhotons_ ) { return; }
  
  if ( ht( *objects ) < totalHt_ ) { return; }

  // -------------------- METs --------------------
  
  std::auto_ptr<Candidates> mhts( new Candidates );
  std::auto_ptr<Candidates> calo_mets( new Candidates );
  std::auto_ptr<Candidates> gen_mets( new Candidates );

  // MHT
  double sum_et = 0.;
  double sum_px = 0.;
  double sum_py = 0.;
  Candidates::const_iterator iobj = objects->begin();
  Candidates::const_iterator jobj = objects->end();
  for ( ; iobj != jobj; ++iobj ) { 
    sum_et +=iobj->p4().Et();
    sum_px +=iobj->p4().Px();
    sum_py +=iobj->p4().Py();
  }
  LorentzVector lv_mht( -1.*sum_px,
			-1.*sum_py,
			0.,
			sqrt( sum_px*sum_px + sum_py*sum_py ) );
  mhts->push_back( reco::MET( sum_et, lv_mht, math::XYZPoint() ) );
  
  // CaloMET  
  if ( !caloMet_.label().empty() ) {
    edm::Handle< edm::View<pat::MET> > handle;
    event.getByLabel(caloMet_,handle);
    if ( handle.isValid() && 
	 !handle->empty() && 
	 handle->front().isCaloMET() ) { 
      calo_mets->push_back( handle->front() );
    } else {
      edm::LogError("TEST") << "Unable to retrieve CaloMET collection!";
    }
  }
  
  // GenMET  
  if ( !genMet_.label().empty() ) {
    edm::Handle< edm::View<reco::GenMET> > handle;
    event.getByLabel(genMet_,handle);
    if ( handle.isValid() && 
	 !handle->empty() ) { 
      gen_mets->push_back( handle->front() );
    } else {
      edm::LogError("TEST") << "Unable to retrieve GenMET collection!";
    }
  }
  
  // -------------------- Fill event -------------------- 

  if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Objects: " << objects->size(); }

  event.put(objects,"PrimaryObjects");
  event.put(jets,"Jets");
  event.put(muons,"Muons");
  event.put(electrons,"Electrons");
  event.put(photons,"Photons");
  event.put(mhts,"MHTs");
  event.put(calo_mets,"CaloMETs");
  event.put(gen_mets,"GenMETs");

}    

// -----------------------------------------------------------------------------
//
void MHT::beginJob( const edm::EventSetup& ) {}

// -----------------------------------------------------------------------------
//
bool MHT::getJets( const edm::Event& event,
		   Candidates& jets,
		   Candidates& gen_jets ) {
  
  jets.clear();
  gen_jets.clear();
  
  if ( !jets_.label().empty() ) {

    edm::Handle< std::vector<pat::Jet> > handle;
    event.getByLabel(jets_,handle);
    
    if ( !handle.isValid() ) { 
      edm::LogError("TEST") << "No Jets for " << jets_; 
      return true;
    }
    
    std::vector<pat::Jet>::const_iterator ijet = handle->begin();
    std::vector<pat::Jet>::const_iterator jjet = handle->end();
    for ( ; ijet != jjet; ++ijet  ) { 
      if ( fabs( ijet->eta() ) < jetEta_ && 
	   ijet->emEnergyFraction() < jetEMfrac_ &&
	    ijet->et() > jetEt_ ) { jets.push_back( *ijet ); }
    }
    
    if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Jets: " << jets.size(); }
    
  }

//   std::sort( jets.begin(), jets.end(), GreaterByEt<Candidate>() );
//   Candidates::const_iterator ii = jets.begin();
//   Candidates::const_iterator jj = jets.end();
//   for ( ; ii != jj; ++ii  ) { 
//     if ( ii->genJet() ) { gen_jets.push_back( *(ii->genJet()) ); }
//     else { gen_jets.push_back( Candidate() ) }
//   }

  return false;
  
}

// -----------------------------------------------------------------------------
//
bool MHT::getMuons( const edm::Event& event,
		    Candidates& muons ) {
  
  muons.clear();
  
  if ( !muons_.label().empty() ) {

    edm::Handle< std::vector<pat::Muon> > handle;
    event.getByLabel(muons_,handle);
    
    if ( !handle.isValid() ) { 
      edm::LogError("TEST") << "No Muons for " << muons_; 
      return true;
    }
    
    std::vector<pat::Muon>::const_iterator imuon = handle->begin();
    std::vector<pat::Muon>::const_iterator jmuon = handle->end();
    for ( ; imuon != jmuon; ++imuon  ) { 
      if ( fabs( imuon->eta() ) < muonEta_ &&
	   imuon->trackIso() < muonTrkIso_ &&
	   imuon->pt() > muonPt_ ) { muons.push_back( *imuon ); }
    }
    
    if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Muons: " << muons.size(); }

  }

  std::sort( muons.begin(), muons.end(), GreaterByEt<Candidate>() );

  return false;
  
}

// -----------------------------------------------------------------------------
//
bool MHT::getElectrons( const edm::Event& event,
			Candidates& electrons ) {
  
  electrons.clear();
  
  if ( !electrons_.label().empty() ) {

    edm::Handle< std::vector<pat::Electron> > handle;
    event.getByLabel(electrons_,handle);

    if ( !handle.isValid() ) { 
      edm::LogError("TEST") << "No Electrons for " << electrons_; 
      return true;
    }
    
    std::vector<pat::Electron>::const_iterator ielectron = handle->begin();
    std::vector<pat::Electron>::const_iterator jelectron = handle->end();
    for ( ; ielectron != jelectron; ++ielectron  ) { 
      if ( fabs( ielectron->eta() ) < electronEta_ &&
	   ielectron->trackIso() < electronTrkIso_ &&
	   ielectron->pt() > electronPt_ ) { electrons.push_back( *ielectron ); }
    }
    
    if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Electrons: " << electrons.size(); }
    
  }

  std::sort( electrons.begin(), electrons.end(), GreaterByEt<Candidate>() );
  
  return false;
  
}

// -----------------------------------------------------------------------------
//
bool MHT::getPhotons( const edm::Event& event,
		      Candidates& photons ) {
  
  photons.clear();
  
  if ( !photons_.label().empty() ) {
    
    edm::Handle< std::vector<pat::Photon> > handle;
    event.getByLabel(photons_,handle);
    
    if ( !handle.isValid() ) { 
      edm::LogError("TEST") << "No Photons for " << photons_; 
      return true;
    }
    
    std::vector<pat::Photon>::const_iterator iphoton = handle->begin();
    std::vector<pat::Photon>::const_iterator jphoton = handle->end();
    for ( ; iphoton != jphoton; ++iphoton  ) { 
      if ( fabs( iphoton->eta() ) < photonEta_ &&
	   iphoton->et() > photonEt_ ) { photons.push_back( *iphoton ); }
    }
    
    if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Photons: " << photons.size(); }

  }

  std::sort( photons.begin(), photons.end(), GreaterByEt<Candidate>() );
  
  return false; 

}

// -----------------------------------------------------------------------------
//
double MHT::ht( const Candidates& input ) {
  LorentzVector lv;
  Candidates::const_iterator ii = input.begin();
  Candidates::const_iterator jj = input.end();
  for ( ; ii != jj; ++ii ) { lv += ii->p4(); }
  double ht = lv.E()*lv.E() - lv.Pz()*lv.Pz();
  ht = ht < 0. ? -1.*sqrt(-1.*ht) : sqrt(ht);
  return ht;
}

// -----------------------------------------------------------------------------
//
double MHT::mht( const Candidates& input ) {
  LorentzVector lv_obj;
  Candidates::const_iterator iobj = input.begin();
  Candidates::const_iterator jobj = input.end();
  for ( ; iobj != jobj; ++iobj ) { lv_obj += iobj->p4(); }
  LorentzVector lv_mht( lv_obj );
  lv_mht.SetPx( -1.*lv_obj.Px() );
  lv_mht.SetPy( -1.*lv_obj.Py() );
  lv_mht.SetPz( -1.*lv_obj.Pz() );
  return sqrt( lv_mht.Perp2() );
}

// -----------------------------------------------------------------------------
// 
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MHT);
