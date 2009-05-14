#include "TestCombination.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/Utilities/interface/EtComparator.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "Combination.h"
#include "TH1D.h"
#include "TH2D.h"
#include <algorithm>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cmath>

// -----------------------------------------------------------------------------
//
TestCombination::TestCombination( const edm::ParameterSet& pset ) 
  : maximum_( pset.getParameter<int>("MaximumObjects") ),
    test_( pset.getParameter<int>("TestObjects") ),
    photons_( pset.getParameter<edm::InputTag>("Photons") ),
    jets_( pset.getParameter<edm::InputTag>("Jets") ),
    muons_( pset.getParameter<edm::InputTag>("Muons") ),
    electrons_( pset.getParameter<edm::InputTag>("Electrons") ),
    photonEt_( pset.getParameter<double>("PhotonEt") ),
    photonEta_( pset.getParameter<double>("PhotonEta") ),
    jetEt_( pset.getParameter<double>("JetEt") ),
    jetEta_( pset.getParameter<double>("JetEta") ),
    jetEMfrac_( pset.getParameter<double>("JetEMfraction") ),
    muonPt_( pset.getParameter<double>("MuonPt") ),
    muonEta_( pset.getParameter<double>("MuonEta") ),
    muonTrkIso_( pset.getParameter<double>("MuonTrkIso") ),
    electronPt_( pset.getParameter<double>("ElectronPt") ),
    electronEta_( pset.getParameter<double>("ElectronEta") ),
    electronTrkIso_( pset.getParameter<double>("ElectronTrkIso") ),
    totalEt_( pset.getParameter<double>("TotalEt") ),
    histos_(),
    histos2d_()
{;}

// -----------------------------------------------------------------------------
//
void TestCombination::analyze( const edm::Event& iEvent, 
			       const edm::EventSetup& iSetup ) {

  if ( test_ < 0 ) {

    // Get photons
    std::vector<reco::Particle> photons;
    if ( getPhotons( iEvent, photons ) ) { return; }

    // Get jets
    std::vector<reco::Particle> jets;
    if ( getJets( iEvent, jets ) ) { return; }

    // Get muons
    std::vector<reco::Particle> muons;
    if ( getMuons( iEvent, muons ) ) { return; }

    // Get electrons
    std::vector<reco::Particle> electrons;
    if ( getElectrons( iEvent, electrons ) ) { return; }
  
    // Event selection based on "common objects"
    if ( photons.size() < 2 ) { return; }
    if ( !muons.empty() || !electrons.empty() ) { return; }
    
    // Build vector of photon and jet objects
    std::vector<reco::Particle> objects;
    if ( !photons.empty() ) { objects.insert( objects.end(), photons.begin(), photons.end() ); }
    if ( !jets.empty() ) { objects.insert( objects.end(), jets.begin(), jets.end() ); }
    if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Objects: " << objects.size(); }
    
    // Order objects by Et
    sort( objects.begin(), objects.end(), GreaterByEt<reco::Particle>() );

    // Calculate sum of Et
    float et = 0.;
    std::vector<reco::Particle>::const_iterator iobject = objects.begin();
    std::vector<reco::Particle>::const_iterator jobject = objects.end();
    for ( ; iobject != jobject; ++iobject  ) { et += iobject->et(); }

    // Calculate alphaT
    std::vector<reco::Particle> jet1;
    std::vector<reco::Particle> jet2;
    float min_et = minDeltaEt( objects, jet1, jet2 );
    float sum_et = sumEt( objects );
    float mass_t = massT( objects );
    float alpha_t = alphaT( min_et, sum_et, mass_t );
    
    if ( et > 250. && et < 350. ) {
      histos_["AlphaT_250"]->Fill( alpha_t ); 
      histos2d_["AlphaT_Vs_NObjects_250"]->Fill( objects.size(), alpha_t ); 
    }

    // Cut on Et sum of "common objects"
    if ( et < totalEt_ ) { return; }
    
    histos_["TotalEt"]->Fill( et ); 
    histos2d_["TotalEt_Vs_NObjects"]->Fill( objects.size(), et ); 
    
    histos_["PostCuts_NObjects"]->Fill( objects.size() ); 
    histos_["PostCuts_NPhotons"]->Fill( photons.size() ); 
    histos_["PostCuts_NJets"]->Fill( jets.size() ); 

    histos_["SumEt"]->Fill( sum_et ); 
    histos2d_["SumEt_Vs_NObjects"]->Fill( objects.size(), sum_et ); 

    histos_["MassT"]->Fill( mass_t ); 
    histos2d_["MassT_Vs_NObjects"]->Fill( objects.size(), mass_t ); 
    
    histos_["MinDeltaEt"]->Fill( min_et ); 
    histos2d_["MinDeltaEt_Vs_NObjects"]->Fill( objects.size(), min_et ); 
    
    histos_["AlphaT"]->Fill( alpha_t ); 
    histos2d_["AlphaT_Vs_NObjects"]->Fill( objects.size(), alpha_t ); 
    
    std::vector<int> phot1;
    std::vector<reco::Particle>::const_iterator ii1 = jet1.begin();
    std::vector<reco::Particle>::const_iterator jj1 = jet1.end();
    for ( ; ii1 != jj1; ++ii1 ) { if ( ii1->pdgId() == 22 ) { phot1.push_back( ii1 - jet1.begin() ); } }
    
    std::vector<int> phot2;
    std::vector<reco::Particle>::const_iterator ii2 = jet2.begin();
    std::vector<reco::Particle>::const_iterator jj2 = jet2.end();
    for ( ; ii2 != jj2; ++ii2 ) { if ( ii2->pdgId() == 22 ) { phot2.push_back( ii2 - jet2.begin() ); } }
  
    histos2d_["PseudoDijets_NObjects"]->Fill( jet1.size(), jet2.size() ); 
    histos2d_["PseudoDijets_NPhotons"]->Fill( phot1.size(), phot2.size() ); 
    histos2d_["PseudoDijets_NJets"]->Fill( jet1.size() - phot1.size(), jet2.size() - phot2.size() ); 
    
    {
      std::vector<int>::const_iterator iii1 = phot1.begin();
      std::vector<int>::const_iterator jjj1 = phot1.end();
      for ( ; iii1 != jjj1; ++iii1 ) {
	std::vector<int>::const_iterator iii2 = phot2.begin();
	std::vector<int>::const_iterator jjj2 = phot2.end();
	for ( ; iii2 != jjj2; ++iii2 ) {
	  histos2d_["PseudoDijets_PhotonEt"]->Fill( jet1[*iii1].et(), jet2[*iii2].et() ); 
	  float photon_et = jet1[*iii1].et() < jet2[*iii2].et() ? jet1[*iii1].et() : jet2[*iii2].et();
	  float delta_phi = reco::deltaPhi<reco::Particle,reco::Particle>( jet1[*iii1], jet2[*iii2] );
	  histos2d_["PseudoDijets_DeltaPhiVsPhotonEt"]->Fill( delta_phi, photon_et ); 
	}
      }
    }

    // Debug
    if ( edm::isDebugEnabled() ) {
      std::stringstream ss;
      ss << "Minimum deltaEt = " << min_et << std::endl;

      {
	ss << "Input objects: id/pdg: ";
	std::vector<reco::Particle>::const_iterator ii = objects.begin();
	std::vector<reco::Particle>::const_iterator jj = objects.end();
	for ( ; ii != jj; ++ii  ) { 
	  ss << ii->pdgId() << "/" << ii->pt() << ", "; 
	}
	ss << std::endl;
      }

      {
	ss << "Pseudo-jet #1: id/pdg: ";
	std::vector<reco::Particle>::const_iterator ii = jet1.begin();
	std::vector<reco::Particle>::const_iterator jj = jet1.end();
	for ( ; ii != jj; ++ii  ) { 
	  ss << ii->pdgId() << "/" << ii->pt() << ", "; 
	}
	ss << std::endl;
      }

      {
	ss << "Pseudo-jet #2: id/pdg: ";
	std::vector<reco::Particle>::const_iterator ii = jet2.begin();
	std::vector<reco::Particle>::const_iterator jj = jet2.end();
	for ( ; ii != jj; ++ii  ) { 
	  ss << ii->pdgId() << "/" << ii->pt() << ", "; 
	}
	ss << std::endl;
      }
    
      LogTrace("TEST") << ss.str();

    }

  } else {

    // Test
    std::vector<reco::Particle> input;
    for ( uint16_t ii = 1; ii <= (iEvent.id().event()-1)%test_; ++ii ) {
      input.push_back( reco::Particle( 0, reco::Particle::LorentzVector( 1./ii, 0, 0, ii ) ) );
    }
    std::vector<reco::Particle> jet1;
    std::vector<reco::Particle> jet2;
    float min_et = minDeltaEt( input, jet1, jet2 );
    
    // Debug
    if ( edm::isDebugEnabled() ) {
      std::stringstream ss;
      ss << "Minimum deltaEt = " << min_et << std::endl;

      {
	ss << "Input objects: pdg/pt: ";
	std::vector<reco::Particle>::const_iterator ii = input.begin();
	std::vector<reco::Particle>::const_iterator jj = input.end();
	for ( ; ii != jj; ++ii  ) { 
	  ss << ii->pdgId() << "/" << ii->pt() << ", "; 
	}
	ss << std::endl;
      }

      {
	ss << "Pseudo-jet #1: pdg/pt: ";
	std::vector<reco::Particle>::const_iterator ii = jet1.begin();
	std::vector<reco::Particle>::const_iterator jj = jet1.end();
	for ( ; ii != jj; ++ii  ) { 
	  ss << ii->pdgId() << "/" << ii->pt() << ", "; 
	}
	ss << std::endl;
      }

      {
	ss << "Pseudo-jet #2: pdg/pt: ";
	std::vector<reco::Particle>::const_iterator ii = jet2.begin();
	std::vector<reco::Particle>::const_iterator jj = jet2.end();
	for ( ; ii != jj; ++ii  ) { 
	  ss << ii->pdgId() << "/" << ii->pt() << ", "; 
	}
	ss << std::endl;
      }
    
      LogTrace("TEST") << ss.str();

    }

  }    
  
}

// -----------------------------------------------------------------------------
//
void TestCombination::beginJob( const edm::EventSetup& ) {

  if ( test_ < 0 ) {

    edm::Service<TFileService> fs;
  
    { // Pre-cuts
      
      TFileDirectory dir = fs->mkdir("PreCuts");
    
      histos_["PreCuts_NPhotons"] = dir.make<TH1D>("PreCuts_NPhotons","",51,-0.5,50.5);
      histos_["PreCuts_NPhotons25"] = dir.make<TH1D>("PreCuts_NPhotons25","",51,-0.5,50.5);
      histos_["PreCuts_NPhotons30"] = dir.make<TH1D>("PreCuts_NPhotons30","",51,-0.5,50.5);
      histos_["PreCuts_NPhotons35"] = dir.make<TH1D>("PreCuts_NPhotons35","",51,-0.5,50.5);
      histos_["PreCuts_NPhotons40"] = dir.make<TH1D>("PreCuts_NPhotons40","",51,-0.5,50.5);
      histos_["PreCuts_NPhotons45"] = dir.make<TH1D>("PreCuts_NPhotons45","",51,-0.5,50.5);
      histos_["PreCuts_NPhotons50"] = dir.make<TH1D>("PreCuts_NPhotons50","",51,-0.5,50.5);
    
      histos_["PreCuts_NJets"] = dir.make<TH1D>("PreCuts_NJets","",51,-0.5,50.5);
      histos_["PreCuts_NJets25"] = dir.make<TH1D>("PreCuts_NJets25","",51,-0.5,50.5);
      histos_["PreCuts_NJets30"] = dir.make<TH1D>("PreCuts_NJets30","",51,-0.5,50.5);
      histos_["PreCuts_NJets35"] = dir.make<TH1D>("PreCuts_NJets35","",51,-0.5,50.5);
      histos_["PreCuts_NJets40"] = dir.make<TH1D>("PreCuts_NJets40","",51,-0.5,50.5);
      histos_["PreCuts_NJets45"] = dir.make<TH1D>("PreCuts_NJets45","",51,-0.5,50.5);
      histos_["PreCuts_NJets50"] = dir.make<TH1D>("PreCuts_NJets50","",51,-0.5,50.5);

      histos_["PreCuts_NMuons"] = dir.make<TH1D>("PreCuts_NMuons","",51,-0.5,50.5);
      histos_["PreCuts_NElectrons"] = dir.make<TH1D>("PreCuts_NElectrons","",51,-0.5,50.5);
    }

    { // Post-cuts

      TFileDirectory dir = fs->mkdir("PostCuts");
      
      histos_["PostCuts_NObjects"] = dir.make<TH1D>("PostCuts_NObjects","",51,-0.5,50.5);
      histos_["PostCuts_NPhotons"] = dir.make<TH1D>("PostCuts_NPhotons","",51,-0.5,50.5);
      histos_["PostCuts_NJets"] = dir.make<TH1D>("PostCuts_NJets","",51,-0.5,50.5);
      
    }

    { // Alpha_T
      
      TFileDirectory dir = fs->mkdir("AlphaT");
      
      histos_["TotalEt"] = dir.make<TH1D>("TotalEt","",100,0.,2000.);
      histos2d_["TotalEt_Vs_NObjects"] = dir.make<TH2D>("TotalEt_Vs_NObjects","",51,-0.5,50.5,100,0.,2000. );
      
      histos_["SumEt"] = dir.make<TH1D>("SumEt","",100,0.,2000.);
      histos2d_["SumEt_Vs_NObjects"] = dir.make<TH2D>("SumEt_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);
      
      histos_["MassT"] = dir.make<TH1D>("MassT","",100,0.,2000.);
      histos2d_["MassT_Vs_NObjects"] = dir.make<TH2D>("MassT_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);
      
      histos_["MinDeltaEt"] = dir.make<TH1D>("MinDeltaEt","",100,0.,100.);
      histos2d_["MinDeltaEt_Vs_NObjects"] = dir.make<TH2D>("MinDeltaEt_Vs_NObjects","",51,-0.5,50.5,100,0.,100.);
      
      histos_["AlphaT"] = dir.make<TH1D>("AlphaT","",300,-0.5,2.5);
      histos2d_["AlphaT_Vs_NObjects"] = dir.make<TH2D>("AlphaT_Vs_NObjects","",51,-0.5,50.5,300,-0.5,2.5);

      histos_["AlphaT_250"] = dir.make<TH1D>("AlphaT_250","",300,-0.5,2.5);
      histos2d_["AlphaT_Vs_NObjects_250"] = dir.make<TH2D>("AlphaT_Vs_NObjects_250","",51,-0.5,50.5,300,-0.5,2.5);
      
    }
    
    { // Pseudo-dijets
      
      TFileDirectory dir = fs->mkdir("PseudoDijets");
      
      histos2d_["PseudoDijets_NObjects"] = dir.make<TH2D>("PseudoDijets_NObjects","",51,-0.5,50.5,51,-0.5,50.5); 
      histos2d_["PseudoDijets_NPhotons"] = dir.make<TH2D>("PseudoDijets_NPhotons","",51,-0.5,50.5,51,-0.5,50.5); 
      histos2d_["PseudoDijets_NJets"] = dir.make<TH2D>("PseudoDijets_NJets","",51,-0.5,50.5,51,-0.5,50.5); 
      
      histos2d_["PseudoDijets_PhotonEt"] = dir.make<TH2D>("PseudoDijets_PhotonEt","",100,0.,1000.,100,0.,1000.);
      histos2d_["PseudoDijets_DeltaPhiVsPhotonEt"] = dir.make<TH2D>("PseudoDijets_DeltaPhiVsPhotonEt","",100,0.,5.,100,0.,1000.);
      
    }

  } else { 
    edm::LogWarning("TEST") << "TEST! Number of objects as input: " << test_; 
  }
  
}

// -----------------------------------------------------------------------------
//
bool TestCombination::getPhotons( const edm::Event& iEvent,
				  std::vector<reco::Particle>& photons ) {

  photons.clear();
  
  if ( !photons_.label().empty() ) {
    
    edm::Handle< std::vector<pat::Photon> > handle;
    iEvent.getByLabel(photons_,handle);
    
    if ( !handle.isValid() ) { 
      edm::LogWarning("TEST") << "No Photons for " << photons_; 
      return true;
    }
    
    uint16_t et25=0, et30=0, et35=0, et40=0, et45=0, et50=0;
    std::vector<pat::Photon>::const_iterator iphoton = handle->begin();
    std::vector<pat::Photon>::const_iterator jphoton = handle->end();
    for ( ; iphoton != jphoton; ++iphoton  ) { 
      if ( fabs( iphoton->eta() ) < photonEta_ ) { 
	if ( iphoton->et() > photonEt_ ) { photons.push_back( *iphoton ); }
	if ( iphoton->et() > 25. ) { et25++; }
	if ( iphoton->et() > 30. ) { et30++; }
	if ( iphoton->et() > 35. ) { et35++; }
	if ( iphoton->et() > 40. ) { et40++; }
	if ( iphoton->et() > 45. ) { et45++; }
	if ( iphoton->et() > 50. ) { et50++; }
      }
    }
    
    histos_["PreCuts_NPhotons25"]->Fill( et25 ); 
    histos_["PreCuts_NPhotons30"]->Fill( et30 ); 
    histos_["PreCuts_NPhotons35"]->Fill( et35 ); 
    histos_["PreCuts_NPhotons40"]->Fill( et40 ); 
    histos_["PreCuts_NPhotons45"]->Fill( et45 ); 
    histos_["PreCuts_NPhotons50"]->Fill( et50 ); 
    
    if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Photons: " << photons.size(); }

    histos_["PreCuts_NPhotons"]->Fill( photons.size() ); 
    
  }
  
  return false; 

}

// -----------------------------------------------------------------------------
//
bool TestCombination::getJets( const edm::Event& iEvent,
			       std::vector<reco::Particle>& jets ) {

  jets.clear();
  
  if ( !jets_.label().empty() ) {

    edm::Handle< std::vector<pat::Jet> > handle;
    iEvent.getByLabel(jets_,handle);
    
    if ( !handle.isValid() ) { 
      edm::LogWarning("TEST") << "No Jets for " << jets_; 
      return true;
    }
    
    uint16_t et25=0, et30=0, et35=0, et40=0, et45=0, et50=0;
    std::vector<pat::Jet>::const_iterator ijet = handle->begin();
    std::vector<pat::Jet>::const_iterator jjet = handle->end();
    for ( ; ijet != jjet; ++ijet  ) { 
      if ( fabs( ijet->eta() ) < jetEta_ && 
	   ijet->emEnergyFraction() < jetEMfrac_ ) { 
	if ( ijet->et() > jetEt_ ) { jets.push_back( *ijet ); }
	if ( ijet->et() > 25. ) { et25++; }
	if ( ijet->et() > 30. ) { et30++; }
	if ( ijet->et() > 35. ) { et35++; }
	if ( ijet->et() > 40. ) { et40++; } 
	if ( ijet->et() > 45. ) { et45++; }
	if ( ijet->et() > 50. ) { et50++; }
      }
    }
    
    histos_["PreCuts_NJets25"]->Fill( et25 ); 
    histos_["PreCuts_NJets30"]->Fill( et30 ); 
    histos_["PreCuts_NJets35"]->Fill( et35 ); 
    histos_["PreCuts_NJets40"]->Fill( et40 ); 
    histos_["PreCuts_NJets45"]->Fill( et45 ); 
    histos_["PreCuts_NJets50"]->Fill( et50 ); 

    if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Jets: " << jets.size(); }

    histos_["PreCuts_NJets"]->Fill( jets.size() ); 

  }
  
  return false;
  
}

// -----------------------------------------------------------------------------
//
bool TestCombination::getMuons( const edm::Event& iEvent,
				std::vector<reco::Particle>& muons ) {
  
  muons.clear();
  
  if ( !muons_.label().empty() ) {

    edm::Handle< std::vector<pat::Muon> > handle;
    iEvent.getByLabel(muons_,handle);
    
    if ( !handle.isValid() ) { 
      edm::LogWarning("TEST") << "No Muons for " << muons_; 
      return true;
    }
    
    std::vector<pat::Muon>::const_iterator imuon = handle->begin();
    std::vector<pat::Muon>::const_iterator jmuon = handle->end();
    for ( ; imuon != jmuon; ++imuon  ) { 
      if ( fabs( imuon->eta() ) < muonEta_ &&
	   imuon->trackIso() < muonTrkIso_ ) { 
	if ( imuon->pt() > muonPt_ ) { muons.push_back( *imuon ); }
      }
    }

    if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Muons: " << muons.size(); }

    histos_["PreCuts_NMuons"]->Fill( muons.size() ); 

  }

  return false;
  
}

// -----------------------------------------------------------------------------
//
bool TestCombination::getElectrons( const edm::Event& iEvent,
				    std::vector<reco::Particle>& electrons ) {
  
  electrons.clear();
  
  if ( !electrons_.label().empty() ) {

    edm::Handle< std::vector<pat::Electron> > handle;
    iEvent.getByLabel(electrons_,handle);

    if ( !handle.isValid() ) { 
      edm::LogWarning("TEST") << "No Electrons for " << electrons_; 
      return true;
    }
    
    std::vector<pat::Electron>::const_iterator ielectron = handle->begin();
    std::vector<pat::Electron>::const_iterator jelectron = handle->end();
    for ( ; ielectron != jelectron; ++ielectron  ) { 
      if ( fabs( ielectron->eta() ) < electronEta_ &&
	   ielectron->trackIso() < electronTrkIso_ ) { 
	if ( ielectron->pt() > electronPt_ ) { electrons.push_back( *ielectron ); }
      }
    }

    if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Electrons: " << electrons.size(); }

    histos_["PreCuts_NElectrons"]->Fill( electrons.size() ); 

  }
  
  return false;
  
}

// -----------------------------------------------------------------------------
//
float TestCombination::minDeltaEt( const std::vector<reco::Particle>& objects, 
				   std::vector<reco::Particle>& jet1, 
				   std::vector<reco::Particle>& jet2  ) {
  
  // Maximum number of objects handled
  static const uint8_t max_size = 50;

  // Init
  jet1.clear();
  jet2.clear();
  float min_delta_et = -1.;
  uint32_t combinations = 0;
  
  // Set number of objects used in combinations (limited to "max_size")
  uint8_t maximum = maximum_; 
  if ( maximum_ < 0 || maximum_ > max_size ) { maximum = max_size; } 
  uint8_t size = objects.size();
  if ( size > maximum ) { 
    edm::LogWarning("TEST")
      << " Number of objects (" << size
      << ") exceeds maximum set (" << maximum
      << ")!";
  }
  
  if ( size == 1 ) { 
    
    // Only one jet, in first pseudo-jet, second is empty
    jet1.push_back( objects.front() );
    jet2.clear();
    min_delta_et = 0.;
    combinations = 1;
    
  } else if ( size > 0 ) {

    // Cached data
    float min_diff_et = -1.;
    std::vector<uint8_t> indices1;
    std::vector<uint8_t> indices2;
    
    // Build char array encoding indices of all jets
    std::vector<char> vc1;
    vc1.reserve(size+1);
    for ( uint8_t ii = 0; ii < size; ++ii ) { vc1.push_back(ii); }
    vc1[vc1.size()] = '\0';
    char c1[256];
    memcpy( c1, &vc1.front(), vc1.size() );

    // Iterate through first half of jets
    for ( uint8_t jj = 0; jj < size/2; ++jj ) { 
      
      // Build char array encoding indices of subset of jets (to maximum of half)
      std::vector<char> vc2;
      vc2.reserve(jj+1+1);
      for ( uint8_t jjj = 0; jjj < jj+1; ++jjj ) { vc2.push_back(jjj); }
      vc2[vc2.size()] = '\0';
      char c2[256];
      memcpy( c2, &vc2.front(), vc2.size() );

      // Size of both char arrays
      uint8_t size1 = vc1.size();
      uint8_t size2 = vc2.size();

      // Iterate through combinations
      do { 

	// Indices of jets assigned to first pseudo-jet and Et sum
	float et1 = 0.;
	std::vector<uint8_t> tmp1;
	tmp1.reserve(size2+1);
	for ( uint8_t ii = 0; ii < size2; ++ii ) { 
	  uint8_t index = static_cast<uint8_t>( c2[ii] );
	  tmp1.push_back( index ); 
	  et1 += objects[index].et();
	}

	// Indices of jets assigned to second pseudo-jet and Et sum
	float et2 = 0.;
	std::vector<uint8_t> tmp2;
	tmp2.reserve(size1+1);
	for ( uint8_t ii = 0; ii < size1; ++ii ) { 
	  if ( std::find( tmp1.begin(), tmp1.end(), ii ) == tmp1.end() ) { 
	    tmp2.push_back(ii); 
	    et2 += objects[ii].et();
	  }
	}
	
	// Calculate difference in Et between two pseudo-jets
	if ( min_diff_et < 0. || fabs( et1 - et2 ) < min_diff_et ) { 
	  min_diff_et = fabs( et1 - et2 ); 
	  indices1.resize( tmp1.size() );
	  indices2.resize( tmp2.size() );
	  std::copy( tmp1.begin(), tmp1.end(), indices1.begin() );
	  std::copy( tmp2.begin(), tmp2.end(), indices2.begin() );
	}
	
	combinations++;

      } while ( stdcomb::next_combination( c1, c1 + size1, 
					   c2, c2 + size2 ) );
      
    }
    
    // Build pseudo-jets
    if ( min_diff_et < 0. ) {
      jet1.clear();
      jet2.clear();
    } else {
      std::vector<uint8_t>::const_iterator ii1 = indices1.begin();
      std::vector<uint8_t>::const_iterator jj1 = indices1.end();
      for ( ; ii1 != jj1; ++ii1 ) { jet1.push_back( objects[*ii1] ); }
      std::vector<uint8_t>::const_iterator ii2 = indices2.begin();
      std::vector<uint8_t>::const_iterator jj2 = indices2.end();
      for ( ; ii2 != jj2; ++ii2 ) { jet2.push_back( objects[*ii2] ); }
    }

    min_delta_et = min_diff_et;

  } else {
    
    jet1.clear();
    jet2.clear();
    min_delta_et = -1.;
    
  }

  return min_delta_et;
  
}

// -----------------------------------------------------------------------------
//
float TestCombination::sumEt( const std::vector<reco::Particle>& input ) {
  float et = 0.;
  std::vector<reco::Particle>::const_iterator ii = input.begin();
  std::vector<reco::Particle>::const_iterator jj = input.end();
  for ( ; ii != jj; ++ii ) { et += ii->et(); }
  return et;
}

// -----------------------------------------------------------------------------
//
float TestCombination::massT( const std::vector<reco::Particle>& input ) {
  float et = 0.; 
  float px = 0.; 
  float py = 0.; 
  std::vector<reco::Particle>::const_iterator ii = input.begin();
  std::vector<reco::Particle>::const_iterator jj = input.end();
  for ( ; ii != jj; ++ii ) { 
    et += ii->et();
    px += ii->px();
    py += ii->py();
  }
  float mt_sq = et*et - px*px - py*py;
  return mt_sq < 0. ? -1.*sqrt(-1.*mt_sq) : sqrt(mt_sq);
}

// -----------------------------------------------------------------------------
//
float TestCombination::alphaT( float min_et, float sum_et, float mass_t ) {
  return 0.5 * ( ( sum_et - min_et ) / mass_t );
}

// -----------------------------------------------------------------------------
//
float TestCombination::alphaT( const std::vector<reco::Particle>& objects,
			       std::vector<reco::Particle>& jet1, 
			       std::vector<reco::Particle>& jet2 ) {
  jet1.clear();
  jet2.clear();
  float min_et = minDeltaEt( objects, jet1, jet2 );
  float sum_et = sumEt( objects );
  float mass_t = massT( objects );
  return alphaT( min_et, sum_et, mass_t );
}
  
// -----------------------------------------------------------------------------
// 
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TestCombination);












//     {
//       std::vector<int>::const_iterator iii1 = phot1.begin();
//       std::vector<int>::const_iterator jjj1 = phot1.end();
//       for ( ; iii1 != jjj1; ++iii1 ) {
// 	reco::Candidate* candidate = 0;
// 	int32_t invalid_id = static_cast<int32_t>(1e8+0.5);
// 	int32_t daughter_id = invalid_id;
// 	int32_t mother_id   = invalid_id;
// 	std::vector<reco::GenParticleRef> gen_particles = iphoton->genParticleRefs();
// 	std::vector<reco::GenParticleRef>::const_iterator igen = gen_particles.begin();
// 	std::vector<reco::GenParticleRef>::const_iterator jgen = gen_particles.end();
// 	for ( ; igen != jgen; ++igen ) {
// 	  if ( igen->isNonnull() ) { 
// 	    reco::GenParticle* part = const_cast<reco::GenParticle*>( &(**igen) );
// 	    reco::Candidate* cand = dynamic_cast<reco::Candidate*>( part );
// 	    if ( cand ) {
// 	      candidate = cand;
// 	      daughter_id = cand->pdgId();
// 	      const reco::Candidate* mother = cand->mother();
// 	      if ( mother && cand->pdgId() == cand->mother()->pdgId() ) { mother = mother->mother(); }
// 	      if ( mother ) { mother_id = mother->pdgId(); }
// 	      break;
// 	    }
// 	  }
// 	}
// 	histos2d_["PseudoDijets_ParentPdgId"]->Fill( static_cast<uint32_t>( iii1 - phot1.begin() ), mother_id ); 
//       }
//     }

//     {
//       std::vector<int>::const_iterator iii2 = phot2.begin();
//       std::vector<int>::const_iterator jjj2 = phot2.end();
//       for ( ; iii2 != jjj2; ++iii2 ) {
// 	reco::Candidate* candidate = 0;
// 	int32_t invalid_id = static_cast<int32_t>(1e8+0.5);
// 	int32_t daughter_id = invalid_id;
// 	int32_t mother_id   = invalid_id;
// 	std::vector<reco::GenParticleRef> gen_particles = iphoton->genParticleRefs();
// 	std::vector<reco::GenParticleRef>::const_iterator igen = gen_particles.begin();
// 	std::vector<reco::GenParticleRef>::const_iterator jgen = gen_particles.end();
// 	for ( ; igen != jgen; ++igen ) {
// 	  if ( igen->isNonnull() ) { 
// 	    reco::GenParticle* part = const_cast<reco::GenParticle*>( &(**igen) );
// 	    reco::Candidate* cand = dynamic_cast<reco::Candidate*>( part );
// 	    if ( cand ) {
// 	      candidate = cand;
// 	      daughter_id = cand->pdgId();
// 	      const reco::Candidate* mother = cand->mother();
// 	      if ( mother && cand->pdgId() == cand->mother()->pdgId() ) { mother = mother->mother(); }
// 	      if ( mother ) { mother_id = mother->pdgId(); }
// 	      break;
// 	    }
// 	  }
// 	}
// 	histos2d_["PseudoDijets_ParentPdgId"]->Fill( static_cast<uint32_t>( iii2 - phot2.begin() ), mother_id ); 
//       }
//     }
    


// // -----------------------------------------------------------------------------
// //
// float TestCombination::constructPseudoJets( const std::vector<reco::Particle>& objects, 
// 					    std::vector<reco::Particle>& pseudo_jet1, 
// 					    std::vector<reco::Particle>& pseudo_jet2  ) {
  
//   // Build combinations
//   std::vector< std::vector<uint8_t> > combination1;
//   std::vector< std::vector<uint8_t> > combination2;
//   combinations( objects.size(), 
// 		combination1, 
// 		combination2 );
  
//   // Calc min Deltaet
//   return minDeltaet( objects, 
// 		     combination1, 
// 		     combination2, 
// 		     pseudo_jet1, 
// 		     pseudo_jet2 );
  
// }

// // -----------------------------------------------------------------------------
// // 
// float TestCombination::minDeltaet( const std::vector<reco::Particle>& objects, 
// 				   const std::vector< std::vector<uint8_t> >& combo1, 
// 				   const std::vector< std::vector<uint8_t> >& combo2,
// 				   std::vector<reco::Particle>& pseudo1, 
// 				   std::vector<reco::Particle>& pseudo2 ) {
  
//   pseudo1.clear();
//   pseudo2.clear();
  
//   uint32_t size = combo1.size() < combo2.size() ? combo1.size() : combo2.size();
  
//   if ( size == 0 ) { return -1.; }
  
//   float min_delta_pt = -1.;
//   uint32_t index = 0;
//   for ( uint32_t ii = 0; ii < size; ++ii ) {
//     float pt1 = 0.;
//     for ( uint32_t jj = 0; jj < combo1[ii].size(); ++jj ) { pt1 += objects[ combo1[ii][jj] ].pt(); }
//     float pt2 = 0.;
//     for ( uint32_t jj = 0; jj < combo2[ii].size(); ++jj ) { pt2 += objects[ combo2[ii][jj] ].pt(); }
//     if ( min_delta_pt < 0. ) { 
//       min_delta_pt = fabs( pt1 - pt2 ); 
//       index = ii;
//     }
//     else { 
//       if ( fabs( pt1 - pt2 ) < min_delta_pt ) { 
// 	min_delta_pt = fabs( pt1 - pt2 ); 
// 	index = ii;
//       } 
//     }
//   }
  
//   for ( uint32_t jj = 0; jj < combo1[index].size(); ++jj ) { pseudo1.push_back( objects[ combo1[index][jj] ] ); }
//   for ( uint32_t jj = 0; jj < combo2[index].size(); ++jj ) { pseudo2.push_back( objects[ combo2[index][jj] ] ); }

//   return min_delta_pt;
  
// }

// // -----------------------------------------------------------------------------
// //
// void TestCombination::combinations( uint8_t size, 
// 				    std::vector< std::vector<uint8_t> >& combo1, 
// 				    std::vector< std::vector<uint8_t> >& combo2  ) {
  
//   static uint8_t max_size = 50;
  
//   combo1.clear();
//   combo2.clear();
  
//   // Set number of objects used in combinations (limited to "max_size")
//   if ( maximum_ < 0 || maximum_ > max_size ) { maximum_ = max_size; } 
//   size = size < static_cast<uint8_t>(maximum_) ? size : static_cast<uint16_t>(maximum_); 

//   uint32_t capacity = static_cast<uint32_t>( pow(2,size-1) + 0.5 ) * size;
  
// //   if ( capacity > combo1.max_size() ||
// //        capacity > combo2.max_size() ) {
// //     size = 10;
// //     capacity = static_cast<uint32_t>( pow(2,size-1) + 0.5 ) * size/2;
// //     edm::LogWarning("TEST")
// //       << " Capacity requested exceeds \"max_size()\" of vector!"
// //       << " Setting size to " << size << " and capacity to " << capacity;
// //   }
  
//   std::cout << std::endl
// 	    << "test1 " << std::endl 
// 	    << " size " << static_cast<uint16_t>(size) << std::endl
// 	    << " capacity " << capacity << std::endl
// 	    << " combo1.max_size() " << combo1.max_size() << std::endl;
  
//   //combo1.reserve( capacity );
//   //combo2.reserve( capacity );
  
//   std::cout << "test2 " << std::endl 
// 	    << " size " << static_cast<uint16_t>(size) << std::endl
// 	    << " combo1.capacity() " << combo1.capacity() << std::endl;
  
//   if ( size == 1 ) { 
    
//     combo1.push_back( std::vector<uint8_t>(1,0) );
//     combo2.push_back( std::vector<uint8_t>() );
    
//   } else if ( size > 0 ) {

//     uint32_t cap1 = 0;
//     uint32_t cap2 = 0;
    
//     std::vector<char> vc1;
//     vc1.reserve(size+1);
//     for ( uint8_t ii = 0; ii < size; ++ii ) { vc1.push_back(ii); }
//     vc1[vc1.size()] = '\0';

//     char c1[256];
//     memcpy( c1, &vc1.front(), vc1.size() );
 
//     for ( uint8_t jj = 0; jj < size/2; ++jj ) { 
      
//       std::vector<char> vc2;
//       vc2.reserve(jj+1+1);
//       for ( uint8_t jjj = 0; jjj < jj+1; ++jjj ) { vc2.push_back(jjj); }
//       vc2[vc2.size()] = '\0';
      
//       char c2[256];
//       memcpy( c2, &vc2.front(), vc2.size() );
      
//       uint8_t size1 = vc1.size();
//       uint8_t size2 = vc2.size();
      
//       do { 
// 	std::vector<uint8_t> tmp1;
// 	tmp1.reserve(size2+1);
// 	for ( uint8_t ii = 0; ii < size2; ++ii ) { tmp1.push_back( static_cast<uint8_t>(c2[ii]) ); }
// 	combo1.push_back( tmp1 );
// 	cap1 += combo1.back().capacity();
// 	std::vector<uint8_t> tmp2;
// 	tmp2.reserve(size1+1);
// 	for ( uint8_t ii = 0; ii < size1; ++ii ) { 
// 	  if ( std::find( tmp1.begin(), tmp1.end(), ii ) == tmp1.end() ) { tmp2.push_back(ii); }
// 	}
// 	combo2.push_back( tmp2 );
// 	cap2 += combo2.back().capacity();
//       }
//       while ( stdcomb::next_combination( c1, c1 + size1, 
// 					 c2, c2 + size2 ) );
      
//     }

//     std::cout << "test3 " << std::endl 
// 	      << " size " << static_cast<uint16_t>(size) << std::endl 
// 	      << " combo1.capacity() " << combo1.capacity() << std::endl 
// 	      << " combo2.capacity() " << combo2.capacity() << std::endl 
// 	      << " combo1.size() " << combo1.size() << std::endl 
// 	      << " combo2.size() " << combo2.size() << std::endl
// 	      << " cap1 " << cap1 << std::endl 
// 	      << " cap2 " << cap2 << std::endl; 

//     if ( size > 0 && size % 2 == 0 ) {
//       uint32_t numpop = factorial(size) / (factorial(size/2) * factorial(size/2) );
//       for ( uint32_t i = 0; i < numpop / 2; i++ ) {
// 	combo1.pop_back();
// 	combo2.pop_back();
//       }
//     }

//   }
  
//   if ( edm::isDebugEnabled() ) {
//     std::stringstream ss;
//     std::stringstream sss;
//     ss << "Number of primary objects = " << static_cast<uint16_t>(size) << std::endl;
//     ss << "Number of combinations for first and second pseudo-jets = " 
//        << combo1.size()
//        << " + "
//        << combo2.size()
//        << std::endl;
//     uint32_t cap1 = 0;
//     uint32_t cap2 = 0;
//     uint32_t combinations = combo1.size() < combo2.size() ? combo1.size() : combo2.size();
//     for ( uint32_t ii = 0; ii < combinations; ++ii ) { 
//       sss << " Combination #" << ii
// 	  << ", size = " << combo1[ii].size()
// 	  << " + " << combo2[ii].size();
//       sss << ", object indices = ";
//       cap1 += combo1[ii].size();
//       cap2 += combo2[ii].size();
//       for ( uint8_t jj = 0; jj < combo1[ii].size(); ++jj ) { 
// 	if ( jj != 0 ) { sss << ","; }
// 	sss << static_cast<uint16_t>(combo1[ii][jj]);
//       }
//       sss << " + ";
//       for ( uint8_t jj = 0; jj < combo2[ii].size(); ++jj ) { 
// 	if ( jj != 0 ) { sss << ","; }
// 	sss << static_cast<uint16_t>(combo2[ii][jj]);
//       }
//       sss << std::endl;
//     }
//     edm::LogVerbatim("TEST") << ss.str();
//     LogTrace("TEST") << sss.str();
//     std::cout << "test4 " << std::endl 
// 	      << " size " << static_cast<uint16_t>(size) << std::endl 
// 	      << " combo1.capacity() " << combo1.capacity() << std::endl 
// 	      << " combo2.capacity() " << combo2.capacity() << std::endl 
// 	      << " combo1.size() " << combo1.size() << std::endl 
// 	      << " combo2.size() " << combo2.size() << std::endl 
// 	      << " cap1 " << cap1 << std::endl 
// 	      << " cap2 " << cap2 << std::endl; 
//   }
  
// }
