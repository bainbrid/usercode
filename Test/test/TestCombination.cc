#include "TestCombination.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/Utilities/interface/PtComparator.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "Combination.h"
#include "TH1D.h"
#include "TH2D.h"
#include <algorithm>
#include <sstream>
#include <cstring>

// -----------------------------------------------------------------------------
//
TestCombination::TestCombination( const edm::ParameterSet& pset ) 
  : jets_( pset.getParameter<edm::InputTag>("JetCollection") ),
    photons_( pset.getParameter<edm::InputTag>("PhotonCollection") ),
    jetPt_( pset.getParameter<double>("JetPt") ),
    jetEta_( pset.getParameter<double>("JetEta") ),
    photonPt_( pset.getParameter<double>("PhotonPt") ),
    photonEta_( pset.getParameter<double>("PhotonEta") ),
    totalPt_( pset.getParameter<double>("TotalPt") ),
    histos_(),
    histos2d_()
{;}

// -----------------------------------------------------------------------------
//
void TestCombination::beginJob( const edm::EventSetup& ) {
  
  edm::Service<TFileService> fs;
  
  histos_["NJets"] = fs->make<TH1D>("NJets","",21,-0.5,20.5);
  histos_["NPhotons"] = fs->make<TH1D>("NPhotons","",21,-0.5,20.5);
  histos_["NObjects"] = fs->make<TH1D>("NObjects","",21,-0.5,20.5);
  histos_["MinDeltaPt"] = fs->make<TH1D>("MinDeltaPt","",100,0.,100.);
  histos_["TotalPt"] = fs->make<TH1D>("TotalPt","",100,0.,1000.);

  histos_["NObjectsInPseudoJet1"] = fs->make<TH1D>("NObjectsInPseudoJet1","",21,-0.5,20.5);
  histos_["NObjectsInPseudoJet2"] = fs->make<TH1D>("NObjectsInPseudoJet2","",21,-0.5,20.5);
  histos_["NObjectsInPseudoJets"] = fs->make<TH1D>("NObjectsInPseudoJets","",21,-0.5,20.5);

  histos2d_["MinDeltaPtVsNObjects"] = fs->make<TH2D>("MinDeltaPtVsNObjects","",21,-0.5,20.5,100,0.,100.);

  histos_["PhotonsInBothPseudoJets_NPhotons1"] = fs->make<TH1D>("PhotonsInBothPseudoJets_NPhotons1","",21,-0.5,20.5);
  histos_["PhotonsInBothPseudoJets_NPhotons2"] = fs->make<TH1D>("PhotonsInBothPseudoJets_NPhotons2","",21,-0.5,20.5);

  histos_["BothPhotonsLeading_PhotonPt1"] = fs->make<TH1D>("BothPhotonsLeading_PhotonPt1","",100,0.,1000.);
  histos_["BothPhotonsLeading_PhotonPt2"] = fs->make<TH1D>("BothPhotonsLeading_PhotonPt2","",100,0.,1000.);
  histos_["BothPhotonsLeading_DeltaPt"] = fs->make<TH1D>("BothPhotonsLeading_DeltaPt","",100,0.,100.);
  histos_["BothPhotonsLeading_DeltaEta"] = fs->make<TH1D>("BothPhotonsLeading_DeltaEta","",100,0.,5.);
  
  histos_["BothPhotonsNotLeading_PhotonPt1"] = fs->make<TH1D>("BothPhotonsNotLeading_PhotonPt1","",100,0.,1000.);
  histos_["BothPhotonsNotLeading_PhotonPt2"] = fs->make<TH1D>("BothPhotonsNotLeading_PhotonPt2","",100,0.,1000.);
  histos_["BothPhotonsNotLeading_DeltaPt"] = fs->make<TH1D>("BothPhotonsNotLeading_DeltaPt","",100,0.,100.);
  
  histos_["BothPhotonsNotLeading_DeltaEta"] = fs->make<TH1D>("BothPhotonsNotLeading_DeltaEta","",100,0.,5.);

  histos_["PhotonsInPseudoJet1_DeltaEta"] = fs->make<TH1D>("PhotonsInPseudoJet1_DeltaEta","",100,0.,5.);
  histos_["PhotonsInPseudoJet2_DeltaEta"] = fs->make<TH1D>("PhotonsInPseudoJet2_DeltaEta","",100,0.,5.);
  
}

// -----------------------------------------------------------------------------
//
void TestCombination::analyze( const edm::Event& iEvent, 
			       const edm::EventSetup& iSetup ) {
  
  // Get jets
  std::vector<reco::Particle> jets;
  if ( !jets_.label().empty() ) {
    edm::Handle< std::vector<pat::Jet> > handle;
    iEvent.getByLabel(jets_,handle);
    if ( !handle.isValid() ) { edm::LogWarning("TEST") << "No Jets for InputTag " << jets_; }
    std::vector<pat::Jet>::const_iterator ijet = handle->begin();
    std::vector<pat::Jet>::const_iterator jjet = handle->end();
    for ( ; ijet != jjet; ++ijet  ) { 
      if ( ijet->pt() > jetPt_ && fabs( ijet->eta() ) < jetEta_ ) { jets.push_back( *ijet ); }
    }
    LogTrace("TEST") << "Number of Jets: " << jets.size();
  }
  
  // Get photons
  std::vector<reco::Particle> photons;
  if ( !photons_.label().empty() ) {
    edm::Handle< std::vector<pat::Photon> > handle;
    iEvent.getByLabel(photons_,handle);
    if ( !handle.isValid() ) { edm::LogWarning("TEST") << "No Photons for InputTag " << photons_; }
    std::vector<pat::Photon>::const_iterator iphoton = handle->begin();
    std::vector<pat::Photon>::const_iterator jphoton = handle->end();
    for ( ; iphoton != jphoton; ++iphoton  ) { 
      if ( iphoton->et() > photonPt_  && fabs( iphoton->eta() ) < photonEta_ ) { photons.push_back( *iphoton ); }
    }
    LogTrace("TEST") << "Number of Photons: " << photons.size();
  }
  
  // All candidates
  std::vector<reco::Particle> particles;
  if ( !jets.empty() ) { particles.insert( particles.end(), jets.begin(), jets.end() ); }
  if ( !photons.empty() ) { particles.insert( particles.end(), photons.begin(), photons.end() ); }
  
  // Check pt total
  float pt = 0.;
  std::vector<reco::Particle>::const_iterator iparticle = particles.begin();
  std::vector<reco::Particle>::const_iterator jparticle = particles.end();
  for ( ; iparticle != jparticle; ++iparticle  ) { pt += iparticle->pt(); }
  if ( pt < totalPt_ ) { return; }
  
  // Build pseudo-jets
  std::vector<reco::Particle> pseudo_jet1;
  std::vector<reco::Particle> pseudo_jet2;
  float min_pt = pseudoJets( particles, pseudo_jet1, pseudo_jet2 );
  
  // Histos
  int size1 = pseudo_jet1.size();
  int size2 = pseudo_jet2.size();
  int size = size1 + size2;

  histos_["NJets"]->Fill( jets.size() ); 
  histos_["NPhotons"]->Fill( photons.size() ); 
  histos_["NObjects"]->Fill( particles.size() ); 
  histos_["MinDeltaPt"]->Fill( min_pt ); 
  histos_["TotalPt"]->Fill( pt ); 
  
  histos_["NObjectsInPseudoJet1"]->Fill( size1 ); 
  histos_["NObjectsInPseudoJet2"]->Fill( size2 ); 
  histos_["NObjectsInPseudoJets"]->Fill( size ); 

  histos2d_["MinDeltaPtVsNObjects"]->Fill( size, min_pt ); 
  
  std::vector<int> phot1;
  std::vector<reco::Particle>::const_iterator ii1 = pseudo_jet1.begin();
  std::vector<reco::Particle>::const_iterator jj1 = pseudo_jet1.end();
  for ( ; ii1 != jj1; ++ii1 ) { if ( ii1->pdgId() == 22 ) { phot1.push_back( ii1 - pseudo_jet1.begin() ); } }
  
  std::vector<int> phot2;
  std::vector<reco::Particle>::const_iterator ii2 = pseudo_jet2.begin();
  std::vector<reco::Particle>::const_iterator jj2 = pseudo_jet2.end();
  for ( ; ii2 != jj2; ++ii2 ) { if ( ii2->pdgId() == 22 ) { phot2.push_back( ii2 - pseudo_jet2.begin() ); } }
  
  if ( !phot1.empty() && !phot2.empty() ) {
    
    histos_["PhotonsInBothPseudoJets_NPhotons1"]->Fill( phot1.size() ); 
    histos_["PhotonsInBothPseudoJets_NPhotons2"]->Fill( phot2.size() ); 
    
    if ( phot1[0] == 0 &&phot2[0] == 0 ) {
      histos_["BothPhotonsLeading_PhotonPt1"]->Fill( pseudo_jet1[ phot1[0] ].pt() ); 
      histos_["BothPhotonsLeading_PhotonPt2"]->Fill( pseudo_jet2[ phot2[0] ].pt() ); 
      histos_["BothPhotonsLeading_DeltaPt"]->Fill( fabs( pseudo_jet1[ phot1[0] ].pt() - 
							 pseudo_jet2[ phot2[0] ].pt() ) ); 
      histos_["BothPhotonsLeading_DeltaEta"]->Fill( fabs( pseudo_jet1[ phot1[0] ].eta() - 
							  pseudo_jet2[ phot2[0] ].eta() ) ); 
    } else {
      histos_["BothPhotonsNotLeading_PhotonPt1"]->Fill( pseudo_jet1[ phot1[0] ].pt() ); 
      histos_["BothPhotonsNotLeading_PhotonPt2"]->Fill( pseudo_jet2[ phot2[0] ].pt() ); 
      histos_["BothPhotonsNotLeading_DeltaPt"]->Fill( fabs( pseudo_jet1[ phot1[0] ].pt() - 
							    pseudo_jet2[ phot2[0] ].pt() ) ); 
      histos_["BothPhotonsNotLeading_DeltaEta"]->Fill( fabs( pseudo_jet1[ phot1[0] ].eta() - 
							     pseudo_jet2[ phot2[0] ].eta() ) ); 
    }

  } else if ( phot1.size() >= 2 ) {
    histos_["PhotonsInPseudoJet1_DeltaEta"]->Fill( fabs( pseudo_jet1[ phot1[0] ].eta() - 
							 pseudo_jet1[ phot1[1] ].eta() ) ); 
  } else if ( phot2.size() >= 2 ) {
    histos_["PhotonsInPseudoJet2_DeltaEta"]->Fill( fabs( pseudo_jet2[ phot2[0] ].eta() - 
							 pseudo_jet2[ phot2[1] ].eta() ) ); 
  }
  
  // Debug
  if ( edm::isDebugEnabled() ) {
    std::stringstream ss;
    ss << "Minimum deltaPt = " << min_pt << std::endl;

    {
      ss << "Input objects: id/pdg: ";
      std::vector<reco::Particle>::const_iterator ii = particles.begin();
      std::vector<reco::Particle>::const_iterator jj = particles.end();
      for ( ; ii != jj; ++ii  ) { 
	ss << ii->pdgId() << "/" << ii->pt() << ", "; 
      }
      ss << std::endl;
    }

    {
      ss << "Pseudo-jet #1: id/pdg: ";
      std::vector<reco::Particle>::const_iterator ii = pseudo_jet1.begin();
      std::vector<reco::Particle>::const_iterator jj = pseudo_jet1.end();
      for ( ; ii != jj; ++ii  ) { 
	ss << ii->pdgId() << "/" << ii->pt() << ", "; 
      }
      ss << std::endl;
    }

    {
      ss << "Pseudo-jet #2: id/pdg: ";
      std::vector<reco::Particle>::const_iterator ii = pseudo_jet2.begin();
      std::vector<reco::Particle>::const_iterator jj = pseudo_jet2.end();
      for ( ; ii != jj; ++ii  ) { 
	ss << ii->pdgId() << "/" << ii->pt() << ", "; 
      }
      ss << std::endl;
    }
    
    LogTrace("TEST") << ss.str();

  }
  
  
}

// -----------------------------------------------------------------------------
//
float TestCombination::pseudoJets( const std::vector<reco::Particle>& all_particles, 
				   std::vector<reco::Particle>& pseudo_jet1, 
				   std::vector<reco::Particle>& pseudo_jet2  ) {
  
  // Make copy of particles
  std::vector<reco::Particle> particles;
  particles.resize( all_particles.size() );
  std::copy( all_particles.begin(), all_particles.end(), particles.begin() );
  sort( particles.begin(), particles.end(), GreaterByPt<reco::Particle>() );
  
  // Build combinations
  std::vector< std::vector<int> > combination1;
  std::vector< std::vector<int> > combination2;
  combinations( particles.size(), combination1, combination2 );

  // Calc min DeltaPt
  pseudo_jet1.clear();
  pseudo_jet2.clear();
  return minDeltaPt( particles, 
		     combination1, 
		     combination2, 
		     pseudo_jet1, 
		     pseudo_jet2 );
  
}

// -----------------------------------------------------------------------------
//
void TestCombination::combinations( int size, 
				    std::vector< std::vector<int> >& combo1, 
				    std::vector< std::vector<int> >& combo2  ) {
  
  size = size < 10 ? size : 10;
  
  combo1.clear();
  combo2.clear();
  
  if ( size == 0 ) { return; }
  if ( size == 1 ) { 
    combo1.push_back( std::vector<int>(1,0) );
    combo2.push_back( std::vector<int>() );
    return; 
  }
  
  std::stringstream ss1;
  for ( uint16_t ii = 0; ii < size; ++ii ) { ss1 << ii; }
  
  std::stringstream ss2;
  for ( uint16_t jj = 0; jj < size/2; ++jj ) { ss2 << jj; }
  
  char cc1[100];
  std::strcpy( cc1, ss1.str().c_str() );
  
  char cc2[100];
  std::strcpy( cc2, ss2.str().c_str() );

  int size1 = ss1.str().size();
  int size2 = ss2.str().size();
  
  do { 
    std::vector<int> tmp1;
    for ( int ii = 0; ii < size2; ++ii ) { 
      std::stringstream sss; sss << cc2[ii];
      int num; sss >> num;
      tmp1.push_back( num ); 
    }
    combo1.push_back( tmp1 );
    std::vector<int> tmp2;
    for ( int ii = 0; ii < size1; ++ii ) { 
      if ( std::find( tmp1.begin(), tmp1.end(), ii ) == tmp1.end() ) { 
	tmp2.push_back(ii); 
      }
    }
    combo2.push_back( tmp2 );
  }
  while ( stdcomb::next_combination( cc1, cc1 + size1, 
				     cc2, cc2 + size2 ) );

  if ( edm::isDebugEnabled() ) {
    std::stringstream ss;
    ss << " size=" << size << std::endl;
    int combo_size = combo1.size() < combo2.size() ? combo1.size() : combo2.size();
    ss << " combo1.size()=" << combo1.size() 
       << " combo2.size()=" << combo2.size()
       << std::endl;
    for ( uint16_t ii = 0; ii < combo_size; ++ii ) { 
      ss << " combo1[" << ii << "].size()=" << combo1[ii].size()
	 << ", values=";
      for ( uint16_t jj = 0; jj < combo1[ii].size(); ++jj ) { 
	ss << combo1[ii][jj] << ",";
      }
      ss << " combo2[" << ii << "].size()=" << combo2[ii].size()
	 << ", values=";
      for ( uint16_t jj = 0; jj < combo2[ii].size(); ++jj ) { 
	ss << combo2[ii][jj] << ",";
      }
      ss << std::endl;
    }
    LogTrace("TEST") << ss.str();
  }
  
}

// -----------------------------------------------------------------------------
// 
float TestCombination::minDeltaPt( const std::vector<reco::Particle>& input, 
				   const std::vector< std::vector<int> >& combo1, 
				   const std::vector< std::vector<int> >& combo2,
				   std::vector<reco::Particle>& pseudo1, 
				   std::vector<reco::Particle>& pseudo2 ) {

  uint16_t size = combo1.size() < combo2.size() ? combo1.size() : combo2.size();

  if ( size == 0 ) { return -1.; }
  
  pseudo1.clear();
  pseudo2.clear();
  
  float min_delta_pt = -1.;
  int index = 0;
  for ( uint16_t ii = 0; ii < size; ++ii ) {
    float pt1 = 0.;
    for ( uint16_t jj = 0; jj < combo1[ii].size(); ++jj ) { pt1 += input[ combo1[ii][jj] ].pt(); }
    float pt2 = 0.;
    for ( uint16_t jj = 0; jj < combo2[ii].size(); ++jj ) { pt2 += input[ combo2[ii][jj] ].pt(); }
    if ( min_delta_pt < 0. ) { 
      min_delta_pt = fabs( pt1 - pt2 ); 
      index = ii;
    }
    else { 
      if ( fabs( pt1 - pt2 ) < min_delta_pt ) { 
	min_delta_pt = fabs( pt1 - pt2 ); 
	index = ii;
      } 
    }
  }
  
  for ( uint16_t jj = 0; jj < combo1[index].size(); ++jj ) { pseudo1.push_back( input[ combo1[index][jj] ] ); }
  for ( uint16_t jj = 0; jj < combo2[index].size(); ++jj ) { pseudo2.push_back( input[ combo2[index][jj] ] ); }

  return min_delta_pt;
  
}

// -----------------------------------------------------------------------------
// 
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TestCombination);
