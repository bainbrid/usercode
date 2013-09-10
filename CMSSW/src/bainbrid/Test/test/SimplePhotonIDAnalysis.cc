#include "bainbrid/Test/test/SimplePhotonIDAnalysis.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
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
int   SimplePhotonIDAnalysis::nbins_ = 30;
float SimplePhotonIDAnalysis::lower_ = 0.;
float SimplePhotonIDAnalysis::upper_ = 300;
int32_t SimplePhotonIDAnalysis::invalidId_  = static_cast<int32_t>(1e8+0.5);

// -----------------------------------------------------------------------------
//
SimplePhotonIDAnalysis::SimplePhotonIDAnalysis( const edm::ParameterSet& pset ) 
  : histos_(),
    histos2d_(),
    labels_(),
    photons_( pset.getParameter<edm::InputTag>("Photons") ),
    others_( pset.getParameter<edm::InputTag>("OtherPhotons") ),
    gen_( pset.getParameter<edm::InputTag>("GenParticles") ),
    genDaughterPdgId_(),
    genMotherPdgId_(),
    matchedDaughterPdgId_(),
    matchedMotherPdgId_(),
    phoIds_( pset.getParameter< std::vector<int32_t> >("PhotonPdgIds") ),
    eleIds_( pset.getParameter< std::vector<int32_t> >("ElectronPdgIds") ),
    allIds_( pset.getParameter< std::vector<int32_t> >("ParticlePdgIds") )
{
  labels_.push_back("SignalPhotons");
  labels_.push_back("BkgdPhotons");
  labels_.push_back("ElectronFakes");
  labels_.push_back("NonElectronFakes");
  labels_.push_back("Unmatched");
  genDaughterPdgId_.resize( labels_.size() );
  genMotherPdgId_.resize( labels_.size() );
  matchedDaughterPdgId_.resize( labels_.size() );
  matchedMotherPdgId_.resize( labels_.size() );
}

// -----------------------------------------------------------------------------
//
void SimplePhotonIDAnalysis::analyze( const edm::Event& event, 
				      const edm::EventSetup& setup ) {
  
//   // Photon counter
//   std::vector<uint16_t> multiplicity_reco( labels_.size(), 0 );
//   std::vector<uint16_t> multiplicity_truth( labels_.size(), 0 );
//   uint16_t              multiplicity_reco_all  = 0;
//   uint16_t              multiplicity_truth_all = 0;
  
//   // Retrieve MC matching maps for photons, electrons, all particles
//   edm::Handle< edm::Association<reco::GenParticleCollection> > matched_photons;
//   event.getByLabel( "matchedPhotons", matched_photons );
//   edm::Handle< edm::Association<reco::GenParticleCollection> > matched_electrons;
//   event.getByLabel( "matchedElectrons", matched_electrons );
//   edm::Handle< edm::Association<reco::GenParticleCollection> > matched_particles;
//   event.getByLabel( "matchedParticles", matched_particles );
  
//   // Retrieve PAT photon collection from Event
//   edm::Handle< edm::View<pat::Photon> > photons;
//   event.getByLabel( photons_, photons );
  
//   // Retrieve PAT PhotonID collection from Event
//   edm::Handle< edm::View<pat::Photon> > others;
//   event.getByLabel( others_, others );
  
//   // Retrieve GenParticle collection from Event
//   edm::Handle< edm::View<reco::GenParticle> > gen;
//   event.getByLabel( gen_, gen );
  
//   // Iterate though categories
//   for ( uint16_t ii = 0; ii < labels_.size(); ++ii ) {
    
//     std::string name = labels_[ii];

//     // Iterate though GenParticles
//     edm::View<reco::GenParticle>::const_iterator igen = gen->begin(); 
//     edm::View<reco::GenParticle>::const_iterator jgen = gen->end(); 
//     for ( ; igen != jgen; ++igen ) {

//       int32_t gen_daughter_id = invalidId_;
//       int32_t gen_mother_id   = invalidId_;
      
//       gen_daughter_id = igen->pdgId(); 
//       const reco::Candidate* mother = igen->mother();
//       if ( mother && igen->pdgId() == igen->mother()->pdgId() ) { mother = mother->mother(); }
//       if ( mother ) { gen_mother_id = mother->pdgId(); }

//       bool pho = phoIds_.empty() || find( phoIds_.begin(), phoIds_.end(), labs(gen_daughter_id) ) != phoIds_.end();
//       bool ele = eleIds_.empty() || find( eleIds_.begin(), eleIds_.end(), labs(gen_daughter_id) ) != eleIds_.end();
//       bool all = allIds_.empty() || find( allIds_.begin(), allIds_.end(), labs(gen_daughter_id) ) != allIds_.end();
      
//       if ( ii == 0 && pho ) { multiplicity_truth_all++; }
      
//       if ( !( ( ii == 0 && pho && gen_mother_id == 1000022 ) || //@@ SignalPhotons
// 	      ( ii == 1 && pho && gen_mother_id != 1000022 ) || //@@ BkgdPhotons
// 	      ( ii == 2 && ele ) ||                             //@@ ElectronFakes
// 	      ( ii == 3 && all && !pho && !ele ) ||             //@@ NonElectronicFakes
// 	      ( ii == 4 && gen_daughter_id == invalidId_ )      //@@ Unmatched
// 	      ) ) { continue; } 
      
//       genDaughterPdgId_[ii][gen_daughter_id]++; 
//       genMotherPdgId_[ii][gen_mother_id]++;
      
//       multiplicity_truth[ii]++;
//       histos_[name+"_GenPhotons"]->Fill( igen->et() ); 
      
//     }
    
//     // Iterate through PAT photons
//     edm::View<pat::Photon>::const_iterator iphoton = photons->begin(); 
//     edm::View<pat::Photon>::const_iterator jphoton = photons->end(); 
//     for ( ; iphoton != jphoton; ++iphoton ) {

//       float et = iphoton->et();
      
//       reco::Candidate* candidate = 0;
//       int32_t daughter_id = invalidId_;
//       int32_t mother_id   = invalidId_;
//       std::vector<reco::GenParticleRef> gen_particles = iphoton->genParticleRefs();
//       std::vector<reco::GenParticleRef>::const_iterator igen = gen_particles.begin();
//       std::vector<reco::GenParticleRef>::const_iterator jgen = gen_particles.end();
//       for ( ; igen != jgen; ++igen ) {
// 	if ( igen->isNonnull() ) { 
// 	  reco::GenParticle* part = const_cast<reco::GenParticle*>( &(**igen) );
//  	  reco::Candidate* cand = dynamic_cast<reco::Candidate*>( part );
//  	  if ( cand ) {
// 	    candidate = cand;
// 	    daughter_id = cand->pdgId();
// 	    const reco::Candidate* mother = cand->mother();
// 	    if ( mother && cand->pdgId() == cand->mother()->pdgId() ) { mother = mother->mother(); }
// 	    if ( mother ) { mother_id = mother->pdgId(); }
// 	    et = cand->et();
// 	    break;
// 	  }
// 	}
//       }
      
//       // Some selection
//       if ( fabs( iphoton->eta() ) > 2. ) { continue; }
//       //if ( et < 100. ) { continue; }

//       // Totals
//       histos_[name+"_TotalEMObjects"]->Fill( et ); 
//       if ( iphoton->isTightPhoton() ) { histos_[name+"_TotalTightPhotons"]->Fill( et ); }
//       if ( iphoton->isLoosePhoton() ) { histos_[name+"_TotalLoosePhotons"]->Fill( et ); }
//       if ( iphoton->isLooseEM() )     { histos_[name+"_TotalLooseEMObjects"]->Fill( et ); }

//       if ( ii == 0 ) { multiplicity_truth_all++; }

//       bool pho = phoIds_.empty() || find( phoIds_.begin(), phoIds_.end(), labs(daughter_id) ) != phoIds_.end();
//       bool ele = eleIds_.empty() || find( eleIds_.begin(), eleIds_.end(), labs(daughter_id) ) != eleIds_.end();
//       bool all = allIds_.empty() || find( allIds_.begin(), allIds_.end(), labs(daughter_id) ) != allIds_.end();

//       // DeltaR and DeltaPt plots for matched photons, electrons and all particles (done once only)
//       if ( ii == 0 && candidate ) { 
// 	if ( pho ) { 
// 	  float dr = deltaR<pat::Photon,reco::Candidate>( *iphoton, *candidate );
// 	  float pt = candidate->pt() != 0. ? fabs( ( iphoton->pt() - candidate->pt() ) / candidate->pt() ) : -1;
// 	  histos_["DeltaR_Photons"]->Fill( dr ); 
// 	  histos_["DeltaPt_Photons"]->Fill( pt ); 
// 	  histos2d_["DeltaRPt_Photons"]->Fill( dr, pt );
// 	} 
// 	if ( ele ) { 
// 	  histos_["DeltaR_Electrons"]->Fill( deltaR<pat::Photon,reco::Candidate>( *iphoton, *candidate ) ); 
// 	  histos_["DeltaPt_Electrons"]->Fill( fabs( iphoton->pt() - candidate->pt() ) ); 
// 	  histos2d_["DeltaRPt_Electrons"]->Fill( deltaR<pat::Photon,reco::Candidate>( *iphoton, *candidate ),
// 					       fabs( iphoton->pt() - candidate->pt() ) ); 
// 	}
// 	if ( all ) { 
// 	  histos_["DeltaR_Particles"]->Fill( deltaR<pat::Photon,reco::Candidate>( *iphoton, *candidate ) ); 
// 	  histos_["DeltaPt_Particles"]->Fill( fabs( iphoton->pt() - candidate->pt() ) ); 
// 	  histos2d_["DeltaRPt_Particles"]->Fill( deltaR<pat::Photon,reco::Candidate>( *iphoton, *candidate ),
// 						 fabs( iphoton->pt() - candidate->pt() ) ); 
// 	}
//       }
      
//       if ( !( ( ii == 0 && pho && mother_id == 1000022 ) || //@@ SignalPhotons
// 	      ( ii == 1 && pho && mother_id != 1000022 ) || //@@ BkgdPhotons
// 	      ( ii == 2 && ele ) ||                         //@@ ElectronFakes
// 	      ( ii == 3 && all && !pho && !ele ) ||         //@@ NonElectronicFakes
// 	      ( ii == 4 && daughter_id == invalidId_ )      //@@ Unmatched
// 	      ) ) { continue; } 
      
//       matchedDaughterPdgId_[ii][daughter_id]++; 
//       matchedMotherPdgId_[ii][mother_id]++; 

//       multiplicity_reco[ii]++; // Multiplicity
      
//       // Matched
//       histos_[name+"_MatchedEMObjects"]->Fill( et ); 
//       if ( iphoton->isTightPhoton() ) { histos_[name+"_MatchedTightPhotons"]->Fill( et ); }
//       if ( iphoton->isLoosePhoton() ) { histos_[name+"_MatchedLoosePhotons"]->Fill( et ); }
//       if ( iphoton->isLooseEM() )     { histos_[name+"_MatchedLooseEMObjects"]->Fill( et ); }

//       // Efficiency
//       histos_[name+"_EfficiencyEMObjects"]->Fill( et ); 
//       if ( iphoton->isTightPhoton() ) { histos_[name+"_EfficiencyTightPhotons"]->Fill( et ); }
//       if ( iphoton->isLoosePhoton() ) { histos_[name+"_EfficiencyLoosePhotons"]->Fill( et ); }
//       if ( iphoton->isLooseEM() )     { histos_[name+"_EfficiencyLooseEMObjects"]->Fill( et ); }
      
//       // PDG ids
//       //histos_[name+"_PdgId"]->Fill( daughter_id );
//       //histos_[name+"_MotherPdgId"]->Fill( mother_id );
      
//       // Kinematics
//       histos_[name+"_Energy"]->Fill( iphoton->energy() );
//       histos_[name+"_Et"]->Fill( iphoton->et() );
//       histos_[name+"_Eta"]->Fill( iphoton->eta() );
      
//       // Track isolation from RECO
//       histos_[name+"_NTrkSolid"]->Fill( iphoton->nTrkSolidCone() );
//       histos_[name+"_NTrkHollow"]->Fill( iphoton->nTrkHollowCone() );
//       histos_[name+"_IsoSolidTrk"]->Fill( iphoton->isolationSolidTrkCone() );
//       histos_[name+"_IsoHollowTrk"]->Fill( iphoton->isolationHollowTrkCone() );
      
//       // ECAL/HCAL isolation from RECO
//       histos_[name+"_IsoEcalRecHit"]->Fill( iphoton->isolationEcalRecHit() );
//       histos_[name+"_IsoHcalRecHit"]->Fill( iphoton->isolationHcalRecHit() );
//       histos_[name+"_HadronicOverEm"]->Fill( iphoton->hadronicOverEm() );
//       if ( fabs( iphoton->isolationEcalRecHit() ) > 0. ) {
// 	histos_[name+"_IsoHcalRecHitOverIsoEcalRecHit"]->Fill( iphoton->isolationHcalRecHit() / 
// 							       iphoton->isolationEcalRecHit() );
//       }
//       histos_[name+"_R9"]->Fill( iphoton->r9() );
      
//       // Isolation from PAT
//       histos_[name+"_IsoEcal"]->Fill( iphoton->ecalIso() );
//       histos_[name+"_IsoHcal"]->Fill( iphoton->hcalIso() );
//       histos_[name+"_IsoTrk"]->Fill( iphoton->trackIso() );
//       histos_[name+"_IsoAll"]->Fill( iphoton->caloIso() );
//       histos_[name+"_HcalIsoOverEcalIso"]->Fill( iphoton->hcalIso() / iphoton->ecalIso() );
      
//       // Photon ID flags
//       histos_[name+"_PhotonID"]->Fill( 0. );
//       if ( iphoton->isLooseEM() )     { histos_[name+"_PhotonID"]->Fill( 1. ); }
//       if ( iphoton->isLoosePhoton() ) { histos_[name+"_PhotonID"]->Fill( 2. ); }
//       if ( iphoton->isTightPhoton() ) { histos_[name+"_PhotonID"]->Fill( 3. ); }
//       if ( !iphoton->isLooseEM() &&
// 	   !iphoton->isLoosePhoton() &&
// 	   !iphoton->isTightPhoton() ) { histos_[name+"_PhotonID"]->Fill( 4. ); }
      
//       // Iterate through OtherPhotonID
//       if ( !others_.label().empty() ) {

// 	edm::View<pat::Photon>::const_iterator iother = others->begin(); 
// 	edm::View<pat::Photon>::const_iterator jother = others->end(); 
// 	for ( ; iother != jother; ++iother ) {

// 	  if ( uint32_t( iphoton - photons->begin() ) !=
// 	       uint32_t( iother - others->begin() ) ) { continue; }

// 	  histos2d_[name+"_CompareEcalIso"]->Fill( iphoton->isolationEcalRecHit(),
// 						   iother->isolationEcalRecHit() );
// 	  histos2d_[name+"_CompareHcalIso"]->Fill( iphoton->isolationHcalRecHit(),
// 						   iother->isolationHcalRecHit() );
// 	  histos2d_[name+"_CompareTrkIso"]->Fill( iphoton->isolationHollowTrkCone(),
// 						  iother->isolationHollowTrkCone() );
	  
// 	  float iii = 0.;
// 	  if ( iphoton->isTightPhoton() ) { iii = 3.; }
// 	  else if ( iphoton->isLoosePhoton() ) { iii = 2.; }
// 	  else if ( iphoton->isLooseEM() ) { iii = 1.; }
// 	  else { iii = 0.; }

// 	  float jjj = 0.;
// 	  if ( iother->isTightPhoton() ) { jjj = 3.; }
// 	  else if ( iother->isLoosePhoton() ) { jjj = 2.; }
// 	  else if ( iother->isLooseEM() ) { jjj = 1.; }
// 	  else { jjj = 0.; }

// 	  histos_[name+"_OtherPhotonID"]->Fill( jjj );
// 	  histos2d_[name+"_ComparePhotonID"]->Fill( iii, jjj );
	  
// 	}
//       }

//     } // photon loop

//     histos_[name+"_Multiplicity_Reco"]->Fill( multiplicity_reco[ii] );
//     histos_[name+"_Multiplicity_Truth"]->Fill( multiplicity_truth[ii] );
    
//   } // category loop
  
//   histos_["AllPhotons_Multiplicity_Reco"]->Fill( multiplicity_reco_all );
//   histos_["AllPhotons_Multiplicity_Truth"]->Fill( multiplicity_truth_all );
  
}

// -----------------------------------------------------------------------------
//
void SimplePhotonIDAnalysis::beginJob( const edm::EventSetup& ) {

  edm::Service<TFileService> fs;

  // Multiplicities  
  histos_["AllPhotons_Multiplicity_Reco"] = fs->make<TH1D>("AllPhotons_Multiplicity_Reco","",21,-0.5,20.5);
  histos_["AllPhotons_Multiplicity_Truth"] = fs->make<TH1D>("AllPhotons_Multiplicity_Truth","",21,-0.5,20.5);
  
  // DeltaR and DeltaPt
  histos_["DeltaR_Photons"] = fs->make<TH1D>("DeltaR_Photons","",100,0.,1.);
  histos_["DeltaR_Electrons"] = fs->make<TH1D>("DeltaR_Electrons","",100,0.,1.);
  histos_["DeltaR_Particles"] = fs->make<TH1D>("DeltaR_Particles","",100,0.,1.);
  histos_["DeltaPt_Photons"] = fs->make<TH1D>("DeltaPt_Photons","",100,0.,10.);
  histos_["DeltaPt_Electrons"] = fs->make<TH1D>("DeltaPt_Electrons","",100,0.,10.);
  histos_["DeltaPt_Particles"] = fs->make<TH1D>("DeltaPt_Particles","",100,0.,10.);
  histos2d_["DeltaRPt_Photons"] = fs->make<TH2D>("DeltaRPt_Photons","",25,0.,1.,25,0.,10.);
  histos2d_["DeltaRPt_Electrons"] = fs->make<TH2D>("DeltaRPt_Electrons","",25,0.,1.,25,0.,10.);
  histos2d_["DeltaRPt_Particles"] = fs->make<TH2D>("DeltaRPt_Particles","",25,0.,1.,25,0.,10.);
  
  for ( uint16_t ii = 0; ii < labels_.size(); ++ii ) {
    
    std::string name = labels_[ii]; 
    
    TFileDirectory dir = fs->mkdir(name);
    
    // PDG ids
    histos_[name+"_GenDaughterPdgId"] = dir.make<TH1D>("GenDaughterPdgId","",1000,0.,1000.);
    histos_[name+"_GenMotherPdgId"] = dir.make<TH1D>("GenMotherPdgId","",1000,0.,1000.);
    histos_[name+"_MatchedDaughterPdgId"] = dir.make<TH1D>("MatchedDaughterPdgId","",1000,0.,1000.);
    histos_[name+"_MatchedMotherPdgId"] = dir.make<TH1D>("MatchedMotherPdgId","",1000,0.,1000.);
    
    // Multiplicities
    histos_[name+"_Multiplicity_Reco"] = dir.make<TH1D>("Multiplicity_Reco","",21,-0.5,20.5);
    histos_[name+"_Multiplicity_Truth"] = dir.make<TH1D>("Multiplicity_Truth","",21,-0.5,20.5);
    
    // Totals 
    histos_[name+"_TotalEMObjects"] = dir.make<TH1D>("TotalEMObjects",";E_{T} [GeV];",nbins_,lower_,upper_);
    histos_[name+"_TotalLooseEMObjects"] = dir.make<TH1D>("TotalLooseEMObjects",";E_{T} [GeV];",nbins_,lower_,upper_);
    histos_[name+"_TotalLoosePhotons"] = dir.make<TH1D>("TotalLoosePhotons",";E_{T} [GeV];",nbins_,lower_,upper_);
    histos_[name+"_TotalTightPhotons"] = dir.make<TH1D>("TotalTightPhotons",";E_{T} [GeV];",nbins_,lower_,upper_);

    // Matched  
    histos_[name+"_MatchedEMObjects"] = dir.make<TH1D>("MatchedEMObjects",";E_{T} [GeV];",nbins_,lower_,upper_);
    histos_[name+"_MatchedLooseEMObjects"] = dir.make<TH1D>("MatchedLooseEMObjects",";E_{T} [GeV];",nbins_,lower_,upper_);
    histos_[name+"_MatchedLoosePhotons"] = dir.make<TH1D>("MatchedLoosePhotons",";E_{T} [GeV];",nbins_,lower_,upper_);
    histos_[name+"_MatchedTightPhotons"] = dir.make<TH1D>("MatchedTightPhotons",";E_{T} [GeV];",nbins_,lower_,upper_);

    // Efficiency
    histos_[name+"_EfficiencyEMObjects"] = dir.make<TH1D>("EfficiencyEMObjects",";E_{T} [GeV];",nbins_,lower_,upper_);
    histos_[name+"_EfficiencyLooseEMObjects"] = dir.make<TH1D>("EfficiencyLooseEMObjects",";E_{T} [GeV];",nbins_,lower_,upper_);
    histos_[name+"_EfficiencyLoosePhotons"] = dir.make<TH1D>("EfficiencyLoosePhotons",";E_{T} [GeV];",nbins_,lower_,upper_);
    histos_[name+"_EfficiencyTightPhotons"] = dir.make<TH1D>("EfficiencyTightPhotons",";E_{T} [GeV];",nbins_,lower_,upper_);

    // GenPhotons
    histos_[name+"_GenPhotons"] = dir.make<TH1D>("GenPhotons",";E_{T} [GeV];",nbins_,lower_,upper_);
    
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
//
void SimplePhotonIDAnalysis::endJob() {

  for ( uint16_t ii = 0; ii < labels_.size(); ++ii ) {
    
    std::string name = labels_[ii];

    if ( !matchedDaughterPdgId_[ii].empty() ) {
      TH1* th1 = histos_[name+"_MatchedDaughterPdgId"];
      if ( th1 ) {
	uint16_t size = matchedDaughterPdgId_[ii].size();
	size = size > 1000 ? 1000 : size;
	th1->SetBins( size, 0., size*1. );
	TAxis* axis = th1->GetXaxis();
	IdMap::const_iterator iii = matchedDaughterPdgId_[ii].begin();
	IdMap::const_iterator jjj = matchedDaughterPdgId_[ii].end();
	uint32_t bin = 0;
	float entries = 0;
	while ( bin < size && iii != jjj ) {
	  bin++;
	  entries += iii->second;
	  th1->SetBinContent( bin, static_cast<double>( iii->second ) );
	  std::stringstream ss; ss << iii->first;
	  if ( axis ) { axis->SetBinLabel( bin, ss.str().c_str() ); }
	  iii++;
	}
	th1->SetEntries( entries );
      }
    }
    
    if ( !matchedMotherPdgId_[ii].empty() ) {
      TH1* th1 = histos_[name+"_MatchedMotherPdgId"];
      if ( th1 ) {
	uint16_t size = matchedMotherPdgId_[ii].size();
	size = size > 1000 ? 1000 : size;
	th1->SetBins( size, 0., size*1. );
	TAxis* axis = th1->GetXaxis();
	IdMap::const_iterator iii = matchedMotherPdgId_[ii].begin();
	IdMap::const_iterator jjj = matchedMotherPdgId_[ii].end();
	uint32_t bin = 0;
	float entries = 0;
	while ( bin < size && iii != jjj ) {
	  bin++;
	  entries += iii->second;
	  th1->SetBinContent( bin, static_cast<double>( iii->second ) );
	  std::stringstream ss; ss << iii->first;
	  if ( axis ) { axis->SetBinLabel( bin, ss.str().c_str() ); }
	  iii++;
	}
	th1->SetEntries( entries );
      }
    }
    
    if ( !genDaughterPdgId_[ii].empty() ) {
      TH1* th1 = histos_[name+"_GenDaughterPdgId"];
      if ( th1 ) {
	uint16_t size = genDaughterPdgId_[ii].size();
	size = size > 1000 ? 1000 : size;
	th1->SetBins( size, 0., size*1. );
	TAxis* axis = th1->GetXaxis();
	IdMap::const_iterator iii = genDaughterPdgId_[ii].begin();
	IdMap::const_iterator jjj = genDaughterPdgId_[ii].end();
	uint32_t bin = 0;
	float entries = 0;
	while ( bin < size && iii != jjj ) {
	  bin++;
	  entries += iii->second;
	  th1->SetBinContent( bin, static_cast<double>( iii->second ) );
	  std::stringstream ss; ss << iii->first;
	  if ( axis ) { axis->SetBinLabel( bin, ss.str().c_str() ); }
	  iii++;
	}
	th1->SetEntries( entries );
      }
    }
    
    if ( !genMotherPdgId_[ii].empty() ) {
      TH1* th1 = histos_[name+"_GenMotherPdgId"];
      if ( th1 ) {
	uint16_t size = genMotherPdgId_[ii].size();
	size = size > 1000 ? 1000 : size;
	th1->SetBins( size, 0., size*1. );
	TAxis* axis = th1->GetXaxis();
	IdMap::const_iterator iii = genMotherPdgId_[ii].begin();
	IdMap::const_iterator jjj = genMotherPdgId_[ii].end();
	uint32_t bin = 0;
	float entries = 0;
	while ( bin < size && iii != jjj ) {
	  bin++;
	  entries += iii->second;
	  th1->SetBinContent( bin, static_cast<double>( iii->second ) );
	  std::stringstream ss; ss << iii->first;
	  if ( axis ) { axis->SetBinLabel( bin, ss.str().c_str() ); }
	  iii++;
	}
	th1->SetEntries( entries );
      }
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
      
//       float bin_width = nbins_ > 0 ? ( upper_ - lower_ ) / static_cast<float>( nbins_ ) : -1.;
//       int bin_number = bin_width > 0. ? static_cast<int>( ( iphoton->et() - lower_ ) / bin_width ) : -1;
//       if ( bin_number >= 0 && bin_number < nbins_ ) {
// 	for ( int bin = bin_number; bin < nbins_; ++bin  ) {
// 	  //for ( int bin = 0; bin < bin_number; ++bin  ) {
//  	  if ( labels_[ii] == name ) { histos_[name+"_AllEMObjects"]->Fill( bin*bin_width ); }
// 	  if ( iphoton->isLooseEM() ) { histos_[name+"_LooseEMObjects"]->Fill( bin*bin_width ); }
// 	  if ( iphoton->isLoosePhoton() ) { histos_[name+"_LoosePhotons"]->Fill( bin*bin_width ); }
// 	  if ( iphoton->isTightPhoton() ) { histos_[name+"_TightPhotons"]->Fill( bin*bin_width ); }
// 	}
//       }

