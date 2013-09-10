#include "TestCombination.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "CommonTools/Utils/interface/EtComparator.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "Combination.h"
#include "TH1D.h"
#include "TH2D.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cmath>

using namespace std;

// -----------------------------------------------------------------------------
//
TestCombination::TestCombination( const edm::ParameterSet& pset ) 
  : maximum_( pset.getParameter<int>("MaximumObjects") ),
    test_( pset.getParameter<int>("TestObjects") ),
    photons_( pset.getParameter<edm::InputTag>("Photons") ),
    jets_( pset.getParameter<edm::InputTag>("Jets") ),
    muons_( pset.getParameter<edm::InputTag>("Muons") ),
    electrons_( pset.getParameter<edm::InputTag>("Electrons") ),
    met_( pset.getParameter<edm::InputTag>("MET") ),
    ccMet_( pset.getParameter<edm::InputTag>("CCMET") ),
    genMet_( pset.getParameter<edm::InputTag>("GenMET") ),
    gen_(),// pset.getUntrackedParameter<edm::InputTag>("GenParticles") ),
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
    totalHt_( pset.getParameter<double>("TotalEt") ),
    nObjects_( pset.getUntrackedParameter<int>("NObjects",4) ),
    nPhotons_( pset.getUntrackedParameter<int>("NPhotons",2) ),
    minJets_(0),
    maxJets_(30),
    histos_(),
    histos2d_()
{
  int min = pset.getUntrackedParameter<int>("MinJets",-1);
  int max = pset.getUntrackedParameter<int>("MaxJets",-1);
  minJets_ = min < 0 ? 0 : (uint32_t)min;
  maxJets_ = max < 30 ? 0 : (uint32_t)max;
}

// -----------------------------------------------------------------------------
//
void TestCombination::analyze( const edm::Event& iEvent, 
			       const edm::EventSetup& iSetup ) {

  uint16_t cut_number = 0;
  histo("CutFlow_Efficiency")->Fill( cut_number ); cut_number++; //@@ increment!
  
  if ( test_ < 0 ) {

    // -------------------- Truth --------------------

    {
    
//       // Get GenParticles
//       edm::Handle< edm::View<reco::GenParticle> > gen;
//       iEvent.getByLabel( gen_, gen );
      
//       edm::View<reco::GenParticle>::const_iterator igen = gen->begin(); 
//       edm::View<reco::GenParticle>::const_iterator jgen = gen->end(); 
//       for ( ; igen != jgen; ++igen ) {
// 	const Candidate* mother = igen->mother();
// 	if ( mother && igen->pdgId() == igen->mother()->pdgId() ) { mother = mother->mother(); }
// 	if ( mother && mother->pdgId() == 1000022 ) { 
// 	  histos_["GenPhotons_Et"]->Fill( igen->et() ); 
// 	}
//       }
      
    }

    // -------------------- Retrieve objects and event selection --------------------

    // Get photons
    std::vector<Candidate> photons;
    if ( getPhotons( iEvent, photons ) ) { return; }

    // Get jets
    std::vector<Candidate> jets;
    if ( getJets( iEvent, jets ) ) { return; }
    
    // Get muons
    std::vector<Candidate> muons;
    if ( getMuons( iEvent, muons ) ) { return; }

    // Get electrons
    std::vector<Candidate> electrons;
    if ( getElectrons( iEvent, electrons ) ) { return; }

    // -------------------- Apply pre-selection --------------------

    histo("PreHtCut_NPhotons")->Fill( photons.size() ); 
    histo("PreHtCut_NJets")->Fill( jets.size() ); 
    histo("PreHtCut_NMuons")->Fill( muons.size() ); 
    histo("PreHtCut_NElectrons")->Fill( electrons.size() ); 

    if ( !muons.empty() ) { return; }
    histo("CutFlow_Efficiency")->Fill( cut_number ); cut_number++; //@@ increment!
    
    if ( !electrons.empty() ) { return; }
    histo("CutFlow_Efficiency")->Fill( cut_number ); cut_number++; //@@ increment!

    if ( jets.size() < minJets_ || jets.size() > maxJets_ ) { return; }
    histo("CutFlow_Efficiency")->Fill( cut_number ); cut_number++; //@@ increment!

    if ( photons.size() < nPhotons_ ) { return; }
    histo("CutFlow_Efficiency")->Fill( cut_number ); cut_number++; //@@ increment!

    // -------------------- "Common objects" --------------------
    
    std::vector<Candidate> objects;
    if ( !photons.empty() ) { objects.insert( objects.end(), photons.begin(), photons.end() ); }
    if ( !jets.empty() ) { objects.insert( objects.end(), jets.begin(), jets.end() ); }
    if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Objects: " << objects.size(); }

    { 
      double ht = HT( objects ); 
      histo("CommonObjectsHT")->Fill( ht ); 
      histo2d("CommonObjectsHT_Vs_NObjects")->Fill( objects.size(), ht ); 
    }
    
    // -------------------- Order by Et --------------------
    
    sort( objects.begin(), objects.end(), GreaterByEt<Candidate>() );
    sort( photons.begin(), photons.end(), GreaterByEt<Candidate>() );
    sort( jets.begin(), jets.end(), GreaterByEt<Candidate>() );

    // -------------------- Minimum number of objects --------------------

    if ( objects.size() < nObjects_ ) { return; }
    histo("CutFlow_Efficiency")->Fill( cut_number ); cut_number++; //@@ increment!

    // -------------------- HT cut --------------------

    {

      double ht = HT( objects );
      if ( ht < totalHt_ ) { return; }

      histo("CutFlow_Efficiency")->Fill( cut_number ); cut_number++; //@@ increment!
    
      histo("PostHtCut_NPhotons")->Fill( photons.size() ); 
      histo("PostHtCut_NJets")->Fill( jets.size() ); 
      histo("PostHtCut_NObjects")->Fill( objects.size() ); 

    }

    // -------------------- Get MET --------------------

    double primary_met = -1.;
    LorentzV lv_primary_met;
    std::vector<Candidate>::const_iterator iobj = objects.begin();
    std::vector<Candidate>::const_iterator jobj = objects.end();
    for ( ; iobj != jobj; ++iobj ) { lv_primary_met += iobj->p4(); }
    lv_primary_met.SetPx( -1.*lv_primary_met.Px() );
    lv_primary_met.SetPy( -1.*lv_primary_met.Py() );
    lv_primary_met.SetPz( -1.*lv_primary_met.Pz() );
    primary_met = sqrt( lv_primary_met.Perp2() );
    
    double calo_met = -1.;
    LorentzV lv_calo_met;
    if ( !met_.label().empty() ) {
      edm::Handle< std::vector<pat::MET> > handle;
      iEvent.getByLabel(met_,handle);
      if ( !handle.isValid() ) { 
	edm::LogWarning("TEST") << "No MET collection!";
	return;
      }
      if ( handle->empty() ) { 
	edm::LogWarning("TEST") << "Empty MET collection!";
	return;
      }
      if ( !handle->front().isCaloMET() ) {
	edm::LogWarning("TEST") << "Not CaloMET!" << endl;
      }
      calo_met = sqrt( handle->front().p4().Perp2() );
      lv_calo_met = handle->front().p4();
    }

    double cc_met = -1.;
    LorentzV lv_cc_met;
    if ( !ccMet_.label().empty() ) {
      edm::Handle< std::vector<pat::MET> > handle;
      iEvent.getByLabel(ccMet_,handle);
      if ( !handle.isValid() ) { 
	edm::LogWarning("TEST") << "No CCMET collection!";
	return;
      }
      if ( handle->empty() ) { 
	edm::LogWarning("TEST") << "Empty CCMET collection!";
	return;
      }
      if ( !handle->front().isCaloMET() ) {
	edm::LogWarning("TEST") << "Not CaloMET!" << endl;
      }
      cc_met = sqrt( handle->front().p4().Perp2() );
      lv_cc_met = handle->front().p4();
    }

    double gen_met = -1.;
    LorentzV lv_gen_met;
    if ( !genMet_.label().empty() ) {
      edm::Handle< std::vector<reco::GenMET> > handle;
      iEvent.getByLabel(genMet_,handle);
      if ( !handle.isValid() ) { 
	edm::LogWarning("TEST") << "No GenMET collection!";
	return;
      }
      if ( handle->empty() ) { 
	edm::LogWarning("TEST") << "Empty GenMET collection!";
	return;
      }
      gen_met = sqrt( handle->front().p4().Perp2() );
      lv_gen_met = handle->front().p4();
    }
    
    cout << " calo_met " << calo_met
	 << " cc_met " << cc_met
	 << " gen_met " << gen_met
	 << endl;

    // -------------------- MET -------------------- 

    {

      double ht = HT( objects ); 

      histo("PrimaryMET")->Fill( primary_met ); 
      histo("CaloMET")->Fill( calo_met ); 
      histo("CcMET")->Fill( cc_met ); 
      histo("GenMET")->Fill( gen_met ); 
      
      {
	double delta_phi = reco::deltaPhi<LorentzV,LorentzV>( lv_primary_met, lv_gen_met );
	histo("DPHI_PrimaryMET_GenMET")->Fill( delta_phi ); 
      }

      {
	double delta_phi = reco::deltaPhi<LorentzV,LorentzV>( lv_calo_met, lv_gen_met );
	histo("DPHI_CaloMET_GenMET")->Fill( delta_phi ); 
      }

      {
	double delta_phi = reco::deltaPhi<LorentzV,LorentzV>( lv_cc_met, lv_gen_met );
	histo("DPHI_CcMET_GenMET")->Fill( delta_phi ); 
      }
      
      histo2d("PrimaryMET_Vs_NObjects")->Fill( objects.size(), primary_met ); 
      histo2d("CaloMET_Vs_NObjects")->Fill( objects.size(), calo_met ); 
      histo2d("CcMET_Vs_NObjects")->Fill( objects.size(), cc_met ); 
      histo2d("GenMET_Vs_NObjects")->Fill( objects.size(), gen_met ); 

      histo2d("PrimaryMET_Vs_NJets")->Fill( jets.size(), primary_met ); 
      histo2d("CaloMET_Vs_NJets")->Fill( jets.size(), calo_met ); 
      histo2d("CcMET_Vs_NJets")->Fill( jets.size(), cc_met ); 
      histo2d("GenMET_Vs_NJets")->Fill( jets.size(), gen_met ); 
      
      histo2d("PrimaryMET_Vs_HT")->Fill( ht, primary_met ); 
      histo2d("CaloMET_Vs_HT")->Fill( ht, calo_met ); 
      histo2d("CcMET_Vs_HT")->Fill( ht, cc_met ); 
      histo2d("GenMET_Vs_HT")->Fill( ht, gen_met ); 
      
      histo2d("PrimaryMET_Vs_GenMET")->Fill( gen_met, primary_met ); 
      histo2d("CaloMET_Vs_GenMET")->Fill( gen_met, calo_met ); 
      histo2d("CcMET_Vs_GenMET")->Fill( gen_met, cc_met ); 
      histo2d("PrimaryMET_Vs_CaloMET")->Fill( calo_met, primary_met ); 
      
      if ( gen_met <= 0. ) { gen_met = 1.e-6; }
      if ( calo_met <= 0. ) { calo_met = 1.e-6; }

      histo2d("PrimaryMET/GenMET_Vs_NObjects")->Fill( objects.size(), primary_met/gen_met ); 
      histo2d("CaloMET/GenMET_Vs_NObjects")->Fill( objects.size(), calo_met/gen_met ); 
      histo2d("CcMET/GenMET_Vs_NObjects")->Fill( objects.size(), cc_met/gen_met ); 
      histo2d("PrimaryMET/CaloMET_Vs_NObjects")->Fill( objects.size(), primary_met/calo_met ); 
      
      histo2d("PrimaryMET/GenMET_Vs_NJets")->Fill( jets.size(), primary_met/gen_met ); 
      histo2d("CaloMET/GenMET_Vs_NJets")->Fill( jets.size(), calo_met/gen_met ); 
      histo2d("CcMET/GenMET_Vs_NJets")->Fill( jets.size(), cc_met/gen_met ); 
      histo2d("PrimaryMET/CaloMET_Vs_NJets")->Fill( jets.size(), primary_met/calo_met ); 
      
      histo2d("PrimaryMET/GenMET_Vs_HT")->Fill( ht, primary_met/gen_met ); 
      histo2d("CaloMET/GenMET_Vs_HT")->Fill( ht, calo_met/gen_met ); 
      histo2d("CcMET/GenMET_Vs_HT")->Fill( ht, cc_met/gen_met ); 
      histo2d("PrimaryMET/CaloMET_Vs_HT")->Fill( ht, primary_met/calo_met ); 
      
      histo("PrimaryMET-GenMET/GenMET")->Fill( (primary_met-gen_met)/gen_met ); 
      histo("CaloMET-GenMET/GenMET")->Fill( (calo_met-gen_met)/gen_met ); 
      histo("CcMET-GenMET/GenMET")->Fill( (cc_met-gen_met)/gen_met ); 
      histo("PrimaryMET-CaloMET/CaloMET")->Fill( (primary_met-calo_met)/calo_met ); 
      
      histo2d("PrimaryMET-GenMET/GenMET_Vs_NObjects")->Fill( objects.size(), (primary_met-gen_met)/gen_met ); 
      histo2d("CaloMET-GenMET/GenMET_Vs_NObjects")->Fill( objects.size(), (calo_met-gen_met)/gen_met ); 
      histo2d("CcMET-GenMET/GenMET_Vs_NObjects")->Fill( objects.size(), (cc_met-gen_met)/gen_met ); 
      histo2d("PrimaryMET-CaloMET/CaloMET_Vs_NObjects")->Fill( objects.size(), (primary_met-calo_met)/calo_met ); 

      histo2d("PrimaryMET-GenMET/GenMET_Vs_HT")->Fill( ht, (primary_met-gen_met)/gen_met ); 
      histo2d("CaloMET-GenMET/GenMET_Vs_HT")->Fill( ht, (calo_met-gen_met)/gen_met ); 
      histo2d("CcMET-GenMET/GenMET_Vs_HT")->Fill( ht, (cc_met-gen_met)/gen_met ); 
      histo2d("PrimaryMET-CaloMET/CaloMET_Vs_HT")->Fill( ht, (primary_met-calo_met)/calo_met ); 

      histo2d("PrimaryMET-GenMET/GenMET_Vs_NObjects_Vs_HT")->Fill( objects.size(), ht, (primary_met-gen_met)/gen_met ); 
      histo2d("CaloMET-GenMET/GenMET_Vs_NObjects_Vs_HT")->Fill( objects.size(), ht, (calo_met-gen_met)/gen_met ); 
      histo2d("CcMET-GenMET/GenMET_Vs_NObjects_Vs_HT")->Fill( objects.size(), ht, (cc_met-gen_met)/gen_met ); 
      histo2d("PrimaryMET-CaloMET/CaloMET_Vs_NObjects_Vs_HT")->Fill( objects.size(), ht, (primary_met-calo_met)/calo_met ); 
      
    }

    // -------------------- MHT vs HT --------------------

    {
      
      double ht_pho = HT( photons ); 
      double ht_jet = HT( jets ); 
      double ht_all = HT( objects ); 
      double mht = MHT( objects );
      
      histo2d("MHT_Vs_PhotonHT")->Fill( ht_pho, mht ); 
      histo2d("MHT_Vs_JetHT")->Fill( ht_jet, mht ); 
      histo2d("MHT_Vs_HT")->Fill( ht_all, mht ); 

    }

    // -------------------- alphaT and betaT --------------------
    
    {
      
      double ht  = HT( objects );
      double mht = MHT( objects );
      double mt  = MT( objects );
      double ht_mht = HT_MHT( objects );
      
      std::vector<Candidate> jet1;
      std::vector<Candidate> jet2;
      double min_dht = minDHT( objects, jet1, jet2 );
      
      histo2d("MHT/HT_Vs_DHT/HT")->Fill( min_dht / ht, mht / ht ); 
      
      histo("HT")->Fill( ht ); 
      histo2d("HT_Vs_NObjects")->Fill( objects.size(), ht ); 

      histo("MT")->Fill( mt ); 
      histo2d("MT_Vs_NObjects")->Fill( objects.size(), mt ); 
      
      histo("MHT")->Fill( mht ); 
      histo2d("MHT_Vs_NObjects")->Fill( objects.size(), mht ); 
      
      histo("MinDHT")->Fill( min_dht ); 
      histo2d("MinDHT_Vs_NObjects")->Fill( objects.size(), min_dht ); 
      histo2d("MinDHT_Vs_MHT")->Fill( mht, min_dht );
      
      histo("HT-MHT")->Fill( ht_mht ); 
      histo2d("HT-MHT_Vs_NObjects")->Fill( objects.size(), ht_mht ); 
      histo2d("HT-MHT_Vs_MT")->Fill( mt, ht_mht );

      double alpha_t = alphaT( min_dht, ht, mt );

      histo("AlphaT")->Fill( alpha_t ); 
      histo2d("AlphaT_Vs_NObjects")->Fill( objects.size(), alpha_t ); 

      double beta_t = alphaT( min_dht, ht, ht_mht );

      histo("BetaT")->Fill( beta_t ); 
      histo2d("BetaT_Vs_NObjects")->Fill( objects.size(), beta_t ); 

      histo2d("BetaT_Vs_AlphaT")->Fill( alpha_t, beta_t ); 

    }

    // -------------------- Biased alphaT and betaT --------------------

    {
      
      double ht  = HT( jets ); //@@ no photons!
      double mht = MHT( objects );
      double mt  = ht*ht - mht*mht; //@@ no photons!
      mt = mt < 0. ? -1.*sqrt(-1.*mt) : sqrt(mt);
      double ht_mht = ht - mht;

      std::vector<Candidate> jet1;
      std::vector<Candidate> jet2;
      double min_dht = minDHT( objects, jet1, jet2 );

      histo("BiasedHT")->Fill( ht ); 
      histo2d("BiasedHT_Vs_NObjects")->Fill( objects.size(), ht ); 

      histo("BiasedMT")->Fill( mt ); 
      histo2d("BiasedMT_Vs_NObjects")->Fill( objects.size(), mt ); 
      
      histo("BiasedMHT")->Fill( mht ); 
      histo2d("BiasedMHT_Vs_NObjects")->Fill( objects.size(), mht ); 
      
      histo("BiasedMinDHT")->Fill( min_dht ); 
      histo2d("BiasedMinDHT_Vs_NObjects")->Fill( objects.size(), min_dht ); 
      histo2d("BiasedMinDHT_Vs_BiasedMHT")->Fill( mht, min_dht );
      
      histo("BiasedHT-BiasedMHT")->Fill( ht_mht ); 
      histo2d("BiasedHT-BiasedMHT_Vs_NObjects")->Fill( objects.size(), ht_mht ); 
      histo2d("BiasedHT-BiasedMHT_Vs_BiasedMT")->Fill( mt, ht_mht );
      
      double alpha_t = alphaT( min_dht, ht, mt );
      
      histo("BiasedAlphaT")->Fill( alpha_t ); 
      histo2d("BiasedAlphaT_Vs_NJets")->Fill( jets.size(), alpha_t ); 
      histo2d("BiasedAlphaT_Vs_AlphaT")->Fill( alpha_t, alpha_t ); 
      
      double beta_t = alphaT( min_dht, ht, ht_mht );
      
      histo("BiasedBetaT")->Fill( beta_t ); 
      histo2d("BiasedBetaT_Vs_NJets")->Fill( jets.size(), beta_t ); 
      histo2d("BiasedBetaT_Vs_BetaT")->Fill( beta_t, beta_t ); 

    }

    // -------------------- Biased DPHI --------------------

    {

      LorentzV lv_obj;
      std::vector<Candidate>::const_iterator iobj = objects.begin();
      std::vector<Candidate>::const_iterator jobj = objects.end();
      for ( ; iobj != jobj; ++iobj ) { lv_obj += iobj->p4(); }
      
      LorentzV lv_recoil( lv_obj );
      lv_recoil.SetPx( -1.*lv_recoil.Px() );
      lv_recoil.SetPy( -1.*lv_recoil.Py() );
      lv_recoil.SetPz( -1.*lv_recoil.Pz() );

      // Min DeltaPhi
      {

	double min_delta_phi = 10.;
	iobj = objects.begin();
	jobj = objects.end();
	for ( ; iobj != jobj; ++iobj ) { 
	  double delta_phi = reco::deltaPhi( iobj->p4().phi(), lv_recoil.phi() );
	  if ( fabs(delta_phi) < fabs(min_delta_phi) ) { min_delta_phi = delta_phi; }
	}
	
	histo("DPHI")->Fill( fabs(min_delta_phi) ); 
	histo2d("DPHI_Vs_NObjects")->Fill( objects.size(), min_delta_phi ); 

      }

      // Min BiasedDeltaPhi
      {
	
	double min_delta_phi = 10.;
	iobj = objects.begin();
	jobj = objects.end();
	for ( ; iobj != jobj; ++iobj ) { 
	  LorentzV biased_recoil = lv_recoil + iobj->p4();
	  double delta_phi = reco::deltaPhi<LorentzV,LorentzV>( biased_recoil, lv_recoil );
	  if ( fabs(delta_phi) < fabs(min_delta_phi) ) { min_delta_phi = delta_phi; }
	}
	
	histo("BiasedDPHI")->Fill( fabs(min_delta_phi) ); 
	histo2d("BiasedDPHI_Vs_NObjects")->Fill( objects.size(), min_delta_phi ); 

      }
      
    }

    // -------------------- Recoil --------------------

    if ( photons.size() >= 2 ) {
      
      double ht  = HT( objects );
      double mht = MHT( objects );
      double mt  = ht*ht - mht*mht; 
      mt = mt < 0. ? -1.*sqrt(-1.*mt) : sqrt(mt);
      double ht_mht = ht - mht;
      std::vector<Candidate> jet1;
      std::vector<Candidate> jet2;
      double dht = recoilDHT( photons, jets, jet1, jet2 );
	
      histo2d("DHT_Vs_NJets")->Fill( jets.size(), dht ); 
      histo2d("MHT_Vs_NJets")->Fill( jets.size(), mht ); 
	
      histo2d("DHT_Vs_MHT")->Fill( mht, dht ); 
      histo2d("DHT_Vs_MHT/HT")->Fill( mht / ht, dht ); 
      histo2d("DHT/MHT_Vs_HT")->Fill( ht, dht / mht ); 
      //double delta_phi = reco::deltaPhi<LorentzV,LorentzV>(lv_jet,lv_recoil);
      //histo2d("DPHI_Vs_MHT/HT")->Fill( mht / ht, delta_phi ); 
	
      double recoil_alpha_t = alphaT( fabs(dht), ht, mt );
	
      histo("RecoilAlphaT")->Fill( recoil_alpha_t ); 
      histo2d("RecoilAlphaT_Vs_NJets")->Fill( jets.size(), recoil_alpha_t ); 
	
      double recoil_beta_t = alphaT( fabs(dht), ht, ht_mht );

      histo("RecoilBetaT")->Fill( recoil_beta_t ); 
      histo2d("RecoilBetaT_Vs_NJets")->Fill( jets.size(), recoil_beta_t ); 
	
    }      

    // -------------------- Projected --------------------

    if ( photons.size() >= 2 ) {
      
      // Photon system
      LorentzV lv_pho;
      std::vector<Candidate>::const_iterator ipho = photons.begin();
      std::vector<Candidate>::const_iterator jpho = photons.begin() + 2; 
      for ( ; ipho != jpho; ++ipho ) { lv_pho += ipho->p4(); }
      
      // "Recoil" axis
      LorentzV lv_recoil( lv_pho );
      lv_recoil.SetPx( -1.*lv_recoil.Px() );
      lv_recoil.SetPy( -1.*lv_recoil.Py() );
      lv_recoil.SetPz( -1.*lv_recoil.Pz() );
      
      // MHT from complete system
      LorentzV lv_obj;
      std::vector<Candidate>::const_iterator iobj = objects.begin();
      std::vector<Candidate>::const_iterator jobj = objects.end();
      for ( ; iobj != jobj; ++iobj ) { lv_obj += iobj->p4(); }
      LorentzV lv_mht( lv_obj );
      lv_mht.SetPx( -1.*lv_obj.Px() );
      lv_mht.SetPy( -1.*lv_obj.Py() );
      lv_mht.SetPz( -1.*lv_obj.Pz() );
      
      // Projection of MHT onto recoil axis
      LorentzV::Scalar dot_product_recoil = lv_recoil.Dot( lv_recoil );
      LorentzV::Scalar dot_product_mht = lv_mht.Dot( lv_recoil );
      LorentzV lv_proj = lv_recoil * ( dot_product_mht / dot_product_recoil );

      double ht  = lv_obj.Et();
      double mht = MHT( objects );
      double mt  = ht*ht - mht*mht;
      mt = mt < 0. ? -1.*sqrt(-1.*mt) : sqrt(mt);
      double ht_mht = ht - mht;
      double dht = sqrt( lv_proj.Perp2() );// - sqrt( lv_mht.Perp2() );
	
      histo2d("ProjectedDHT_Vs_NJets")->Fill( jets.size(), dht ); 
	
      histo2d("ProjectedDHT_Vs_MHT")->Fill( mht, dht ); 
      histo2d("ProjectedDHT_Vs_MHT/HT")->Fill( mht / ht, dht ); 
      histo2d("ProjectedDHT/MHT_Vs_HT")->Fill( ht, dht / mht ); 
	
      double new_alpha_t = alphaT( fabs(dht), ht, mt );
	
      histo("ProjectedAlphaT")->Fill( new_alpha_t ); 
      histo2d("ProjectedAlphaT_Vs_NJets")->Fill( jets.size(), new_alpha_t ); 
	
      double new_beta_t = alphaT( fabs(dht), ht, ht_mht );

      histo("ProjectedBetaT")->Fill( new_beta_t ); 
      histo2d("ProjectedBetaT_Vs_NJets")->Fill( jets.size(), new_beta_t ); 

    }
    
    // -------------------- Super --------------------
    
    {
      
      double ht  = HT( objects );
      double mht = MHT( objects );
      double mt  = ht*ht - mht*mht; 
      mt = mt < 0. ? -1.*sqrt(-1.*mt) : sqrt(mt);
      double ht_mht = ht - mht;
      std::vector<Candidate> jet1;
      std::vector<Candidate> jet2;
      double super_min_dht = minSuperDHT( objects, jet1, jet2 );
      
      double super_alpha_t = alphaT( super_min_dht, ht, mt );
      
      histo("SuperMinDHT")->Fill( super_min_dht ); 
      histo2d("SuperMinDHT_Vs_NObjects")->Fill( objects.size(), super_min_dht ); 
      
      histo("SuperAlphaT")->Fill( super_alpha_t ); 
      histo2d("SuperAlphaT_Vs_NObjects")->Fill( objects.size(), super_alpha_t ); 

      edm::LogVerbatim("TEST")
	<< " Number of objects = " << objects.size()
	<< ", in jet #1 = " << jet1.size()
	<< ", in jet #2 = " << jet2.size()
	<< std::endl;

      uint16_t cntr = 0;
      std::vector<Candidate>::const_iterator ii = jet1.begin();
      std::vector<Candidate>::const_iterator jj = jet1.end();
      for ( ; ii != jj; ++ii ) { if ( ii->pdgId() == 22 ) { cntr++; } }
      
      histo("SuperAlphaT_NPhotons")->Fill( cntr ); 
      histo2d("SuperAlphaT_NPhotons_Vs_NObjects")->Fill( objects.size(), cntr ); 

      double super_beta_t = alphaT( super_min_dht, ht, ht_mht );
      
      histo("SuperBetaT")->Fill( super_beta_t ); 
      histo2d("SuperBetaT_Vs_NObjects")->Fill( objects.size(), super_beta_t ); 
      
      histo2d("SuperBetaT_Vs_SuperAlphaT")->Fill( super_alpha_t, super_beta_t ); 
      
    }

    // -------------------- Photons in pseudo-jets --------------------

    {    
      
      std::vector<Candidate> jet1;
      std::vector<Candidate> jet2;
      double min_dht = minDHT( objects, jet1, jet2 );
      
      std::vector<int> phot1;
      std::vector<Candidate>::const_iterator ii1 = jet1.begin();
      std::vector<Candidate>::const_iterator jj1 = jet1.end();
      for ( ; ii1 != jj1; ++ii1 ) { if ( ii1->pdgId() == 22 ) { phot1.push_back( ii1 - jet1.begin() ); } }
    
      std::vector<int> phot2;
      std::vector<Candidate>::const_iterator ii2 = jet2.begin();
      std::vector<Candidate>::const_iterator jj2 = jet2.end();
      for ( ; ii2 != jj2; ++ii2 ) { if ( ii2->pdgId() == 22 ) { phot2.push_back( ii2 - jet2.begin() ); } }
  
      histo2d("PseudoDijets_NObjects")->Fill( jet1.size(), jet2.size() ); 
      histo2d("PseudoDijets_NPhotons")->Fill( phot1.size(), phot2.size() ); 
      histo2d("PseudoDijets_NJets")->Fill( jet1.size() - phot1.size(), jet2.size() - phot2.size() ); 
    
      {
	std::vector<int>::const_iterator iii1 = phot1.begin();
	std::vector<int>::const_iterator jjj1 = phot1.end();
	for ( ; iii1 != jjj1; ++iii1 ) {
	  std::vector<int>::const_iterator iii2 = phot2.begin();
	  std::vector<int>::const_iterator jjj2 = phot2.end();
	  for ( ; iii2 != jjj2; ++iii2 ) {
	    histo2d("PseudoDijets_PhotonEt")->Fill( jet1[*iii1].et(), jet2[*iii2].et() ); 
	    double photon_et = jet1[*iii1].et() < jet2[*iii2].et() ? jet1[*iii1].et() : jet2[*iii2].et();
	    double delta_phi = reco::deltaPhi<Candidate,Candidate>( jet1[*iii1], jet2[*iii2] );
	    histo2d("PseudoDijets_DeltaPhi_Vs_PhotonEt")->Fill( photon_et, delta_phi ); 
	    histo2d("PseudoDijets_DeltaPhi_Vs_NObjects")->Fill( objects.size(), delta_phi ); 
	  }
	}
      }

      // -------------------- Debug -------------------- 
      
      if ( edm::isDebugEnabled() ) {
	std::stringstream ss;
	ss << "Minimum deltaEt = " << min_dht << std::endl;

	{
	  ss << "Input objects: id/pdg: ";
	  std::vector<Candidate>::const_iterator ii = objects.begin();
	  std::vector<Candidate>::const_iterator jj = objects.end();
	  for ( ; ii != jj; ++ii  ) { 
	    ss << ii->pdgId() << "/" << ii->pt() << ", "; 
	  }
	  ss << std::endl;
	}

	{
	  ss << "Pseudo-jet #1: id/pdg: ";
	  std::vector<Candidate>::const_iterator ii = jet1.begin();
	  std::vector<Candidate>::const_iterator jj = jet1.end();
	  for ( ; ii != jj; ++ii  ) { 
	    ss << ii->pdgId() << "/" << ii->pt() << ", "; 
	  }
	  ss << std::endl;
	}

	{
	  ss << "Pseudo-jet #2: id/pdg: ";
	  std::vector<Candidate>::const_iterator ii = jet2.begin();
	  std::vector<Candidate>::const_iterator jj = jet2.end();
	  for ( ; ii != jj; ++ii  ) { 
	    ss << ii->pdgId() << "/" << ii->pt() << ", "; 
	  }
	  ss << std::endl;
	}
    
	LogTrace("TEST") << ss.str();

      }

    }

  } else { // -------------------- Test --------------------

    std::vector<Candidate> input;
    for ( uint16_t ii = 1; ii <= (iEvent.id().event()-1)%test_; ++ii ) {
      input.push_back( Candidate( 0, LorentzV( 1./ii, 0, 0, ii ) ) );
    }
    std::vector<Candidate> jet1;
    std::vector<Candidate> jet2;
    double min_dht = minDHT( input, jet1, jet2 );
    
    // -------------------- Debug --------------------

    if ( edm::isDebugEnabled() ) {
      std::stringstream ss;
      ss << "Minimum deltaEt = " << min_dht << std::endl;

      {
	ss << "Input objects: pdg/pt: ";
	std::vector<Candidate>::const_iterator ii = input.begin();
	std::vector<Candidate>::const_iterator jj = input.end();
	for ( ; ii != jj; ++ii  ) { 
	  ss << ii->pdgId() << "/" << ii->pt() << ", "; 
	}
	ss << std::endl;
      }

      {
	ss << "Pseudo-jet #1: pdg/pt: ";
	std::vector<Candidate>::const_iterator ii = jet1.begin();
	std::vector<Candidate>::const_iterator jj = jet1.end();
	for ( ; ii != jj; ++ii  ) { 
	  ss << ii->pdgId() << "/" << ii->pt() << ", "; 
	}
	ss << std::endl;
      }

      {
	ss << "Pseudo-jet #2: pdg/pt: ";
	std::vector<Candidate>::const_iterator ii = jet2.begin();
	std::vector<Candidate>::const_iterator jj = jet2.end();
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

    // -------------------- "MET" --------------------
  
    { 
      
      TFileDirectory dir = fs->mkdir("MET");
      
      histos_["PrimaryMET"] = dir.make<TH1D>("PrimaryMET","",100,0.,2000.);
      histos_["CaloMET"] = dir.make<TH1D>("CaloMET","",100,0.,2000.);
      histos_["CcMET"] = dir.make<TH1D>("CcMET","",100,0.,2000.);
      histos_["GenMET"] = dir.make<TH1D>("GenMET","",100,0.,2000.);

      histos_["DPHI_PrimaryMET_GenMET"] = dir.make<TH1D>("DPHI_PrimaryMET_GenMET","",160,-3.2,3.2);
      histos_["DPHI_CaloMET_GenMET"] = dir.make<TH1D>("DPHI_CaloMET_GenMET","",160,-3.2,3.2);
      histos_["DPHI_CcMET_GenMET"] = dir.make<TH1D>("DPHI_CcMET_GenMET","",160,-3.2,3.2);

      histos2d_["PrimaryMET_Vs_NObjects"] = dir.make<TH2D>("PrimaryMET_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);
      histos2d_["CaloMET_Vs_NObjects"] = dir.make<TH2D>("CaloMET_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);
      histos2d_["CcMET_Vs_NObjects"] = dir.make<TH2D>("CcMET_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);
      histos2d_["GenMET_Vs_NObjects"] = dir.make<TH2D>("GenMET_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);

      histos2d_["PrimaryMET_Vs_NJets"] = dir.make<TH2D>("PrimaryMET_Vs_NJets","",51,-0.5,50.5,100,0.,2000.);
      histos2d_["CaloMET_Vs_NJets"] = dir.make<TH2D>("CaloMET_Vs_NJets","",51,-0.5,50.5,100,0.,2000.);
      histos2d_["CcMET_Vs_NJets"] = dir.make<TH2D>("CcMET_Vs_NJets","",51,-0.5,50.5,100,0.,2000.);
      histos2d_["GenMET_Vs_NJets"] = dir.make<TH2D>("GenMET_Vs_NJets","",51,-0.5,50.5,100,0.,2000.);
      
      histos2d_["PrimaryMET_Vs_HT"] = dir.make<TH2D>("PrimaryMET_Vs_HT","",100,0.,2000.,100,0.,2000.);
      histos2d_["CaloMET_Vs_HT"] = dir.make<TH2D>("CaloMET_Vs_HT","",100,0.,2000.,100,0.,2000.);
      histos2d_["CcMET_Vs_HT"] = dir.make<TH2D>("CcMET_Vs_HT","",100,0.,2000.,100,0.,2000.);
      histos2d_["GenMET_Vs_HT"] = dir.make<TH2D>("GenMET_Vs_HT","",100,0.,2000.,100,0.,2000.);
      
      histos2d_["PrimaryMET_Vs_GenMET"] = dir.make<TH2D>("PrimaryMET_Vs_GenMET","",100,0.,2000.,100,0.,2000.);
      histos2d_["CaloMET_Vs_GenMET"] = dir.make<TH2D>("CaloMET_Vs_GenMET","",100,0.,2000.,100,0.,2000.);
      histos2d_["CcMET_Vs_GenMET"] = dir.make<TH2D>("CcMET_Vs_GenMET","",100,0.,2000.,100,0.,2000.);
      histos2d_["PrimaryMET_Vs_CaloMET"] = dir.make<TH2D>("PrimaryMET_Vs_CaloMET","",100,0.,2000.,100,0.,2000.);
      
      histos2d_["PrimaryMET/GenMET_Vs_NObjects"] = dir.make<TH2D>("PrimaryMET/GenMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);
      histos2d_["CaloMET/GenMET_Vs_NObjects"] = dir.make<TH2D>("CaloMET/GenMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);
      histos2d_["CcMET/GenMET_Vs_NObjects"] = dir.make<TH2D>("CcMET/GenMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);
      histos2d_["PrimaryMET/CaloMET_Vs_NObjects"] = dir.make<TH2D>("PrimaryMET/CaloMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);

      histos2d_["PrimaryMET/GenMET_Vs_NJets"] = dir.make<TH2D>("PrimaryMET/GenMET_Vs_NJets","",51,-0.5,50.5,200,-10.,10.);
      histos2d_["CaloMET/GenMET_Vs_NJets"] = dir.make<TH2D>("CaloMET/GenMET_Vs_NJets","",51,-0.5,50.5,200,-10.,10.);
      histos2d_["CcMET/GenMET_Vs_NJets"] = dir.make<TH2D>("CcMET/GenMET_Vs_NJets","",51,-0.5,50.5,200,-10.,10.);
      histos2d_["PrimaryMET/CaloMET_Vs_NJets"] = dir.make<TH2D>("PrimaryMET/CaloMET_Vs_NJets","",51,-0.5,50.5,200,-10.,10.);
      
      histos2d_["PrimaryMET/GenMET_Vs_HT"] = dir.make<TH2D>("PrimaryMET/GenMET_Vs_HT","",100,0.,2000.,200,-10.,10.);
      histos2d_["CaloMET/GenMET_Vs_HT"] = dir.make<TH2D>("CaloMET/GenMET_Vs_HT","",100,0.,2000.,200,-10.,10.);
      histos2d_["CcMET/GenMET_Vs_HT"] = dir.make<TH2D>("CcMET/GenMET_Vs_HT","",100,0.,2000.,200,-10.,10.);
      histos2d_["PrimaryMET/CaloMET_Vs_HT"] = dir.make<TH2D>("PrimaryMET/CaloMET_Vs_HT","",100,0.,2000.,200,-10.,10.);

      histos_["PrimaryMET-GenMET/GenMET"] = dir.make<TH1D>("PrimaryMET-GenMET/GenMET","",51,-0.5,50.5);
      histos_["CaloMET-GenMET/GenMET"] = dir.make<TH1D>("CaloMET-GenMET/GenMET","",51,-0.5,50.5);
      histos_["CcMET-GenMET/GenMET"] = dir.make<TH1D>("CcMET-GenMET/GenMET","",51,-0.5,50.5);
      histos_["PrimaryMET-CaloMET/CaloMET"] = dir.make<TH1D>("PrimaryMET-CaloMET/CaloMET","",51,-0.5,50.5);

      histos2d_["PrimaryMET-GenMET/GenMET_Vs_NObjects"] = 
	dir.make<TH2D>("PrimaryMET-GenMET/GenMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);
      histos2d_["CaloMET-GenMET/GenMET_Vs_NObjects"] = 
	dir.make<TH2D>("CaloMET-GenMET/GenMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);
      histos2d_["CcMET-GenMET/GenMET_Vs_NObjects"] = 
	dir.make<TH2D>("CcMET-GenMET/GenMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);
      histos2d_["PrimaryMET-CaloMET/CaloMET_Vs_NObjects"] = 
	dir.make<TH2D>("PrimaryMET-CaloMET/CaloMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);
      
      histos2d_["PrimaryMET-GenMET/GenMET_Vs_HT"] = 
	dir.make<TH2D>("PrimaryMET-GenMET/GenMET_Vs_HT","",100,0.,2000.,200,-10.,10.);
      histos2d_["CaloMET-GenMET/GenMET_Vs_HT"] = 
	dir.make<TH2D>("CaloMET-GenMET/GenMET_Vs_HT","",100,0.,2000.,200,-10.,10.);
      histos2d_["CcMET-GenMET/GenMET_Vs_HT"] = 
	dir.make<TH2D>("CcMET-GenMET/GenMET_Vs_HT","",100,0.,2000.,200,-10.,10.);
      histos2d_["PrimaryMET-CaloMET/CaloMET_Vs_HT"] = 
	dir.make<TH2D>("PrimaryMET-CaloMET/CaloMET_Vs_HT","",100,0.,2000.,200,-10.,10.);

      histos2d_["PrimaryMET-GenMET/GenMET_Vs_NObjects_Vs_HT"] = 
	dir.make<TH2D>("PrimaryMET-GenMET/GenMET_Vs_NObjects_Vs_HT","",51,-0.5,50.5,100,0.,2000.);
      histos2d_["CaloMET-GenMET/GenMET_Vs_NObjects_Vs_HT"] = 
	dir.make<TH2D>("CaloMET-GenMET/GenMET_Vs_NObjects_Vs_HT","",51,-0.5,50.5,100,0.,2000.);
      histos2d_["CcMET-GenMET/GenMET_Vs_NObjects_Vs_HT"] = 
	dir.make<TH2D>("CcMET-GenMET/GenMET_Vs_NObjects_Vs_HT","",51,-0.5,50.5,100,0.,2000.);
      histos2d_["PrimaryMET-CaloMET/CaloMET_Vs_NObjects_Vs_HT"] = 
	dir.make<TH2D>("PrimaryMET-CaloMET/CaloMET_Vs_NObjects_Vs_HT","",51,-0.5,50.5,100,0.,2000.);

    }      

    // -------------------- "Common objects" --------------------
  
    { 
      
      TFileDirectory dir = fs->mkdir("CommonObjects");
      
      histos_["PreHtCut_NPhotons"] = dir.make<TH1D>("PreHtCut_NPhotons","",51,-0.5,50.5);
      histos_["PreHtCut_NJets"] = dir.make<TH1D>("PreHtCut_NJets","",51,-0.5,50.5);
      histos_["PreHtCut_NMuons"] = dir.make<TH1D>("PreHtCut_NMuons","",51,-0.5,50.5);
      histos_["PreHtCut_NElectrons"] = dir.make<TH1D>("PreHtCut_NElectrons","",51,-0.5,50.5);
      
      histos_["PostHtCut_NPhotons"] = dir.make<TH1D>("PostHtCut_NPhotons","",51,-0.5,50.5);
      histos_["PostHtCut_NJets"] = dir.make<TH1D>("PostHtCut_NJets","",51,-0.5,50.5);
      histos_["PostHtCut_NObjects"] = dir.make<TH1D>("PostHtCut_NObjects","",51,-0.5,50.5);
      
      histos_["GenPhotons_Et"] = dir.make<TH1D>("GenPhotons_Et","",200,0.,1000.);
      
    }

    // -------------------- Event weight --------------------

    { 

      TFileDirectory dir = fs->mkdir("Common");

      histos_["CommonObjectsHT"] = dir.make<TH1D>("CommonObjectsHT","",100,0.,2000.);
      histos2d_["CommonObjectsHT_Vs_NObjects"] = dir.make<TH2D>("CommonObjectsHT_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);

      histos_["CutFlow_Efficiency"] = dir.make<TH1D>("CutFlow_Efficiency","",11,-0.5,10.5);
      //histos_["EventWeight_Coefficient"] = dir.make<TH1D>("EventWeight_Coefficient","",1000,0.,10.);
      //histos_["EventWeight_Base"] = dir.make<TH1D>("EventWeight_Base","",21,-10.5,10.5);
    
    }

    // -------------------- MHT --------------------

    { 
      
      TFileDirectory dir = fs->mkdir("MHT");
      
      histos2d_["MHT_Vs_PhotonHT"] = dir.make<TH2D>("MHT_Vs_PhotonHT","",100,0.,2000.,100,0.,1000.);
      histos2d_["MHT_Vs_JetHT"] = dir.make<TH2D>("MHT_Vs_JetHT","",100,0.,2000.,100,0.,1000.);
      histos2d_["MHT_Vs_HT"] = dir.make<TH2D>("MHT_Vs_HT","",100,0.,2000.,100,0.,1000.);
	
    }

    // -------------------- AlphaT --------------------

    { 
      
      TFileDirectory dir = fs->mkdir("AlphaT");

      histos2d_["MHT/HT_Vs_DHT/HT"] = dir.make<TH2D>("MHT/HT_Vs_DHT/HT","",100,0.,1.,100,0.,1.);
      
      histos_["HT"] = dir.make<TH1D>("HT","",100,0.,2000.);
      histos2d_["HT_Vs_NObjects"] = dir.make<TH2D>("HT_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);
      
      histos_["MT"] = dir.make<TH1D>("MT","",100,0.,2000.);
      histos2d_["MT_Vs_NObjects"] = dir.make<TH2D>("MT_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);
      
      histos_["MHT"] = dir.make<TH1D>("MHT","",100,0.,1000.);
      histos2d_["MHT_Vs_NObjects"] = dir.make<TH2D>("MHT_Vs_NObjects","",51,-0.5,50.5,100,0.,1000.);
      
      histos_["MinDHT"] = dir.make<TH1D>("MinDHT","",100,0.,100.);
      histos2d_["MinDHT_Vs_NObjects"] = dir.make<TH2D>("MinDHT_Vs_NObjects","",51,-0.5,50.5,100,0.,100.);
      histos2d_["MinDHT_Vs_MHT"] = dir.make<TH2D>("MinDHT_Vs_MHT","",100,0.,1000.,100,0.,100.);
      
      histos_["HT-MHT"] = dir.make<TH1D>("HT-MHT","",100,0.,2000.);
      histos2d_["HT-MHT_Vs_NObjects"] = dir.make<TH2D>("HT-MHT_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);
      histos2d_["HT-MHT_Vs_MT"] = dir.make<TH2D>("HT-MHT_Vs_MT","",100,0.,2000.,100,0.,2000.);
      
      histos_["AlphaT"] = dir.make<TH1D>("AlphaT","",300,-0.5,2.5);
      histos2d_["AlphaT_Vs_NObjects"] = dir.make<TH2D>("AlphaT_Vs_NObjects","",51,-0.5,50.5,300,-0.5,2.5);
      
      histos_["BetaT"] = dir.make<TH1D>("BetaT","",300,-0.5,2.5);
      histos2d_["BetaT_Vs_NObjects"] = dir.make<TH2D>("BetaT_Vs_NObjects","",51,-0.5,50.5,300,-0.5,2.5);
      
      histos2d_["BetaT_Vs_AlphaT"] = dir.make<TH2D>("BetaT_Vs_AlphaT","",300,-0.5,2.5,300,-0.5,2.5);

    }

    // -------------------- DPHI --------------------
    
    { 
      
      TFileDirectory dir = fs->mkdir("DeltaPhi");
      
      histos_["DPHI"] = dir.make<TH1D>("DPHI","",200,-5.,5.);
      histos2d_["DPHI_Vs_NObjects"] = dir.make<TH2D>("DPHI_Vs_NObjects","",51,-0.5,50.5,200,-5.,5.);

      histos_["BiasedDPHI"] = dir.make<TH1D>("BiasedDPHI","",200,-5.,5.);
      histos2d_["BiasedDPHI_Vs_NObjects"] = dir.make<TH2D>("BiasedDPHI_Vs_NObjects","",51,-0.5,50.5,200,-5.,5.);
      
    }
    
    // -------------------- Biased AlphaT --------------------

    { 
      
      TFileDirectory dir = fs->mkdir("BiasedAlphaT");

      histos_["BiasedHT"] = dir.make<TH1D>("BiasedHT","",100,0.,2000.);
      histos2d_["BiasedHT_Vs_NObjects"] = dir.make<TH2D>("BiasedHT_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);
      
      histos_["BiasedMT"] = dir.make<TH1D>("BiasedMT","",100,0.,2000.);
      histos2d_["BiasedMT_Vs_NObjects"] = dir.make<TH2D>("BiasedMT_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);
      
      histos_["BiasedMHT"] = dir.make<TH1D>("BiasedMHT","",100,0.,1000.);
      histos2d_["BiasedMHT_Vs_NObjects"] = dir.make<TH2D>("BiasedMHT_Vs_NObjects","",51,-0.5,50.5,100,0.,1000.);
      
      histos_["BiasedMinDHT"] = dir.make<TH1D>("BiasedMinDHT","",100,0.,100.);
      histos2d_["BiasedMinDHT_Vs_NObjects"] = dir.make<TH2D>("BiasedMinDHT_Vs_NObjects","",51,-0.5,50.5,100,0.,100.);
      histos2d_["BiasedMinDHT_Vs_BiasedMHT"] = dir.make<TH2D>("BiasedMinDHT_Vs_BiasedMHT","",100,0.,1000.,100,0.,100.);
      
      histos_["BiasedHT-BiasedMHT"] = dir.make<TH1D>("BiasedHT-BiasedMHT","",100,0.,2000.);
      histos2d_["BiasedHT-BiasedMHT_Vs_NObjects"] = dir.make<TH2D>("BiasedHT-BiasedMHT_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);
      histos2d_["BiasedHT-BiasedMHT_Vs_BiasedMT"] = dir.make<TH2D>("BiasedHT-BiasedMHT_Vs_BiasedMT","",100,0.,2000.,100,0.,2000.);

      histos_["BiasedAlphaT"] = dir.make<TH1D>("BiasedAlphaT","",300,-0.5,2.5);
      histos2d_["BiasedAlphaT_Vs_NJets"] = dir.make<TH2D>("BiasedAlphaT_Vs_NJets","",51,-0.5,50.5,300,-0.5,2.5);

      histos_["BiasedBetaT"] = dir.make<TH1D>("BiasedBetaT","",300,-0.5,2.5);
      histos2d_["BiasedBetaT_Vs_NJets"] = dir.make<TH2D>("BiasedBetaT_Vs_NJets","",51,-0.5,50.5,300,-0.5,2.5);

      histos2d_["BiasedAlphaT_Vs_AlphaT"] = dir.make<TH2D>("BiasedAlphaT_Vs_AlphaT","",300,-0.5,2.5,300,-0.5,2.5);
      histos2d_["BiasedBetaT_Vs_BetaT"] = dir.make<TH2D>("BiasedBetaT_Vs_BetaT","",300,-0.5,2.5,300,-0.5,2.5);

    }

    // -------------------- Recoil --------------------
    
    { 
      
      TFileDirectory dir = fs->mkdir("Recoil");
      
      histos_["DPHI_BetweenJetAndPhotonSystems"] = dir.make<TH1D>("DPHI_BetweenJetAndPhotonSystems","",100,-5.,5.); 
      histos_["DPHI_BetweenJetSystemAndRecoilAxis"] = dir.make<TH1D>("DPHI_BetweenJetSystemAndRecoilAxis","",100,-5.,5.); 
      histos_["DPHI_BetweenMhtAndRecoilAxis"] = dir.make<TH1D>("DPHI_BetweenMhtAndRecoilAxis","",100,-5.,5.); 

      histos_["RecoilAlphaT"] = dir.make<TH1D>("RecoilAlphaT","",300,-0.5,2.5);
      histos_["RecoilBetaT"] = dir.make<TH1D>("RecoilBetaT","",300,-0.5,2.5);
      
      histos2d_["HT_Photons_Vs_Jets"] = dir.make<TH2D>("HT_Photons_Vs_Jets","",100,0.,2000.,100,0.,2000.);
      
      histos2d_["DHT_Vs_MHT"] = dir.make<TH2D>("DHT_Vs_MHT","",100,0.,2000.,100,-1000.,1000.);

      histos2d_["DHT_Vs_NJets"] = dir.make<TH2D>("DHT_Vs_NJets","",51,-0.5,50.5,100,-1000.,1000.);
      histos2d_["MHT_Vs_NJets"] = dir.make<TH2D>("MHT_Vs_NJets","",51,-0.5,50.5,100,0.,2000.);

      histos2d_["DHT_Vs_MHT/HT"] = dir.make<TH2D>("DHT_Vs_MHT/HT","",100,0.,1.,100,-1000.,1000.);
      histos2d_["DHT/MHT_Vs_HT"] = dir.make<TH2D>("DHT/MHT_Vs_HT","",100,0.,1000.,200,-20.,20.);
      //histos2d_["DPHI_Vs_MHT/HT"] = dir.make<TH2D>("DPHI_Vs_MHT/HT","",100,0.,1.,100,-5.,5.);

      histos2d_["RecoilAlphaT_Vs_NJets"] = dir.make<TH2D>("RecoilAlphaT_Vs_NJets","",51,-0.5,50.5,300,-0.5,2.5);
      histos2d_["RecoilBetaT_Vs_NJets"] = dir.make<TH2D>("RecoilBetaT_Vs_NJets","",51,-0.5,50.5,300,-0.5,2.5);
      
    }

    // -------------------- Projected --------------------
    
    { 
      
      TFileDirectory dir = fs->mkdir("Projected");

      histos_["ProjectedAlphaT"] = dir.make<TH1D>("ProjectedAlphaT","",300,-0.5,2.5);
      histos_["ProjectedBetaT"] = dir.make<TH1D>("ProjectedBetaT","",300,-0.5,2.5);
      
      histos2d_["ProjectedDHT_Vs_MHT"] = dir.make<TH2D>("ProjectedDHT_Vs_MHT","",100,0.,2000.,100,-1000.,1000.);

      histos2d_["ProjectedDHT_Vs_NJets"] = dir.make<TH2D>("ProjectedDHT_Vs_NJets","",51,-0.5,50.5,100,-1000.,1000.);

      histos2d_["ProjectedDHT_Vs_MHT/HT"] = dir.make<TH2D>("ProjectedDHT_Vs_MHT/HT","",100,0.,1.,100,-1000.,1000.);
      histos2d_["ProjectedDHT/MHT_Vs_HT"] = dir.make<TH2D>("ProjectedDHT/MHT_Vs_HT","",100,0.,1000.,200,-20.,20.);

      histos2d_["ProjectedAlphaT_Vs_NJets"] = dir.make<TH2D>("ProjectedAlphaT_Vs_NJets","",51,-0.5,50.5,300,-0.5,2.5);
      histos2d_["ProjectedBetaT_Vs_NJets"] = dir.make<TH2D>("ProjectedBetaT_Vs_NJets","",51,-0.5,50.5,300,-0.5,2.5);
      
    }
    
    // -------------------- Super --------------------

    { 
      
      TFileDirectory dir = fs->mkdir("Super");
      
      histos_["SuperMinDHT"] = dir.make<TH1D>("SuperMinDHT","",100,0.,100.);
      histos2d_["SuperMinDHT_Vs_NObjects"] = dir.make<TH2D>("SuperMinDHT_Vs_NObjects","",51,-0.5,50.5,100,0.,100.);
      
      histos_["SuperAlphaT"] = dir.make<TH1D>("SuperAlphaT","",300,-0.5,2.5);
      histos2d_["SuperAlphaT_Vs_NObjects"] = dir.make<TH2D>("SuperAlphaT_Vs_NObjects","",51,-0.5,50.5,300,-0.5,2.5);
      
      histos_["SuperAlphaT_NPhotons"] = dir.make<TH1D>("SuperAlphaT_NPhotons","",3,-0.5,2.5);
      histos2d_["SuperAlphaT_NPhotons_Vs_NObjects"] = dir.make<TH2D>("SuperAlphaT_NPhotons_Vs_NObjects","",51,-0.5,50.5,3,-0.5,2.5);
      
      histos_["SuperBetaT"] = dir.make<TH1D>("SuperBetaT","",300,-0.5,2.5);
      histos2d_["SuperBetaT_Vs_NObjects"] = dir.make<TH2D>("SuperBetaT_Vs_NObjects","",51,-0.5,50.5,300,-0.5,2.5);
      
      histos2d_["SuperBetaT_Vs_SuperAlphaT"] = dir.make<TH2D>("SuperBetaT_Vs_SuperAlphaT","",300,-0.5,2.5,300,-0.5,2.5);

    }
    
    // -------------------- Pseudo-dijets --------------------
    
    { 
      
      TFileDirectory dir = fs->mkdir("PseudoDijets");
      
      histos2d_["PseudoDijets_NObjects"] = dir.make<TH2D>("PseudoDijets_NObjects","",51,-0.5,50.5,51,-0.5,50.5); 
      histos2d_["PseudoDijets_NPhotons"] = dir.make<TH2D>("PseudoDijets_NPhotons","",51,-0.5,50.5,51,-0.5,50.5); 
      histos2d_["PseudoDijets_NJets"] = dir.make<TH2D>("PseudoDijets_NJets","",51,-0.5,50.5,51,-0.5,50.5); 
      
      histos2d_["PseudoDijets_PhotonEt"] = dir.make<TH2D>("PseudoDijets_PhotonEt","",100,0.,1000.,100,0.,1000.);
      histos2d_["PseudoDijets_DeltaPhi_Vs_PhotonEt"] = dir.make<TH2D>("PseudoDijets_DeltaPhi_Vs_PhotonEt","",100,0.,1000.,200,-5.,5.);
      histos2d_["PseudoDijets_DeltaPhi_Vs_NObjects"] = dir.make<TH2D>("PseudoDijets_DeltaPhi_Vs_NObjects","",51,-0.5,50.5,200,-5.,5.);
      
    }

  } else { 
    edm::LogWarning("TEST") << "TEST! Number of objects as input: " << test_; 
  }
  
}

// -----------------------------------------------------------------------------
//
bool TestCombination::getPhotons( const edm::Event& iEvent,
				  std::vector<Candidate>& photons,
				  double threshold ) {
  
  photons.clear();
  
  if ( !photons_.label().empty() ) {
    
    edm::Handle< std::vector<pat::Photon> > handle;
    iEvent.getByLabel(photons_,handle);
    
    if ( !handle.isValid() ) { 
      edm::LogWarning("TEST") << "No Photons for " << photons_; 
      return true;
    }

    double photon_et = threshold < 0. ? photonEt_ : threshold;
    
    uint16_t et25=0, et30=0, et35=0, et40=0, et45=0, et50=0;
    std::vector<pat::Photon>::const_iterator iphoton = handle->begin();
    std::vector<pat::Photon>::const_iterator jphoton = handle->end();
    for ( ; iphoton != jphoton; ++iphoton  ) { 
      if ( fabs( iphoton->eta() ) < photonEta_ ) { 
	if ( iphoton->et() > photon_et ) { photons.push_back( *iphoton ); }
	if ( iphoton->et() > 25. ) { et25++; }
	if ( iphoton->et() > 30. ) { et30++; }
	if ( iphoton->et() > 35. ) { et35++; }
	if ( iphoton->et() > 40. ) { et40++; }
	if ( iphoton->et() > 45. ) { et45++; }
	if ( iphoton->et() > 50. ) { et50++; }
      }
    }
    
    //     histo("PreCuts_NPhotons25")->Fill( et25 ); 
    //     histo("PreCuts_NPhotons30")->Fill( et30 ); 
    //     histo("PreCuts_NPhotons35")->Fill( et35 ); 
    //     histo("PreCuts_NPhotons40")->Fill( et40 ); 
    //     histo("PreCuts_NPhotons45")->Fill( et45 ); 
    //     histo("PreCuts_NPhotons50")->Fill( et50 ); 
    
    if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Photons: " << photons.size(); }

  }
  
  return false; 

}

// -----------------------------------------------------------------------------
//
bool TestCombination::getJets( const edm::Event& iEvent,
			       std::vector<Candidate>& jets,
			       double threshold ) {
  
  jets.clear();
  
  if ( !jets_.label().empty() ) {

    edm::Handle< std::vector<pat::Jet> > handle;
    iEvent.getByLabel(jets_,handle);
    
    if ( !handle.isValid() ) { 
      edm::LogWarning("TEST") << "No Jets for " << jets_; 
      return true;
    }

    double jet_et = threshold < 0. ? jetEt_ : threshold;
    
    uint16_t et25=0, et30=0, et35=0, et40=0, et45=0, et50=0;
    std::vector<pat::Jet>::const_iterator ijet = handle->begin();
    std::vector<pat::Jet>::const_iterator jjet = handle->end();
    for ( ; ijet != jjet; ++ijet  ) { 
      if ( fabs( ijet->eta() ) < jetEta_ && 
	   ijet->emEnergyFraction() < jetEMfrac_ ) { 
	if ( ijet->et() > jet_et ) { jets.push_back( *ijet ); }
	if ( ijet->et() > 25. ) { et25++; }
	if ( ijet->et() > 30. ) { et30++; }
	if ( ijet->et() > 35. ) { et35++; }
	if ( ijet->et() > 40. ) { et40++; } 
	if ( ijet->et() > 45. ) { et45++; }
	if ( ijet->et() > 50. ) { et50++; }
      }
    }
    
    //     histo("PreCuts_NJets25")->Fill( et25 ); 
    //     histo("PreCuts_NJets30")->Fill( et30 ); 
    //     histo("PreCuts_NJets35")->Fill( et35 ); 
    //     histo("PreCuts_NJets40")->Fill( et40 ); 
    //     histo("PreCuts_NJets45")->Fill( et45 ); 
    //     histo("PreCuts_NJets50")->Fill( et50 ); 

    if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Jets: " << jets.size(); }

  }
  
  return false;
  
}

// -----------------------------------------------------------------------------
//
bool TestCombination::getMuons( const edm::Event& iEvent,
				std::vector<Candidate>& muons ) {
  
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

  }

  return false;
  
}

// -----------------------------------------------------------------------------
//
bool TestCombination::getElectrons( const edm::Event& iEvent,
				    std::vector<Candidate>& electrons ) {
  
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

  }
  
  return false;
  
}

// -----------------------------------------------------------------------------
//
double TestCombination::HT( const std::vector<Candidate>& input ) {
  LorentzV lv;
  double ht = 0.;
  double et = 0.;
  std::vector<Candidate>::const_iterator ii = input.begin();
  std::vector<Candidate>::const_iterator jj = input.end();
  for ( ; ii != jj; ++ii ) { 
    lv += ii->p4();
    et += ii->p4().Et(); 
    double tmp = ii->p4().E()*ii->p4().E() - ii->p4().Pz()*ii->p4().Pz();
    ht += tmp < 0. ? -1.*sqrt(-1.*tmp) : sqrt(tmp);
  }
  double ht_ = lv.E()*lv.E() - lv.Pz()*lv.Pz();
  ht_ = ht_ < 0. ? -1.*sqrt(-1.*ht_) : sqrt(ht_);
  cout << " TEMP ht= " << ht << " et= " << et << " ht_= " << ht_ << endl; 
  return ht;
}

// -----------------------------------------------------------------------------
//
double TestCombination::MHT( const std::vector<Candidate>& input ) {
  LorentzV lv_obj;
  std::vector<Candidate>::const_iterator iobj = input.begin();
  std::vector<Candidate>::const_iterator jobj = input.end();
  for ( ; iobj != jobj; ++iobj ) { lv_obj += iobj->p4(); }
  LorentzV lv_mht( lv_obj );
  lv_mht.SetPx( -1.*lv_obj.Px() );
  lv_mht.SetPy( -1.*lv_obj.Py() );
  lv_mht.SetPz( -1.*lv_obj.Pz() );
  return sqrt( lv_mht.Perp2() );
}

// -----------------------------------------------------------------------------
//
double TestCombination::MT( const std::vector<Candidate>& input ) {
  double ht  = HT( input );
  double mht = MHT( input );
  double mt  = ht*ht - mht*mht;
  return mt < 0. ? -1.*sqrt(-1.*mt) : sqrt(mt);
}

// -----------------------------------------------------------------------------
//
double TestCombination::HT_MHT( const std::vector<Candidate>& input ) {
  double ht  = HT( input );
  double mht = MHT( input );
  double ht_mht = ht - mht;
  return ht_mht;
}

// -----------------------------------------------------------------------------
//
double TestCombination::minDHT( const std::vector<Candidate>& objects, 
				std::vector<Candidate>& jet1, 
				std::vector<Candidate>& jet2  ) {
  
  // Maximum number of objects handled
  static const uint8_t max_size = 50;

  // Init
  jet1.clear();
  jet2.clear();
  double min_delta_ht = -1.;
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
    min_delta_ht = 0.;
    combinations = 1;
    
  } else if ( size > 0 ) {

    // Cached data
    double min_diff_ht = -1.;
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
	double et1 = 0.;
	std::vector<uint8_t> tmp1;
	tmp1.reserve(size2+1);
	for ( uint8_t ii = 0; ii < size2; ++ii ) { 
	  uint8_t index = static_cast<uint8_t>( c2[ii] );
	  tmp1.push_back( index ); 
	  et1 += Et( objects[index].p4() );
	}

	// Indices of jets assigned to second pseudo-jet and Et sum
	double et2 = 0.;
	std::vector<uint8_t> tmp2;
	tmp2.reserve(size1+1);
	for ( uint8_t ii = 0; ii < size1; ++ii ) { 
	  if ( std::find( tmp1.begin(), tmp1.end(), ii ) == tmp1.end() ) { 
	    tmp2.push_back(ii); 
	    et2 += Et( objects[ii].p4() );
	  }
	}
	
	// Calculate difference in Et between two pseudo-jets
	if ( min_diff_ht < 0. || fabs( et1 - et2 ) < min_diff_ht ) { 
	  min_diff_ht = fabs( et1 - et2 ); 
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
    if ( min_diff_ht < 0. ) {
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

    min_delta_ht = min_diff_ht;

  } else {
    
    jet1.clear();
    jet2.clear();
    min_delta_ht = -1.;
    
  }

  return min_delta_ht;
  
}

// -----------------------------------------------------------------------------
//
double TestCombination::recoilDHT( const std::vector<Candidate>& photons, 
				   const std::vector<Candidate>& jets, 
				   std::vector<Candidate>& jet1, 
				   std::vector<Candidate>& jet2 ) {

  std::stringstream ss;
  ss << " TEST " << endl;
  
  // Init
  jet1.clear();
  jet2.clear();
  double delta_ht = -1.;

  // Check objects used  
  if ( photons.size() < 2 || jets.empty() ) { 
    delta_ht = -1.;
  } else {
    
    // Photon system
    double et_pho = 0.;
    LorentzV lv_pho;
    std::vector<Candidate>::const_iterator ipho = photons.begin();
    std::vector<Candidate>::const_iterator jpho = photons.begin() + 2; 
    ss << " photons: " << endl;
    for ( ; ipho != jpho; ++ipho ) { 
      std::stringstream sss; sss << ipho->p4(); 
      ss << " #" << uint32_t(ipho-photons.begin()) << " " << sss.str().substr(sss.str().find_first_of("(")+1,sss.str().find_first_of(")")-1) << endl;
      lv_pho += ipho->p4(); 
      et_pho += Et( ipho->p4() );
    }
    ss << endl;
    
    // "Recoil" axis
    LorentzV lv_recoil( lv_pho );
    lv_recoil.SetPx( -1.*lv_recoil.Px() );
    lv_recoil.SetPy( -1.*lv_recoil.Py() );
    lv_recoil.SetPz( -1.*lv_recoil.Pz() );
    
    // Self dot product of recoil 4-momentum
    LorentzV::Scalar self_product = lv_recoil.Dot( lv_recoil );
    
    // Jet system
    double et_jet = 0.;
    LorentzV lv_jet;
    std::vector<Candidate>::const_iterator ijet = jets.begin();
    std::vector<Candidate>::const_iterator jjet = jets.end();
    ss << " jets: " << endl;
    for ( ; ijet != jjet; ++ijet ) { 
      std::stringstream sss; sss << ijet->p4(); 
      ss << " #" << uint32_t(ijet-jets.begin()) << " " << sss.str().substr(sss.str().find_first_of("(")+1,sss.str().find_first_of(")")-1) << endl;
      lv_jet += ijet->p4(); 
      LorentzV::Scalar dot_product = ijet->p4().Dot( lv_recoil );
      LorentzV lv_proj = lv_recoil * ( dot_product / self_product );
      et_jet += Et( lv_proj );
    }
    ss << endl;
    
    // Add "other" lower-Et photons to jet system
    if ( photons.size() > 2 ) {
      std::vector<Candidate>::const_iterator iother = photons.begin() + 2;
      std::vector<Candidate>::const_iterator jother = photons.end();
      ss << " other photons: " << endl;
      for ( ; iother < jother; ++iother ) { 
	std::stringstream sss; sss << iother->p4(); 
	ss << " #" << uint32_t(iother-(photons.begin()+2)) << " " << sss.str().substr(sss.str().find_first_of("(")+1,sss.str().find_first_of(")")-1) << endl;
	lv_jet += iother->p4(); 
	LorentzV::Scalar dot_product = iother->p4().Dot( lv_recoil );
	LorentzV lv_proj = lv_recoil * ( dot_product / self_product );
	et_jet += Et( lv_proj );
// 	ss << " projection: " << endl;
// 	ss << lv_proj.Px() << " " 
// 	   << lv_proj.Py() << " " 
// 	   << lv_proj.Pz() << " " 
// 	   << lv_proj.E() << " " 
// 	   << lv_proj.Perp2();
// 	ss << endl;
//     ss << " recoil.recoil = " << self_product
//        << " recoil.jet = " << dot_product
//        << endl;
      }
    }
    ss << endl;
  
    histo2d("HT_Photons_Vs_Jets")->Fill( HT( jets ), HT( photons ) ); 
    
    // MHT from complete system
    LorentzV lv_mht;
    lv_mht = lv_pho + lv_jet;
    lv_mht.SetPx( -1.*lv_mht.Px() );
    lv_mht.SetPy( -1.*lv_mht.Py() );
    lv_mht.SetPz( -1.*lv_mht.Pz() );
    
    ss << " recoil: " << endl;
    ss << lv_recoil.Px() << " " 
       << lv_recoil.Py() << " " 
       << lv_recoil.Pz() << " " 
       << lv_recoil.E() << " " 
       << lv_recoil.Perp2() << " "
       << lv_recoil.Et() << " "
       << sqrt( lv_recoil.E()*lv_recoil.E() - lv_recoil.Pz()*lv_recoil.Pz() );
    ss << endl;
      
    ss << " mht: " << endl;
    ss << lv_mht.Px() << " " 
       << lv_mht.Py() << " " 
       << lv_mht.Pz() << " " 
       << lv_mht.E() << " " 
       << lv_mht.Perp2();
    ss << endl;

    cout << ss.str() << endl;

    // DPHI b/w two systems
    {
      //double pi = 4. * atan(1.);
      double delta_phi = reco::deltaPhi( lv_pho.phi(), lv_jet.phi() );
      histo("DPHI_BetweenJetAndPhotonSystems")->Fill( fabs(delta_phi) ); 
    }
      
    // DPHI b/w jet system and recoil axis
    {
      double delta_phi = reco::deltaPhi( lv_jet.phi(), lv_recoil.phi() );
      histo("DPHI_BetweenJetSystemAndRecoilAxis")->Fill( fabs(delta_phi) ); 
    }

    // DPHI b/w MHT and recoil axis
    {
      double delta_phi = reco::deltaPhi( lv_recoil.phi(), lv_mht.phi() );
      histo("DPHI_BetweenMhtAndRecoilAxis")->Fill( fabs(delta_phi) ); 
    }

    delta_ht = fabs( et_pho - et_jet );
    
    //     // Store pseudo-jets
    //     if ( !photons.empty() )   { jet1.insert( jet1.end(), photons.begin(), photons.begin()+2 ); }
    //     if ( !jets.empty() )      { jet2.insert( jet2.end(), jets.begin(), jets.end() ); }
    //     if ( photons.size() > 2 ) { jet2.insert( jet2.end(), photons.begin()+2, photons.end() ); }
    
  }
  
  return delta_ht;
  
}

// -----------------------------------------------------------------------------
//
double TestCombination::minSuperDHT( const std::vector<Candidate>& objects, 
				     std::vector<Candidate>& jet1, 
				     std::vector<Candidate>& jet2 ) {
  
  // Maximum number of objects handled
  static const uint8_t max_size = 50;

  // Init
  jet1.clear();
  jet2.clear();
  double min_delta_ht = -1.;
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
  
  if ( size == 0 ) { 
    
    jet1.clear();
    jet2.clear();
    combinations = 0;
    min_delta_ht = -1.;
    
  } else if ( size == 1 ) {
    
    jet1.push_back( objects[0] );
    jet2.clear();
    min_delta_ht = -1.;
    combinations = 1;
    
  } else if ( size == 2 ) {
    
    jet1.push_back( objects[0] );
    jet1.push_back( objects[1] );
    jet2.clear();
    min_delta_ht = -1.;
    combinations = 1;
    
  } else {
    
    edm::LogVerbatim("TEST") << "GET HERE! size = " << (uint16_t)size;
    
    // Cached data
    double min_diff_ht = -1.;
    std::vector<uint8_t> indices1;
    std::vector<uint8_t> indices2;
    
    // Build char array encoding indices of all objects
    std::vector<char> vc1;
    vc1.reserve(size+1);
    for ( uint8_t ii = 0; ii < size; ++ii ) { vc1.push_back(ii); }
    vc1[vc1.size()] = '\0';
    char c1[256];
    memcpy( c1, &vc1.front(), vc1.size() );

    // Build char array encoding indices of subset of objects
    int number = 2;
    std::vector<char> vc2;
    vc2.reserve(number+1);
    for ( uint8_t jjj = 0; jjj < number; ++jjj ) { vc2.push_back(jjj); }
    vc2[vc2.size()] = '\0';
    char c2[256];
    memcpy( c2, &vc2.front(), vc2.size() );
      
    // Size of both char arrays
    uint8_t size1 = vc1.size();
    uint8_t size2 = vc2.size();

    // Iterate through combinations
    do { 

      // Indices of jets assigned to first pseudo-jet and LorentzV sum
      double et1 = 0.;
      LorentzV lv1;
      std::vector<uint8_t> tmp1;
      tmp1.reserve(size2+1);
      for ( uint8_t ii = 0; ii < size2; ++ii ) { 
	uint8_t index = static_cast<uint8_t>( c2[ii] );
	tmp1.push_back( index ); 
	lv1 += objects[index].p4();
	et1 += Et( objects[index].p4() );
      }

      // "Recoil" axis
      LorentzV lv_recoil( lv1 );
      lv_recoil.SetPx( -1.*lv_recoil.Px() );
      lv_recoil.SetPy( -1.*lv_recoil.Py() );
      lv_recoil.SetPz( -1.*lv_recoil.Pz() );

      // Self dot product of recoil 4-momentum
      LorentzV::Scalar self_product = lv_recoil.Dot( lv_recoil );
      
      // Indices of jets assigned to second pseudo-jet and LorentzV sum
      double et2 = 0.;
      LorentzV lv2;
      std::vector<uint8_t> tmp2;
      tmp2.reserve(size1+1);
      for ( uint8_t ii = 0; ii < size1; ++ii ) { 
	if ( std::find( tmp1.begin(), tmp1.end(), ii ) == tmp1.end() ) { 
	  tmp2.push_back(ii); 
	  lv2 += objects[ii].p4();
	  LorentzV::Scalar dot_product = objects[ii].p4().Dot( lv_recoil );
	  LorentzV lv_proj = lv_recoil * ( dot_product / self_product );
	  et1 += Et( lv_proj );
	}
      }

      { // Debug
	std::stringstream ss;
	ss << "Combination #" << combinations << " = {";
	std::vector<uint8_t>::const_iterator ii1 = tmp1.begin();
	std::vector<uint8_t>::const_iterator jj1 = tmp1.end();
	for ( ; ii1 != jj1; ++ii1 ) { ss << (uint16_t)*ii1; }
	ss << ",";
	std::vector<uint8_t>::const_iterator ii2 = tmp2.begin();
	std::vector<uint8_t>::const_iterator jj2 = tmp2.end();
	for ( ; ii2 != jj2; ++ii2 ) { ss << (uint16_t)*ii2; }
	ss << "}";
	edm::LogVerbatim("TEST") << ss.str();
      }
      
      // Calculate difference in Et between two pseudo-jets
      if ( min_diff_ht < 0. || fabs( et1 - et2 ) < min_diff_ht ) { 
	min_diff_ht = fabs( et1 - et2 ); 
	indices1.resize( tmp1.size() );
	indices2.resize( tmp2.size() );
	std::copy( tmp1.begin(), tmp1.end(), indices1.begin() );
	std::copy( tmp2.begin(), tmp2.end(), indices2.begin() );
      }
	
      combinations++;

    } while ( stdcomb::next_combination( c1, c1 + size1, 
					 c2, c2 + size2 ) );
    
    // Build pseudo-jets
    if ( min_diff_ht < 0. ) {
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

    min_delta_ht = min_diff_ht;    
    
  }

  return min_delta_ht;
  
}

// -----------------------------------------------------------------------------
//
double TestCombination::superDPHI( const std::vector<Candidate>& objects, 
				   std::vector<Candidate>& jet1, 
				   std::vector<Candidate>& jet2 ) {
  
  // Maximum number of objects handled
  static const uint8_t max_size = 50;

  // Init
  jet1.clear();
  jet2.clear();
  double min_delta_phi = -1.;
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
  
  if ( size == 0 ) { 
    
    jet1.clear();
    jet2.clear();
    combinations = 0;
    min_delta_phi = -1.;
    
  } else if ( size == 1 ) {
    
    jet1.push_back( objects[0] );
    jet2.clear();
    min_delta_phi = -1.;
    combinations = 1;
    
  } else if ( size == 2 ) {
    
    jet1.push_back( objects[0] );
    jet1.push_back( objects[1] );
    jet2.clear();
    min_delta_phi = -1.;
    combinations = 1;
    
  } else {
    
    edm::LogVerbatim("TEST") << "GET HERE! size = " << (uint16_t)size;
    
    // Cached data
    std::vector<uint8_t> indices1;
    std::vector<uint8_t> indices2;
    
    // Build char array encoding indices of all objects
    std::vector<char> vc1;
    vc1.reserve(size+1);
    for ( uint8_t ii = 0; ii < size; ++ii ) { vc1.push_back(ii); }
    vc1[vc1.size()] = '\0';
    char c1[256];
    memcpy( c1, &vc1.front(), vc1.size() );

    // Build char array encoding indices of subset of objects
    int number = 2;
    std::vector<char> vc2;
    vc2.reserve(number+1);
    for ( uint8_t jjj = 0; jjj < number; ++jjj ) { vc2.push_back(jjj); }
    vc2[vc2.size()] = '\0';
    char c2[256];
    memcpy( c2, &vc2.front(), vc2.size() );
      
    // Size of both char arrays
    uint8_t size1 = vc1.size();
    uint8_t size2 = vc2.size();

    // Iterate through combinations
    do { 

      // Indices of jets assigned to first pseudo-jet and LorentzV sum
      LorentzV lv1;
      std::vector<uint8_t> tmp1;
      tmp1.reserve(size2+1);
      for ( uint8_t ii = 0; ii < size2; ++ii ) { 
	uint8_t index = static_cast<uint8_t>( c2[ii] );
	tmp1.push_back( index ); 
	lv1 += objects[index].p4();
      }
      
      // Indices of jets assigned to second pseudo-jet and LorentzV sum
      LorentzV lv2;
      std::vector<uint8_t> tmp2;
      tmp2.reserve(size1+1);
      for ( uint8_t ii = 0; ii < size1; ++ii ) { 
	if ( std::find( tmp1.begin(), tmp1.end(), ii ) == tmp1.end() ) { 
	  tmp2.push_back(ii); 
	  lv2 += objects[ii].p4();
	}
      }

      { // Debug
	std::stringstream ss;
	ss << "Combination #" << combinations << " = {";
	std::vector<uint8_t>::const_iterator ii1 = tmp1.begin();
	std::vector<uint8_t>::const_iterator jj1 = tmp1.end();
	for ( ; ii1 != jj1; ++ii1 ) { ss << (uint16_t)*ii1; }
	ss << ",";
	std::vector<uint8_t>::const_iterator ii2 = tmp2.begin();
	std::vector<uint8_t>::const_iterator jj2 = tmp2.end();
	for ( ; ii2 != jj2; ++ii2 ) { ss << (uint16_t)*ii2; }
	ss << "}";
	edm::LogVerbatim("TEST") << ss.str();
      }

      // Calculate missing HT
      LorentzV lv_mht;
      lv_mht = lv1 + lv2;
      lv_mht.SetPx( -1.*lv_mht.Px() );
      lv_mht.SetPy( -1.*lv_mht.Py() );
      lv_mht.SetPz( -1.*lv_mht.Pz() );

      // Calculate delta phi b/w MHT and pseudo-jet #1
      { 
	double delta_phi = reco::deltaPhi<LorentzV,LorentzV>( lv_mht, lv1 );
	if ( fabs(delta_phi) < fabs(min_delta_phi) ) { min_delta_phi = delta_phi; }
      }

      // Calculate delta phi b/w MHT and pseudo-jet #2
      {
	double delta_phi = reco::deltaPhi<LorentzV,LorentzV>( lv_mht, lv2 );
	if ( fabs(delta_phi) < fabs(min_delta_phi) ) { min_delta_phi = delta_phi; }
      }

      combinations++;

    } while ( stdcomb::next_combination( c1, c1 + size1, 
					 c2, c2 + size2 ) );
    
    // Build pseudo-jets
    if ( min_delta_phi < 0. ) {
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
    
  }

  return min_delta_phi;
  
}

// -----------------------------------------------------------------------------
//
double TestCombination::alphaT( const std::vector<Candidate>& objects,
				std::vector<Candidate>& jet1, 
				std::vector<Candidate>& jet2 ) {
  jet1.clear();
  jet2.clear();
  double min_dht = minDHT( objects, jet1, jet2 );
  double ht = HT( objects );
  double mt = MT( objects );
  return alphaT( min_dht, ht, mt );
}

// -----------------------------------------------------------------------------
//
double TestCombination::alphaT( double min_dht, double ht, double mt ) {
  return 0.5 * ( ( ht - min_dht ) / mt );
}
  
// -----------------------------------------------------------------------------
//
TH1D* TestCombination::histo( const std::string& histogram_name ) {
  std::map<std::string,TH1D*>::const_iterator ii = histos_.find(histogram_name);
  if ( ii != histos_.end() ) { return ii->second; }
  edm::LogWarning("TEST") << "Cannot find string: " << histogram_name;
  return 0;
}
  
// -----------------------------------------------------------------------------
//
TH2D* TestCombination::histo2d( const std::string& histogram_name ) {
  std::map<std::string,TH2D*>::const_iterator jj = histos2d_.find(histogram_name);
  if ( jj != histos2d_.end() ) { return jj->second; }
  edm::LogWarning("TEST") << "Cannot find string: " << histogram_name;
  return 0;
}

// -----------------------------------------------------------------------------
// 
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TestCombination);













//       // DPHI b/w two systems
//       {
// 	double delta_phi = reco::deltaPhi<LorentzV,LorentzV>( lv_pho, lv_jet );
// 	histo("DPHI_BetweenJetAndPhotonSystems")->Fill( fabs(delta_phi) ); 
//       }
      
//       // DPHI b/w jet system and recoil axis
//       {
// 	double delta_phi = reco::deltaPhi<LorentzV,LorentzV>(lv_jet,lv_recoil);
// 	histo("DPHI_BetweenJetSystemAndRecoilAxis")->Fill( fabs(delta_phi) ); 
//       }

//       // DPHI b/w MHT and recoil axis
//       {
// 	double delta_phi = reco::deltaPhi<LorentzV,LorentzV>(lv_mht,lv_recoil);
// 	histo("DPHI_BetweenMhtAndRecoilAxis")->Fill( fabs(delta_phi) ); 
//       }

//       // "Biased" alphaT
//       {

// 	double ht  = HT(objects);
// 	double mht = sqrt( lv_mht.Perp2() );
// 	double mt  = ht*ht - mht*mht;
// 	mt = mt < 0. ? -1.*sqrt(-1.*mt) : sqrt(mt);
// 	double ht_mht = ht - mht;
// 	double dht = 









//     {
//       std::vector<int>::const_iterator iii1 = phot1.begin();
//       std::vector<int>::const_iterator jjj1 = phot1.end();
//       for ( ; iii1 != jjj1; ++iii1 ) {
// 	Candidate* candidate = 0;
// 	int32_t invalid_id = static_cast<int32_t>(1e8+0.5);
// 	int32_t daughter_id = invalid_id;
// 	int32_t mother_id   = invalid_id;
// 	std::vector<reco::GenParticleRef> gen_particles = iphoton->genParticleRefs();
// 	std::vector<reco::GenParticleRef>::const_iterator igen = gen_particles.begin();
// 	std::vector<reco::GenParticleRef>::const_iterator jgen = gen_particles.end();
// 	for ( ; igen != jgen; ++igen ) {
// 	  if ( igen->isNonnull() ) { 
// 	    reco::GenParticle* part = const_cast<reco::GenParticle*>( &(**igen) );
// 	    Candidate* cand = dynamic_cast<Candidate*>( part );
// 	    if ( cand ) {
// 	      candidate = cand;
// 	      daughter_id = cand->pdgId();
// 	      const Candidate* mother = cand->mother();
// 	      if ( mother && cand->pdgId() == cand->mother()->pdgId() ) { mother = mother->mother(); }
// 	      if ( mother ) { mother_id = mother->pdgId(); }
// 	      break;
// 	    }
// 	  }
// 	}
// 	histo2d("PseudoDijets_ParentPdgId")->Fill( static_cast<uint32_t>( iii1 - phot1.begin() ), mother_id ); 
//       }
//     }

//     {
//       std::vector<int>::const_iterator iii2 = phot2.begin();
//       std::vector<int>::const_iterator jjj2 = phot2.end();
//       for ( ; iii2 != jjj2; ++iii2 ) {
// 	Candidate* candidate = 0;
// 	int32_t invalid_id = static_cast<int32_t>(1e8+0.5);
// 	int32_t daughter_id = invalid_id;
// 	int32_t mother_id   = invalid_id;
// 	std::vector<reco::GenParticleRef> gen_particles = iphoton->genParticleRefs();
// 	std::vector<reco::GenParticleRef>::const_iterator igen = gen_particles.begin();
// 	std::vector<reco::GenParticleRef>::const_iterator jgen = gen_particles.end();
// 	for ( ; igen != jgen; ++igen ) {
// 	  if ( igen->isNonnull() ) { 
// 	    reco::GenParticle* part = const_cast<reco::GenParticle*>( &(**igen) );
// 	    Candidate* cand = dynamic_cast<Candidate*>( part );
// 	    if ( cand ) {
// 	      candidate = cand;
// 	      daughter_id = cand->pdgId();
// 	      const Candidate* mother = cand->mother();
// 	      if ( mother && cand->pdgId() == cand->mother()->pdgId() ) { mother = mother->mother(); }
// 	      if ( mother ) { mother_id = mother->pdgId(); }
// 	      break;
// 	    }
// 	  }
// 	}
// 	histo2d("PseudoDijets_ParentPdgId")->Fill( static_cast<uint32_t>( iii2 - phot2.begin() ), mother_id ); 
//       }
//     }
    


// // -----------------------------------------------------------------------------
// //
// double TestCombination::constructPseudoJets( const std::vector<Candidate>& objects, 
// 					    std::vector<Candidate>& pseudo_jet1, 
// 					    std::vector<Candidate>& pseudo_jet2  ) {
  
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
// double TestCombination::minDeltaet( const std::vector<Candidate>& objects, 
// 				   const std::vector< std::vector<uint8_t> >& combo1, 
// 				   const std::vector< std::vector<uint8_t> >& combo2,
// 				   std::vector<Candidate>& pseudo1, 
// 				   std::vector<Candidate>& pseudo2 ) {
  
//   pseudo1.clear();
//   pseudo2.clear();
  
//   uint32_t size = combo1.size() < combo2.size() ? combo1.size() : combo2.size();
  
//   if ( size == 0 ) { return -1.; }
  
//   double min_delta_pt = -1.;
//   uint32_t index = 0;
//   for ( uint32_t ii = 0; ii < size; ++ii ) {
//     double pt1 = 0.;
//     for ( uint32_t jj = 0; jj < combo1[ii].size(); ++jj ) { pt1 += objects[ combo1[ii][jj] ].pt(); }
//     double pt2 = 0.;
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









//       typedef Candidate::LorentzV LorentzV;
      
//       LorentzV lv_pho;
//       std::vector<Candidate>::const_iterator ipho = photons.begin();
//       std::vector<Candidate>::const_iterator jpho = photons.end();
//       for ( ; ipho != jpho; ++ipho ) { lv_pho += ipho->p4(); }
//       Candidate pho( 22, lv_pho );
      
//       LorentzV lv_jet;
//       std::vector<Candidate>::const_iterator ijet = jets.begin();
//       std::vector<Candidate>::const_iterator jjet = jets.end();
//       for ( ; ijet != jjet; ++ijet ) { lv_jet += ijet->p4(); }
//       Candidate jet( 0, lv_jet );

//       LorentzV lv_all;
//       std::vector<Candidate>::const_iterator iall = objects.begin();
//       std::vector<Candidate>::const_iterator jall = objects.end();
//       for ( ; iall != jall; ++iall ) { lv_all += iall->p4(); }
//       Candidate all( 0, lv_all );

//       {
// 	histo2d("PhiPhotons_Vs_PhiJets")->Fill( lv_jet.phi(), lv_pho.phi() ); 
// 	double delta_phi = reco::deltaPhi<Candidate,Candidate>( pho, jet );
// 	histo("DPHI_PhotonsJets")->Fill( delta_phi ); 

// 	histo2d("EtPhotons_Vs_EtJets")->Fill( lv_jet.Et(), lv_pho.Et() ); 
// 	histo("DeltaEt_PhotonsJets")->Fill( lv_jet.Et() - lv_pho.Et() ); 

// 	histo2d("PtPhotons_Vs_PtJets")->Fill( lv_jet.Pt(), lv_pho.Pt() ); 
// 	histo("DeltaPt_PhotonsJets")->Fill( lv_jet.Pt() - lv_pho.Pt() ); 
//       }

//       {
// 	histo2d("PhiPhotons_Vs_PhiObjects")->Fill( lv_all.phi(), lv_pho.phi() ); 
// 	double delta_phi = reco::deltaPhi<Candidate,Candidate>( pho, all );
// 	histo("DPHI_PhotonsObjects")->Fill( delta_phi ); 

// 	histo2d("EtPhotons_Vs_EtObjects")->Fill( lv_all.Et(), lv_pho.Et() ); 
// 	histo("DeltaEt_PhotonsObjects")->Fill( lv_all.Et() - lv_pho.Et() ); 

// 	histo2d("PtPhotons_Vs_PtObjects")->Fill( lv_all.Pt(), lv_pho.Pt() ); 
// 	histo("DeltaPt_PhotonsObjects")->Fill( lv_all.Pt() - lv_pho.Pt() ); 
//       }

//       {

// 	//LorentzV recoil( -1.*lv_pho.Px(), -1.*lv_pho.Py(), -1.*lv_pho.Pz(), lv_pho.E() );
// 	LorentzV recoil( lv_pho );
// 	recoil.SetPx( -1.*lv_pho.Px() );
// 	recoil.SetPy( -1.*lv_pho.Py() );
// 	recoil.SetPz( -1.*lv_pho.Pz() );
// 	LorentzV::Scalar mag2 = lv_pho.Dot(lv_pho);
	
// 	LorentzV::Scalar dot_product_all = recoil.Dot(lv_all);
// 	LorentzV::Scalar dot_product_jet = recoil.Dot(lv_jet);
	
// 	histo("DotProduct_All")->Fill( dot_product_all ); 
// 	histo("DotProduct_Jets")->Fill( dot_product_jet ); 
// 	histo("DotProduct_Magnitude2")->Fill( mag2 ); 

// 	LorentzV proj_all = recoil * ( dot_product_all / mag2 );
// 	LorentzV proj_jet = recoil * ( dot_product_jet / mag2 );
	
// 	double delta_phi_all = reco::deltaPhi<LorentzV,LorentzV>(lv_all,recoil);
// 	double delta_phi_jet = reco::deltaPhi<LorentzV,LorentzV>(lv_jet,recoil);
	
// 	histo("DeltaEt_All")->Fill( proj_all.Et() - recoil.Et() ); 
// 	histo("DeltaEt_Jets")->Fill( proj_jet.Et() - recoil.Et() ); 
// 	histo2d("Et_Photons_Vs_All")->Fill( proj_all.Et(), recoil.Et() ); 
// 	histo2d("Et_Photons_Vs_Jets")->Fill( proj_jet.Et(), recoil.Et() ); 

// 	histo2d("Et_Proj_Vs_All")->Fill( lv_all.Et(), proj_all.Et() ); 
// 	histo2d("Et_Proj_Vs_Jets")->Fill( lv_jet.Et(), proj_jet.Et() ); 

// 	histo("DeltaPt_All")->Fill( proj_all.Pt() - recoil.Pt() ); 
// 	histo("DeltaPt_Jets")->Fill( proj_jet.Pt() - recoil.Pt() ); 
// 	histo2d("Pt_Photons_Vs_All")->Fill( proj_all.Pt(), recoil.Pt() ); 
// 	histo2d("Pt_Photons_Vs_Jets")->Fill( proj_jet.Pt(), recoil.Pt() ); 

// 	histo2d("Pt_Proj_Vs_All")->Fill( lv_all.Pt(), proj_all.Pt() ); 
// 	histo2d("Pt_Proj_Vs_Jets")->Fill( lv_jet.Pt(), proj_jet.Pt() ); 
	
// 	histo("DPHI_All")->Fill( delta_phi_all ); 
// 	histo("DPHI_Jets")->Fill( delta_phi_jet ); 
// 	histo2d("Phi_Photons_Vs_All")->Fill( lv_all.phi(), recoil.phi() ); 
// 	histo2d("Phi_Photons_Vs_Jets")->Fill( lv_jet.phi(), recoil.phi() ); 
	
//       }




//     eventWeight_( pset.getUntrackedParameter<bool>("ApplyEventWeight",false) ),
//     normLumi_( pset.getParameter<double>("NormalisedLumi") ),
//     xSec_( pset.getParameter<double>("CrossSection") ),
//     nEvents_( pset.getParameter<int>("NumberOfEvents") ),
//     weight_(0.),
    
//     if ( eventWeight_ ) {
//       double lumi = xSec_ > 0. ? ( nEvents_ / xSec_ ) : -1.;
//       weight_ = lumi > 0. ? ( normLumi_ / lumi ) : -1.;
//       if ( weight_ < 0. ) {
// 	edm::LogWarning("TEST") 
// 	  << " Problem with calculating lumi! "
// 	  << " CrossSection [pb]: " << xSec_
// 	  << " NumberOfEvents: " << nEvents_
// 	  << " Luminosity [pb-1]: " << lumi
// 	  << " NormalisedLumi [pb-1]: " << normLumi_
// 	  << " EventWeight: " << weight_;
// 	return;
//       } 
//     } else { weight_ = 1.; }
    
//     double sign = weight_ < 0. ? -1. : 1.;
//     double base = floor( log10( fabs(weight_) ) );
//     double coeff = fabs(weight_) / pow(10,base);
//     histo("EventWeight_Coefficient")->Fill( coeff * sign ); 
//     histo("EventWeight_Base")->Fill( base ); 
