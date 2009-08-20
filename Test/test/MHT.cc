#include "MHT.h"
#include "CommonTools/Utils/interface/EtComparator.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
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
MHT::MHT( const edm::ParameterSet& pset ) : 
  // Primary Objects
  jets_( pset.getParameter<edm::InputTag>("Jets") ),
  muons_( pset.getParameter<edm::InputTag>("Muons") ),
  electrons_( pset.getParameter<edm::InputTag>("Electrons") ),
  photons_( pset.getParameter<edm::InputTag>("Photons") ),
  // MET 
  met_( pset.getParameter<edm::InputTag>("MET") ),
  ccMet_( pset.getParameter<edm::InputTag>("CCMET") ),
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
  minPhotons_( pset.getParameter<int>("MinPhotons") ),
  // Histogram containers
  histos_(),
  histos2d_()
{
  //produces< std::vector<Candidate> >(tag).setBranchAlias(tag);
}


// -----------------------------------------------------------------------------
//
void MHT::analyze( const edm::Event& iEvent, 
		   const edm::EventSetup& iSetup ) {

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
  
  // -------------------- "Common objects" --------------------
    
  std::vector<Candidate> objects;
  if ( !photons.empty() ) { objects.insert( objects.end(), photons.begin(), photons.end() ); }
  if ( !jets.empty() ) { objects.insert( objects.end(), jets.begin(), jets.end() ); }
  if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Objects: " << objects.size(); }
  sort( objects.begin(), objects.end(), GreaterByEt<Candidate>() );

  // -------------------- Minimum number of objects --------------------
  
  if ( objects.size() < minObjects_ ) { return; }
  if ( jets.size() < minJets_ ) { return; }
  if ( muons.size() < minMuons_ ) { return; }
  if ( electrons.size() < minElectrons_ ) { return; }
  if ( photons.size() < minPhotons_ ) { return; }
    
  // -------------------- HT cut --------------------

  if ( ht( objects ) < totalHt_ ) { return; }

//   // -------------------- Get MET --------------------

//   double primary_met = -1.;
//   LorentzVector lv_primary_met;
//   std::vector<Candidate>::const_iterator iobj = objects.begin();
//   std::vector<Candidate>::const_iterator jobj = objects.end();
//   for ( ; iobj != jobj; ++iobj ) { lv_primary_met += iobj->p4(); }
//   lv_primary_met.SetPx( -1.*lv_primary_met.Px() );
//   lv_primary_met.SetPy( -1.*lv_primary_met.Py() );
//   lv_primary_met.SetPz( -1.*lv_primary_met.Pz() );
//   primary_met = sqrt( lv_primary_met.Perp2() );
    
//   double calo_met = -1.;
//   LorentzVector lv_calo_met;
//   if ( !met_.label().empty() ) {
//     edm::Handle< std::vector<pat::MET> > handle;
//     iEvent.getByLabel(met_,handle);
//     if ( !handle.isValid() ) { 
//       edm::LogWarning("TEST") << "No MET collection!";
//       return;
//     }
//     if ( handle->empty() ) { 
//       edm::LogWarning("TEST") << "Empty MET collection!";
//       return;
//     }
//     if ( !handle->front().isCaloMET() ) {
//       edm::LogWarning("TEST") << "Not CaloMET!" << endl;
//     }
//     calo_met = sqrt( handle->front().p4().Perp2() );
//     lv_calo_met = handle->front().p4();
//   }

//   double cc_met = -1.;
//   LorentzVector lv_cc_met;
//   if ( !ccMet_.label().empty() ) {
//     edm::Handle< std::vector<pat::MET> > handle;
//     iEvent.getByLabel(ccMet_,handle);
//     if ( !handle.isValid() ) { 
//       edm::LogWarning("TEST") << "No CCMET collection!";
//       return;
//     }
//     if ( handle->empty() ) { 
//       edm::LogWarning("TEST") << "Empty CCMET collection!";
//       return;
//     }
//     if ( !handle->front().isCaloMET() ) {
//       edm::LogWarning("TEST") << "Not CaloMET!" << endl;
//     }
//     cc_met = sqrt( handle->front().p4().Perp2() );
//     lv_cc_met = handle->front().p4();
//   }

//   double gen_met = -1.;
//   LorentzVector lv_gen_met;
//   if ( !genMet_.label().empty() ) {
//     edm::Handle< std::vector<reco::GenMET> > handle;
//     iEvent.getByLabel(genMet_,handle);
//     if ( !handle.isValid() ) { 
//       edm::LogWarning("TEST") << "No GenMET collection!";
//       return;
//     }
//     if ( handle->empty() ) { 
//       edm::LogWarning("TEST") << "Empty GenMET collection!";
//       return;
//     }
//     gen_met = sqrt( handle->front().p4().Perp2() );
//     lv_gen_met = handle->front().p4();
//   }
    
//   cout << " calo_met " << calo_met
//        << " cc_met " << cc_met
//        << " gen_met " << gen_met
//        << endl;

//   // -------------------- MET -------------------- 

//   {

//     double ht = ht( objects ); 

//     histo("PrimaryMET")->Fill( primary_met ); 
//     histo("CaloMET")->Fill( calo_met ); 
//     histo("CcMET")->Fill( cc_met ); 
//     histo("GenMET")->Fill( gen_met ); 
      
//     {
//       double delta_phi = reco::deltaPhi<LorentzVector,LorentzVector>( lv_primary_met, lv_gen_met );
//       histo("DPHI_PrimaryMET_GenMET")->Fill( delta_phi ); 
//     }

//     {
//       double delta_phi = reco::deltaPhi<LorentzVector,LorentzVector>( lv_calo_met, lv_gen_met );
//       histo("DPHI_CaloMET_GenMET")->Fill( delta_phi ); 
//     }

//     {
//       double delta_phi = reco::deltaPhi<LorentzVector,LorentzVector>( lv_cc_met, lv_gen_met );
//       histo("DPHI_CcMET_GenMET")->Fill( delta_phi ); 
//     }
      
//     histo2d("PrimaryMET_Vs_NObjects")->Fill( objects.size(), primary_met ); 
//     histo2d("CaloMET_Vs_NObjects")->Fill( objects.size(), calo_met ); 
//     histo2d("CcMET_Vs_NObjects")->Fill( objects.size(), cc_met ); 
//     histo2d("GenMET_Vs_NObjects")->Fill( objects.size(), gen_met ); 

//     histo2d("PrimaryMET_Vs_NJets")->Fill( jets.size(), primary_met ); 
//     histo2d("CaloMET_Vs_NJets")->Fill( jets.size(), calo_met ); 
//     histo2d("CcMET_Vs_NJets")->Fill( jets.size(), cc_met ); 
//     histo2d("GenMET_Vs_NJets")->Fill( jets.size(), gen_met ); 
      
//     histo2d("PrimaryMET_Vs_HT")->Fill( ht, primary_met ); 
//     histo2d("CaloMET_Vs_HT")->Fill( ht, calo_met ); 
//     histo2d("CcMET_Vs_HT")->Fill( ht, cc_met ); 
//     histo2d("GenMET_Vs_HT")->Fill( ht, gen_met ); 
      
//     histo2d("PrimaryMET_Vs_GenMET")->Fill( gen_met, primary_met ); 
//     histo2d("CaloMET_Vs_GenMET")->Fill( gen_met, calo_met ); 
//     histo2d("CcMET_Vs_GenMET")->Fill( gen_met, cc_met ); 
//     histo2d("PrimaryMET_Vs_CaloMET")->Fill( calo_met, primary_met ); 
      
//     if ( gen_met <= 0. ) { gen_met = 1.e-6; }
//     if ( calo_met <= 0. ) { calo_met = 1.e-6; }

//     histo2d("PrimaryMET/GenMET_Vs_NObjects")->Fill( objects.size(), primary_met/gen_met ); 
//     histo2d("CaloMET/GenMET_Vs_NObjects")->Fill( objects.size(), calo_met/gen_met ); 
//     histo2d("CcMET/GenMET_Vs_NObjects")->Fill( objects.size(), cc_met/gen_met ); 
//     histo2d("PrimaryMET/CaloMET_Vs_NObjects")->Fill( objects.size(), primary_met/calo_met ); 
      
//     histo2d("PrimaryMET/GenMET_Vs_NJets")->Fill( jets.size(), primary_met/gen_met ); 
//     histo2d("CaloMET/GenMET_Vs_NJets")->Fill( jets.size(), calo_met/gen_met ); 
//     histo2d("CcMET/GenMET_Vs_NJets")->Fill( jets.size(), cc_met/gen_met ); 
//     histo2d("PrimaryMET/CaloMET_Vs_NJets")->Fill( jets.size(), primary_met/calo_met ); 
      
//     histo2d("PrimaryMET/GenMET_Vs_HT")->Fill( ht, primary_met/gen_met ); 
//     histo2d("CaloMET/GenMET_Vs_HT")->Fill( ht, calo_met/gen_met ); 
//     histo2d("CcMET/GenMET_Vs_HT")->Fill( ht, cc_met/gen_met ); 
//     histo2d("PrimaryMET/CaloMET_Vs_HT")->Fill( ht, primary_met/calo_met ); 
      
//     histo("PrimaryMET-GenMET/GenMET")->Fill( (primary_met-gen_met)/gen_met ); 
//     histo("CaloMET-GenMET/GenMET")->Fill( (calo_met-gen_met)/gen_met ); 
//     histo("CcMET-GenMET/GenMET")->Fill( (cc_met-gen_met)/gen_met ); 
//     histo("PrimaryMET-CaloMET/CaloMET")->Fill( (primary_met-calo_met)/calo_met ); 
      
//     histo2d("PrimaryMET-GenMET/GenMET_Vs_NObjects")->Fill( objects.size(), (primary_met-gen_met)/gen_met ); 
//     histo2d("CaloMET-GenMET/GenMET_Vs_NObjects")->Fill( objects.size(), (calo_met-gen_met)/gen_met ); 
//     histo2d("CcMET-GenMET/GenMET_Vs_NObjects")->Fill( objects.size(), (cc_met-gen_met)/gen_met ); 
//     histo2d("PrimaryMET-CaloMET/CaloMET_Vs_NObjects")->Fill( objects.size(), (primary_met-calo_met)/calo_met ); 

//     histo2d("PrimaryMET-GenMET/GenMET_Vs_HT")->Fill( ht, (primary_met-gen_met)/gen_met ); 
//     histo2d("CaloMET-GenMET/GenMET_Vs_HT")->Fill( ht, (calo_met-gen_met)/gen_met ); 
//     histo2d("CcMET-GenMET/GenMET_Vs_HT")->Fill( ht, (cc_met-gen_met)/gen_met ); 
//     histo2d("PrimaryMET-CaloMET/CaloMET_Vs_HT")->Fill( ht, (primary_met-calo_met)/calo_met ); 

//     histo2d("PrimaryMET-GenMET/GenMET_Vs_NObjects_Vs_HT")->Fill( objects.size(), ht, (primary_met-gen_met)/gen_met ); 
//     histo2d("CaloMET-GenMET/GenMET_Vs_NObjects_Vs_HT")->Fill( objects.size(), ht, (calo_met-gen_met)/gen_met ); 
//     histo2d("CcMET-GenMET/GenMET_Vs_NObjects_Vs_HT")->Fill( objects.size(), ht, (cc_met-gen_met)/gen_met ); 
//     histo2d("PrimaryMET-CaloMET/CaloMET_Vs_NObjects_Vs_HT")->Fill( objects.size(), ht, (primary_met-calo_met)/calo_met ); 
      
//   }

//   // -------------------- MHT vs HT --------------------

//   {
      
//     double ht_pho = ht( photons ); 
//     double ht_jet = ht( jets ); 
//     double ht_all = ht( objects ); 
//     double mht = mht( objects );
      
//     histo2d("MHT_Vs_PhotonHT")->Fill( ht_pho, mht ); 
//     histo2d("MHT_Vs_JetHT")->Fill( ht_jet, mht ); 
//     histo2d("MHT_Vs_HT")->Fill( ht_all, mht ); 

//   }

}    

// -----------------------------------------------------------------------------
//
void MHT::beginJob( const edm::EventSetup& ) {

//   edm::Service<TFileService> fs;

//   // -------------------- "MET" --------------------
  
//   { 
      
//     TFileDirectory dir = fs->mkdir("MET");
      
//     histos_["PrimaryMET"] = dir.make<TH1D>("PrimaryMET","",100,0.,2000.);
//     histos_["CaloMET"] = dir.make<TH1D>("CaloMET","",100,0.,2000.);
//     histos_["CcMET"] = dir.make<TH1D>("CcMET","",100,0.,2000.);
//     histos_["GenMET"] = dir.make<TH1D>("GenMET","",100,0.,2000.);

//     histos_["DPHI_PrimaryMET_GenMET"] = dir.make<TH1D>("DPHI_PrimaryMET_GenMET","",160,-3.2,3.2);
//     histos_["DPHI_CaloMET_GenMET"] = dir.make<TH1D>("DPHI_CaloMET_GenMET","",160,-3.2,3.2);
//     histos_["DPHI_CcMET_GenMET"] = dir.make<TH1D>("DPHI_CcMET_GenMET","",160,-3.2,3.2);

//     histos2d_["PrimaryMET_Vs_NObjects"] = dir.make<TH2D>("PrimaryMET_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);
//     histos2d_["CaloMET_Vs_NObjects"] = dir.make<TH2D>("CaloMET_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);
//     histos2d_["CcMET_Vs_NObjects"] = dir.make<TH2D>("CcMET_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);
//     histos2d_["GenMET_Vs_NObjects"] = dir.make<TH2D>("GenMET_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);

//     histos2d_["PrimaryMET_Vs_NJets"] = dir.make<TH2D>("PrimaryMET_Vs_NJets","",51,-0.5,50.5,100,0.,2000.);
//     histos2d_["CaloMET_Vs_NJets"] = dir.make<TH2D>("CaloMET_Vs_NJets","",51,-0.5,50.5,100,0.,2000.);
//     histos2d_["CcMET_Vs_NJets"] = dir.make<TH2D>("CcMET_Vs_NJets","",51,-0.5,50.5,100,0.,2000.);
//     histos2d_["GenMET_Vs_NJets"] = dir.make<TH2D>("GenMET_Vs_NJets","",51,-0.5,50.5,100,0.,2000.);
      
//     histos2d_["PrimaryMET_Vs_HT"] = dir.make<TH2D>("PrimaryMET_Vs_HT","",100,0.,2000.,100,0.,2000.);
//     histos2d_["CaloMET_Vs_HT"] = dir.make<TH2D>("CaloMET_Vs_HT","",100,0.,2000.,100,0.,2000.);
//     histos2d_["CcMET_Vs_HT"] = dir.make<TH2D>("CcMET_Vs_HT","",100,0.,2000.,100,0.,2000.);
//     histos2d_["GenMET_Vs_HT"] = dir.make<TH2D>("GenMET_Vs_HT","",100,0.,2000.,100,0.,2000.);
      
//     histos2d_["PrimaryMET_Vs_GenMET"] = dir.make<TH2D>("PrimaryMET_Vs_GenMET","",100,0.,2000.,100,0.,2000.);
//     histos2d_["CaloMET_Vs_GenMET"] = dir.make<TH2D>("CaloMET_Vs_GenMET","",100,0.,2000.,100,0.,2000.);
//     histos2d_["CcMET_Vs_GenMET"] = dir.make<TH2D>("CcMET_Vs_GenMET","",100,0.,2000.,100,0.,2000.);
//     histos2d_["PrimaryMET_Vs_CaloMET"] = dir.make<TH2D>("PrimaryMET_Vs_CaloMET","",100,0.,2000.,100,0.,2000.);
      
//     histos2d_["PrimaryMET/GenMET_Vs_NObjects"] = dir.make<TH2D>("PrimaryMET/GenMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);
//     histos2d_["CaloMET/GenMET_Vs_NObjects"] = dir.make<TH2D>("CaloMET/GenMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);
//     histos2d_["CcMET/GenMET_Vs_NObjects"] = dir.make<TH2D>("CcMET/GenMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);
//     histos2d_["PrimaryMET/CaloMET_Vs_NObjects"] = dir.make<TH2D>("PrimaryMET/CaloMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);

//     histos2d_["PrimaryMET/GenMET_Vs_NJets"] = dir.make<TH2D>("PrimaryMET/GenMET_Vs_NJets","",51,-0.5,50.5,200,-10.,10.);
//     histos2d_["CaloMET/GenMET_Vs_NJets"] = dir.make<TH2D>("CaloMET/GenMET_Vs_NJets","",51,-0.5,50.5,200,-10.,10.);
//     histos2d_["CcMET/GenMET_Vs_NJets"] = dir.make<TH2D>("CcMET/GenMET_Vs_NJets","",51,-0.5,50.5,200,-10.,10.);
//     histos2d_["PrimaryMET/CaloMET_Vs_NJets"] = dir.make<TH2D>("PrimaryMET/CaloMET_Vs_NJets","",51,-0.5,50.5,200,-10.,10.);
      
//     histos2d_["PrimaryMET/GenMET_Vs_HT"] = dir.make<TH2D>("PrimaryMET/GenMET_Vs_HT","",100,0.,2000.,200,-10.,10.);
//     histos2d_["CaloMET/GenMET_Vs_HT"] = dir.make<TH2D>("CaloMET/GenMET_Vs_HT","",100,0.,2000.,200,-10.,10.);
//     histos2d_["CcMET/GenMET_Vs_HT"] = dir.make<TH2D>("CcMET/GenMET_Vs_HT","",100,0.,2000.,200,-10.,10.);
//     histos2d_["PrimaryMET/CaloMET_Vs_HT"] = dir.make<TH2D>("PrimaryMET/CaloMET_Vs_HT","",100,0.,2000.,200,-10.,10.);

//     histos_["PrimaryMET-GenMET/GenMET"] = dir.make<TH1D>("PrimaryMET-GenMET/GenMET","",51,-0.5,50.5);
//     histos_["CaloMET-GenMET/GenMET"] = dir.make<TH1D>("CaloMET-GenMET/GenMET","",51,-0.5,50.5);
//     histos_["CcMET-GenMET/GenMET"] = dir.make<TH1D>("CcMET-GenMET/GenMET","",51,-0.5,50.5);
//     histos_["PrimaryMET-CaloMET/CaloMET"] = dir.make<TH1D>("PrimaryMET-CaloMET/CaloMET","",51,-0.5,50.5);

//     histos2d_["PrimaryMET-GenMET/GenMET_Vs_NObjects"] = 
//       dir.make<TH2D>("PrimaryMET-GenMET/GenMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);
//     histos2d_["CaloMET-GenMET/GenMET_Vs_NObjects"] = 
//       dir.make<TH2D>("CaloMET-GenMET/GenMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);
//     histos2d_["CcMET-GenMET/GenMET_Vs_NObjects"] = 
//       dir.make<TH2D>("CcMET-GenMET/GenMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);
//     histos2d_["PrimaryMET-CaloMET/CaloMET_Vs_NObjects"] = 
//       dir.make<TH2D>("PrimaryMET-CaloMET/CaloMET_Vs_NObjects","",51,-0.5,50.5,200,-10.,10.);
      
//     histos2d_["PrimaryMET-GenMET/GenMET_Vs_HT"] = 
//       dir.make<TH2D>("PrimaryMET-GenMET/GenMET_Vs_HT","",100,0.,2000.,200,-10.,10.);
//     histos2d_["CaloMET-GenMET/GenMET_Vs_HT"] = 
//       dir.make<TH2D>("CaloMET-GenMET/GenMET_Vs_HT","",100,0.,2000.,200,-10.,10.);
//     histos2d_["CcMET-GenMET/GenMET_Vs_HT"] = 
//       dir.make<TH2D>("CcMET-GenMET/GenMET_Vs_HT","",100,0.,2000.,200,-10.,10.);
//     histos2d_["PrimaryMET-CaloMET/CaloMET_Vs_HT"] = 
//       dir.make<TH2D>("PrimaryMET-CaloMET/CaloMET_Vs_HT","",100,0.,2000.,200,-10.,10.);

//     histos2d_["PrimaryMET-GenMET/GenMET_Vs_NObjects_Vs_HT"] = 
//       dir.make<TH2D>("PrimaryMET-GenMET/GenMET_Vs_NObjects_Vs_HT","",51,-0.5,50.5,100,0.,2000.);
//     histos2d_["CaloMET-GenMET/GenMET_Vs_NObjects_Vs_HT"] = 
//       dir.make<TH2D>("CaloMET-GenMET/GenMET_Vs_NObjects_Vs_HT","",51,-0.5,50.5,100,0.,2000.);
//     histos2d_["CcMET-GenMET/GenMET_Vs_NObjects_Vs_HT"] = 
//       dir.make<TH2D>("CcMET-GenMET/GenMET_Vs_NObjects_Vs_HT","",51,-0.5,50.5,100,0.,2000.);
//     histos2d_["PrimaryMET-CaloMET/CaloMET_Vs_NObjects_Vs_HT"] = 
//       dir.make<TH2D>("PrimaryMET-CaloMET/CaloMET_Vs_NObjects_Vs_HT","",51,-0.5,50.5,100,0.,2000.);

//   }      

//   // -------------------- "Common objects" --------------------
  
//   { 
      
//     TFileDirectory dir = fs->mkdir("CommonObjects");
      
//     histos_["PreHtCut_NPhotons"] = dir.make<TH1D>("PreHtCut_NPhotons","",51,-0.5,50.5);
//     histos_["PreHtCut_NJets"] = dir.make<TH1D>("PreHtCut_NJets","",51,-0.5,50.5);
//     histos_["PreHtCut_NMuons"] = dir.make<TH1D>("PreHtCut_NMuons","",51,-0.5,50.5);
//     histos_["PreHtCut_NElectrons"] = dir.make<TH1D>("PreHtCut_NElectrons","",51,-0.5,50.5);
      
//     histos_["PostHtCut_NPhotons"] = dir.make<TH1D>("PostHtCut_NPhotons","",51,-0.5,50.5);
//     histos_["PostHtCut_NJets"] = dir.make<TH1D>("PostHtCut_NJets","",51,-0.5,50.5);
//     histos_["PostHtCut_NObjects"] = dir.make<TH1D>("PostHtCut_NObjects","",51,-0.5,50.5);
      
//     histos_["GenPhotons_Et"] = dir.make<TH1D>("GenPhotons_Et","",200,0.,1000.);
      
//   }

//   // -------------------- Event weight --------------------

//   { 

//     TFileDirectory dir = fs->mkdir("Common");

//     histos_["CommonObjectsHT"] = dir.make<TH1D>("CommonObjectsHT","",100,0.,2000.);
//     histos2d_["CommonObjectsHT_Vs_NObjects"] = dir.make<TH2D>("CommonObjectsHT_Vs_NObjects","",51,-0.5,50.5,100,0.,2000.);

//     histos_["CutFlow_Efficiency"] = dir.make<TH1D>("CutFlow_Efficiency","",11,-0.5,10.5);
//     //histos_["EventWeight_Coefficient"] = dir.make<TH1D>("EventWeight_Coefficient","",1000,0.,10.);
//     //histos_["EventWeight_Base"] = dir.make<TH1D>("EventWeight_Base","",21,-10.5,10.5);
    
//   }

}

// -----------------------------------------------------------------------------
//
bool MHT::getJets( const edm::Event& iEvent,
		   std::vector<Candidate>& jets ) {
  
  jets.clear();
  
  if ( !jets_.label().empty() ) {

    edm::Handle< std::vector<pat::Jet> > handle;
    iEvent.getByLabel(jets_,handle);
    
    if ( !handle.isValid() ) { 
      edm::LogWarning("TEST") << "No Jets for " << jets_; 
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

  sort( jets.begin(), jets.end(), GreaterByEt<Candidate>() );
  
  return false;
  
}

// -----------------------------------------------------------------------------
//
bool MHT::getMuons( const edm::Event& iEvent,
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
	   imuon->trackIso() < muonTrkIso_ &&
	   imuon->pt() > muonPt_ ) { muons.push_back( *imuon ); }
    }
    
    if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Muons: " << muons.size(); }

  }

  sort( muons.begin(), muons.end(), GreaterByEt<Candidate>() );

  return false;
  
}

// -----------------------------------------------------------------------------
//
bool MHT::getElectrons( const edm::Event& iEvent,
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
	   ielectron->trackIso() < electronTrkIso_ &&
	   ielectron->pt() > electronPt_ ) { electrons.push_back( *ielectron ); }
    }
    
    if ( edm::isDebugEnabled() ) { edm::LogVerbatim("TEST") << "Number of Electrons: " << electrons.size(); }
    
  }

  sort( electrons.begin(), electrons.end(), GreaterByEt<Candidate>() );
  
  return false;
  
}

// -----------------------------------------------------------------------------
//
bool MHT::getPhotons( const edm::Event& iEvent,
		      std::vector<Candidate>& photons ) {
  
  photons.clear();
  
  if ( !photons_.label().empty() ) {
    
    edm::Handle< std::vector<pat::Photon> > handle;
    iEvent.getByLabel(photons_,handle);
    
    if ( !handle.isValid() ) { 
      edm::LogWarning("TEST") << "No Photons for " << photons_; 
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

  sort( photons.begin(), photons.end(), GreaterByEt<Candidate>() );
  
  return false; 

}

// -----------------------------------------------------------------------------
//
double MHT::ht( const std::vector<Candidate>& input ) {
  LorentzVector lv;
  std::vector<Candidate>::const_iterator ii = input.begin();
  std::vector<Candidate>::const_iterator jj = input.end();
  for ( ; ii != jj; ++ii ) { lv += ii->p4(); }
  double ht = lv.E()*lv.E() - lv.Pz()*lv.Pz();
  ht = ht < 0. ? -1.*sqrt(-1.*ht) : sqrt(ht);
  return ht;
}

// -----------------------------------------------------------------------------
//
double MHT::mht( const std::vector<Candidate>& input ) {
  LorentzVector lv_obj;
  std::vector<Candidate>::const_iterator iobj = input.begin();
  std::vector<Candidate>::const_iterator jobj = input.end();
  for ( ; iobj != jobj; ++iobj ) { lv_obj += iobj->p4(); }
  LorentzVector lv_mht( lv_obj );
  lv_mht.SetPx( -1.*lv_obj.Px() );
  lv_mht.SetPy( -1.*lv_obj.Py() );
  lv_mht.SetPz( -1.*lv_obj.Pz() );
  return sqrt( lv_mht.Perp2() );
}

// -----------------------------------------------------------------------------
//
TH1D* MHT::histo( const std::string& histogram_name ) {
  std::map<std::string,TH1D*>::const_iterator ii = histos_.find(histogram_name);
  if ( ii != histos_.end() ) { return ii->second; }
  edm::LogWarning("TEST") << "Cannot find string: " << histogram_name;
  return 0;
}
  
// -----------------------------------------------------------------------------
//
TH2D* MHT::histo2d( const std::string& histogram_name ) {
  std::map<std::string,TH2D*>::const_iterator jj = histos2d_.find(histogram_name);
  if ( jj != histos2d_.end() ) { return jj->second; }
  edm::LogWarning("TEST") << "Cannot find string: " << histogram_name;
  return 0;
}

// -----------------------------------------------------------------------------
// 
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MHT);
