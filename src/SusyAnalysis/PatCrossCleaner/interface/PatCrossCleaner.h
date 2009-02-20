// -*- C++ -*-
//
// Package:    PatCrossCleaner
// Class:      PatCrossCleaner
// 
/**\class PatCrossCleaner PatCrossCleaner.h SusyAnalysis/PatCrossCleaner/interface/PatCrossCleaner.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Christian AUTERMANN
//         Created:  Sat Mar 22 12:58:04 CET 2008
// $Id: PatCrossCleaner.h,v 1.1 2009/02/20 12:42:08 bainbrid Exp $
//
//

#ifndef PatCrossCleaner_h
#define PatCrossCleaner_h

// system include files
#include <memory>
#include <map>
#include <utility>//pair
#include <cmath>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/UtilAlgos/interface/ParameterAdapter.h"
#include "SusyAnalysis/PatCrossCleaner/interface/ElectronJetCrossCleaner.h"
#include "SusyAnalysis/PatCrossCleaner/interface/PhotonJetCrossCleaner.h"
#include "SusyAnalysis/PatCrossCleaner/interface/MuonJetCrossCleaner.h"
#include "SusyAnalysis/PatCrossCleaner/interface/ElectronPhotonCrossCleaner.h"
#include "SusyAnalysis/PatCrossCleaner/interface/CrossCleanerResult.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"

//
// class decleration
//

class PatCrossCleaner : public edm::EDProducer {
   public:
      explicit PatCrossCleaner(const edm::ParameterSet&);
      ~PatCrossCleaner();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
      /// Input: All PAT objects that are to cross-clean  or needed for that
      edm::InputTag _patJets;
      edm::InputTag _patMets;
      edm::InputTag _patMuons;
      edm::InputTag _patElectrons;
      edm::InputTag _patPhotons;
      edm::InputTag _patTaus;
      // The following are no PAT objects (yet?)
      edm::InputTag _patCaloTowers;
      edm::InputTag _patTracks;
      edm::InputTag _patVertices;

      // Service for Jet corrections after cleaning
      std::string L1JetCorrService_;
      std::string L2JetCorrService_;
      std::string L3JetCorrService_;
      std::string L4JetCorrService_;
      std::string L6JetCorrService_;
      std::string L5udsJetCorrService_;
      std::string L5gluJetCorrService_;
      std::string L5cJetCorrService_;
      std::string L5bJetCorrService_;
      std::string L7udsJetCorrService_;
      std::string L7gluJetCorrService_;
      std::string L7cJetCorrService_;
      std::string L7bJetCorrService_;
      
      bool bl1_,bl2_,bl3_,bl4_,bl6_,bl5uds_,bl5g_,bl5c_,bl5b_,bl7uds_,bl7g_,bl7c_,bl7b_;
      
      ///The actual cross-cleaners:
      //Electron-Jet
      bool _doElectronJetCC;
      edm::ParameterSet _EJselectionCfg;
      pat::ElectronJetCrossCleaner _ElectronJetCC;
      
      //Photon-Jet
      bool _doPhotonJetCC;
      edm::ParameterSet _PJselectionCfg;
      pat::PhotonJetCrossCleaner _PhotonJetCC;

      //Muon-Jet
      bool _doMuonJetCC;
      edm::ParameterSet _MJselectionCfg;
      pat::MuonJetCrossCleaner _MuonJetCC;

      //Electron-Photon
      bool _doElectronPhotonCC;
      pat::ElectronPhotonCrossCleaner _ElectronPhotonCC;

  template <class Object> void putObjects(CrossCleanerResult& , const edm::View<Object>& initial,edm::Event&, const edm::EventSetup&, std::string);

  // Cross-cleaning corrections
  template <class Object> void correctObject( Object&, CrossCleanerModifier&, const edm::EventSetup& );
  void correctObject( pat::Electron&, CrossCleanerModifier&, const edm::EventSetup& );
  void correctObject( pat::Jet&, CrossCleanerModifier&, const edm::EventSetup& );
  
  // Jet energy corrections
  template <class Object> void applyJEC( Object&, const edm::Event&, const edm::EventSetup& );
  void applyJEC( pat::Jet&, const edm::Event&, const edm::EventSetup& );

  // MET corrections
  template <class Object> CrossCleanerMETCorrection metCorr( Object&, const Object& );
  CrossCleanerMETCorrection metCorr( pat::Jet&, const pat::Jet& );
  
  template <class Object> CrossCleanerMETCorrection droppingObject(const Object & o);
  CrossCleanerMETCorrection droppingObject(const pat::Jet & j);

      ///Helper function for debugging: print out all objects that should be dropped/corrected
      void printConflicts(CrossCleanerMap& conflicts);
      void printCorrected(std::map<edm::RefToBase<reco::Candidate>, double>& conflicts) const;
      void printDropped(std::set<edm::RefToBase<reco::Candidate> >) const;
      ///statistical information
      std::map<int,std::pair<int,int> > _statistics;
  std::string ccinstance_;
  std::string droppedinstance_;
};

//two methods of sorting
template <class Object> bool sortObjectByEt(const Object & lhs,const Object & rhs){ return (lhs.et()>rhs.et());}
template <class Object> bool sortObjectByPt(const Object & lhs,const Object & rhs){ return (lhs.pt()>rhs.pt());}
//generic object are sorted by Pt
template <class Object> bool sortObject(const Object & lhs,const Object & rhs){ return sortObjectByPt(lhs,rhs);}
//specification: pat::Jet are sorted by Et
template<> bool sortObject<pat::Jet>(const pat::Jet & j_lhs, const pat::Jet & j_rhs){return sortObjectByEt(j_lhs,j_rhs);}

template <class Object> 
void PatCrossCleaner::correctObject(Object & o, CrossCleanerModifier & modifier, const edm::EventSetup&){
  //@@ do nothing
}

void PatCrossCleaner::correctObject(pat::Jet & o, CrossCleanerModifier & modifier, const edm::EventSetup& es){

  if (modifier.keepKeyObj){
    if (abs(modifier.object->pdgId()==13)){
      LogDebug("PatCrossCleaner")<<"I am correcting a jet of momentum: "<<o.momentum()
				 <<" because of a muon of momentum: "<<modifier.object->momentum();
    }
      
    if (o.energy()!=0){
      LogDebug("PatCrossCleaner") << "Original Jet p4 " << o.p4()
                                  << "\n jet correction is " << o.jetCorrFactors().scaleDefault()
                                  << "\n Original uncorrected Jet p4" << o.correctedJet("RAW").p4();
      
      const pat::Jet& oNoCorr = o.correctedJet("RAW"); // change uncorrected jet
      double theSharedE = sqrt( modifier.modEnergy.Mag2() );
      
      pat::Jet::LorentzVector newP4( ( oNoCorr.px() - modifier.modEnergy.X() ), 
				     ( oNoCorr.py() - modifier.modEnergy.Y() ), 
				     ( oNoCorr.pz() - modifier.modEnergy.Z() ),
				     ( oNoCorr.energy() - theSharedE ) );
      o.setP4( newP4 ); 
      
      LogDebug("PatCrossCleaner")<< "changing jet uncorrected energy from: " << oNoCorr.energy()
				 <<" to: " << ( oNoCorr.energy() - theSharedE )
				 << "\n the modifier energy is: "<< theSharedE
				 << "\n modifier id: "<< modifier.object->pdgId();
      
    }
  }
}

// correct an electron due to a modifier.
void PatCrossCleaner::correctObject(pat::Electron & o, CrossCleanerModifier & modifier, const edm::EventSetup&){
  if (modifier.keepKeyObj){
    if (o.energy()!=0){
      double oldEt = o.et();
      double newE = o.energy() + sqrt( modifier.modEnergy.Mag2() );
      double k = newE/o.energy();
      pat::Electron::LorentzVector newP4(k*o.px(), k*o.py(), k*o.pz(), newE );
      o.setP4( newP4 );
      LogDebug("PatCrossCleaner")<<"changing an electron energy from: "<<oldEt<<" to: "<<o.et()
				 <<"\n the k factor is: "<<k
				 <<"\n the modifier energy is: "<< sqrt( modifier.modEnergy.Mag2() )
				 <<"\n modifier id: "<<modifier.object->pdgId();
    }
  }
}

template <class Object> 
CrossCleanerMETCorrection PatCrossCleaner::droppingObject(const Object & o) { return CrossCleanerMETCorrection();}

CrossCleanerMETCorrection PatCrossCleaner::droppingObject(const pat::Jet & j){
  return CrossCleanerMETCorrection(j.correctedJet("RAW").px() - j.px(), j.correctedJet("RAW").py() - j.py());
}


template <class Object>
void PatCrossCleaner::putObjects(CrossCleanerResult& result, const edm::View<Object> & initial,
                                edm::Event & iEvent, const edm::EventSetup& es, std::string Label)
{
  CrossCleanerMap & assMap = result.map;
  std::auto_ptr<std::vector<Object> >  cleanObj(new std::vector<Object>);
  std::auto_ptr<std::vector<Object> >  droppedObj(new std::vector<Object>);
  
  cleanObj->reserve(initial.size()); //reserve enough memory there
  droppedObj->reserve(initial.size()); //reserve enough memory there
  
  CrossCleanerMETCorrection & totalMETcorr = result.metCorrection;

  for ( uint i = 0; i < initial.size(); ++i ) {
    
    const Object& initialO = initial[i];
    cleanObj->push_back(initialO);
    droppedObj->push_back(initialO);
    Object& newObject = cleanObj->back();
    
    edm::RefToBase<reco::Candidate> candRef(initial.refAt(i) );
    
    bool keepMe=true;
    CrossCleanerMap::iterator action = assMap.find(candRef);
    if (action != assMap.end()) {
      if (action->second.keepFromModifiers()){
	for (uint m=0;m!=action->second.modifiers.size();++m){
	  correctObject( newObject, action->second.modifiers[m], es );
	}
      } else { keepMe = false; }
    }
    
    applyJEC( newObject, iEvent, es );
    if ( newObject.energy() < 0.0 ) { keepMe = false; }
    
    CrossCleanerMETCorrection thisMETcorr = metCorr( newObject, initialO );
    totalMETcorr+=thisMETcorr;
    
    if (!keepMe){
      //dropping the object
      cleanObj->pop_back();
      ++_statistics[candRef->pdgId()].second;
      totalMETcorr+=droppingObject(initialO);
    }
    else {
      //keeping the object
      droppedObj->pop_back();
      ++_statistics[candRef->pdgId()].first;
    }

  }

  LogDebug("PatCrossCleaner")<<" putting ("<< Label <<") collection in the event"
			     <<" \n keep: "<<cleanObj->size()
			     <<" \n drop: "<<droppedObj->size();

  //sort them
  sort(cleanObj->begin(), cleanObj->end(), sortObject<Object>);

  iEvent.put(cleanObj, ccinstance_+Label);
  iEvent.put(droppedObj, droppedinstance_+Label);
  //  return res;
}


namespace reco {
  namespace modules {
    /// Helper struct to convert from ParameterSet to ElectronSelection
    template<> 
    struct ParameterAdapter<pat::ElectronJetCrossCleaner> { 
      static pat::ElectronJetCrossCleaner make(const edm::ParameterSet & cfg) {
        pat::ElectronJetCleaning config_;
	config_.SusyAnalyzerCleaning = cfg.getParameter<bool>("SusyAnalyzerCleaning");
        config_.deltaR_min = cfg.getParameter<double>("deltaR_min");
        config_.SharedEtoJetE= cfg.getParameter<double>("SharedEtoJetE");
        config_.IsoValueCut= cfg.getParameter<double>("IsoValueCut");
        config_.SharedEForNIsoEle= cfg.getParameter<double>("SharedEForNIsoEle");
        config_.IsolationKey= cfg.getParameter<std::string>("IsolationKey");
	config_.ElectronID = cfg.getParameter<std::string>("ElectronID");
	//... read in further parameters for Electron-Jet here
	
	return pat::ElectronJetCrossCleaner( config_ );
      }
    };

    /// Helper struct to convert from ParameterSet to PhotonSelection
    template<> 
    struct ParameterAdapter<pat::PhotonJetCrossCleaner> { 
      static pat::PhotonJetCrossCleaner make(const edm::ParameterSet & cfg) {
        pat::PhotonJetCleaning config_;
        config_.deltaR_min = cfg.getParameter<double>("deltaR_min");
        config_.IsoValueCut= cfg.getParameter<double>("IsoValueCut");
        config_.IsolationKey= cfg.getParameter<std::string>("IsolationKey");
	config_.PhotonID = cfg.getParameter<std::string>("PhotonID");
	//... read in further parameters for Photon-Jet here
	
	return pat::PhotonJetCrossCleaner( config_ );
      }
    };
    
    /// Helper struct to convert from ParameterSet to MuonSelection
    template<> 
    struct ParameterAdapter<pat::MuonJetCrossCleaner> { 
      static pat::MuonJetCrossCleaner make(const edm::ParameterSet & cfg) {
        pat::MuonJetCleaning config_;
        config_.deltaR_min = cfg.getParameter<double>("deltaR_min");
	config_.caloIso_max = cfg.getParameter<double>("caloIso_max");
	config_.trackIso_max = cfg.getParameter<double>("trackIso_max");
	config_.muonID = cfg.getParameter<std::string>("MuonID");
	config_.modifyJetEnergy = cfg.getParameter<bool>("modifyJetEnergy");
	//... read in further parameters for Muon-Jet here
	
	return pat::MuonJetCrossCleaner( config_ );
      }
    };

  }
}

template <class Object> 
void PatCrossCleaner::applyJEC( Object& obj, const edm::Event& event, const edm::EventSetup& es ) {
  // nothing yet
}

void PatCrossCleaner::applyJEC( pat::Jet& o, const edm::Event& event, const edm::EventSetup& es ) {

//   std::cout << "JET0: " << o.jetCorrFactors().scaleDefault() << std::endl;
  
  // get jet correction services
  const JetCorrector 
    *L1JetCorr=0, *L2JetCorr=0, *L3JetCorr=0, *L4JetCorr=0, 
    *L5udsJetCorr=0,*L5gluJetCorr=0,*L5cJetCorr=0,*L5bJetCorr=0,
    *L6JetCorr=0,
    *L7udsJetCorr=0,*L7gluJetCorr=0,*L7cJetCorr=0,*L7bJetCorr=0;

  if (bl1_)    L1JetCorr     = JetCorrector::getJetCorrector(L1JetCorrService_, es);
  if (bl2_)    L2JetCorr     = JetCorrector::getJetCorrector(L2JetCorrService_, es);
  if (bl3_)    L3JetCorr     = JetCorrector::getJetCorrector(L3JetCorrService_, es);
  if (bl4_)    L4JetCorr     = JetCorrector::getJetCorrector(L4JetCorrService_, es);
  if (bl6_)    L6JetCorr     = JetCorrector::getJetCorrector(L6JetCorrService_, es);
  if (bl5uds_) L5udsJetCorr  = JetCorrector::getJetCorrector(L5udsJetCorrService_, es);
  if (bl5g_)   L5gluJetCorr  = JetCorrector::getJetCorrector(L5gluJetCorrService_, es);
  if (bl5c_)   L5cJetCorr    = JetCorrector::getJetCorrector(L5cJetCorrService_, es);
  if (bl5b_)   L5bJetCorr    = JetCorrector::getJetCorrector(L5bJetCorrService_, es);
  if (bl7uds_) L7udsJetCorr  = JetCorrector::getJetCorrector(L7udsJetCorrService_, es);
  if (bl7g_)   L7gluJetCorr  = JetCorrector::getJetCorrector(L7gluJetCorrService_, es);
  if (bl7c_)   L7cJetCorr    = JetCorrector::getJetCorrector(L7cJetCorrService_, es);
  if (bl7b_)   L7bJetCorr    = JetCorrector::getJetCorrector(L7bJetCorrService_, es);

  // retrieve the energy correction factors
  float l1=-1, l2=-1, l3=-1, l4=-1, l6=-1;
  pat::JetCorrFactors::FlavourCorrections l5, l7; 
      
  // start from RAW; then nest all subsequent correction step,
  // i.e. input is the 4-vector of jet correted up to level n-1
  reco::Jet refJet = (o);
  if (bl1_){
    l1 = L1JetCorr->correction( refJet, event, es );
    refJet.setP4(l1*refJet.p4());
  }
  if (bl2_){
    l2 = L2JetCorr->correction( refJet, event, es );  
    refJet.setP4(l2*refJet.p4());
  }
  if (bl3_){
    l3 = L3JetCorr->correction( refJet, event, es );
    refJet.setP4(l3*refJet.p4());
  }
  if (bl4_){
    l4 = L4JetCorr->correction( refJet );
    refJet.setP4(l4*refJet.p4());
  }
  if (bl6_){
    l6 = L6JetCorr->correction( refJet );
    refJet.setP4(l6*refJet.p4());
  }

  // from here on 4 variations of jet corrections co-exist, which 
  // start from the up to level 4 full corrected jet
  reco::Jet refJetuds = refJet;
  reco::Jet refJetg   = refJet;
  reco::Jet refJetc   = refJet;
  reco::Jet refJetb   = refJet;

  // ----------------------------------
  // L5 hadron level corrections
  // ----------------------------------
  if (bl5uds_){
    l5.uds = L5udsJetCorr->correction( refJetuds );
    refJetuds.setP4(l5.uds*refJet.p4());
  }
  if (bl5g_){
    l5.g   = L5gluJetCorr->correction( refJetg   );
    refJetg  .setP4(l5.g  *refJet.p4());
  }
  if (bl5c_){
    l5.c   = L5cJetCorr->correction  ( refJetc   );
    refJetc  .setP4(l5.c  *refJet.p4());
  }
  if (bl5b_){
    l5.b   = L5bJetCorr->correction  ( refJetb   );
    refJetb  .setP4(l5.b  *refJet.p4());
  }
  // ----------------------------------
  // L7 parton level corrections
  // ----------------------------------
  if (bl7uds_){
    l7.uds = L7udsJetCorr->correction( refJetuds );
    refJetuds.setP4(l7.uds*refJet.p4());
  }
  if (bl7g_){
    l7.g   = L7gluJetCorr->correction( refJetg   );
    refJetg  .setP4(l7.g  *refJet.p4());
  }
  if (bl7c_){
    l7.c   = L7cJetCorr->correction  ( refJetc   );
    refJetc  .setP4(l7.c  *refJet.p4());
  }
  if (bl7b_){
    l7.b   = L7bJetCorr->correction  ( refJetb   );
    refJetb  .setP4(l7.b  *refJet.p4());
  }

  // create the actual object with scalefactors
  pat::JetCorrFactors aJetCorr( l1, l2, l3, l4, l5, l6, l7 );
  o.setJetCorrFactors(aJetCorr);
  o.setP4(fabs(aJetCorr.correction(pat::JetCorrFactors::L3)) * o.p4());

//   std::cout << "JET1: " << o.jetCorrFactors().scaleDefault() << std::endl;

//   std::cout << "JET1a: " 
// 	    << std::endl
// 	    << l1 << " "
// 	    << l2 << " "
// 	    << l3 << " "
// 	    << l4 << " "
// 	    << std::endl
// 	    << l5.g << " "
// 	    << l5.uds << " "
// 	    << l5.c << " "
// 	    << l5.b << " "
// 	    << std::endl
// 	    << l6 << " "
// 	    << std::endl
// 	    << l7.g << " "
// 	    << l7.uds << " "
// 	    << l7.c << " "
// 	    << l7.b << " "
// 	    << std::endl;

//   std::cout << "JET1b: " 
// 	    << std::endl
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L1) << " " 
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L2) << " " 
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L3) << " " 
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L4) << " " 
// 	    << std::endl
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L5g) << " " 
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L5uds) << " " 
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L5c) << " " 
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L5b) << " " 
// 	    << std::endl
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L6g) << " " 
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L6uds) << " " 
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L6c) << " " 
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L6b) << " " 
// 	    << std::endl
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L7g) << " " 
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L7uds) << " " 
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L7c) << " " 
// 	    << o.jetCorrFactors().correction(pat::JetCorrFactors::L7b) << " " 
// 	    << std::endl;

  LogDebug("PatCrossCleaner") << "L3 jet correction factor for changed p4 is "
			      << fabs(aJetCorr.correction(pat::JetCorrFactors::L3))
			      << "\n Changed jet corrected p4 is " << o.p4()
			      << "\n for correction step " << o.jetCorrName()
			      << "\n Changed jet 'raw' p4 is     " << o.correctedJet("RAW").p4();

}

template <class Object> 
CrossCleanerMETCorrection PatCrossCleaner::metCorr( Object& o, const Object& orig ) {
  return CrossCleanerMETCorrection();
}

CrossCleanerMETCorrection PatCrossCleaner::metCorr( pat::Jet& o, const pat::Jet& orig ) {
  
  CrossCleanerMETCorrection metCorr( orig.correctedJet("RAW").px() - orig.px(), 
				     orig.correctedJet("RAW").py() - orig.py() );
  
  metCorr += CrossCleanerMETCorrection( o.px() - o.correctedJet("RAW").px(), 
					o.py() - o.correctedJet("RAW").py() );
  
  return metCorr;
  
}

#endif
