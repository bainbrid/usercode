// -*- C++ -*-
//
// Package:    PatCrossCleaner
// Class:      PatCrossCleaner
// 
/**\class PatCrossCleaner PatCrossCleaner.cc SusyAnalysis/PatCrossCleaner/src/PatCrossCleaner.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Christian AUTERMANN
//         Created:  Sat Mar 22 12:58:04 CET 2008
// $Id: PatCrossCleaner.cc,v 1.5 2008/11/25 12:52:02 bmura Exp $
//
//


//system
#include <vector>
#include <memory>
//PAT
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
//DataFormats
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
//User
#include "SusyAnalysis/PatCrossCleaner/interface/PatCrossCleaner.h"

// new
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include <DataFormats/Provenance/interface/Provenance.h>
#include <DataFormats/Provenance/interface/BranchDescription.h>


using namespace pat;


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PatCrossCleaner::PatCrossCleaner(const edm::ParameterSet& iConfig) :
  _patJets       ( iConfig.getParameter<edm::InputTag>( "patJets" ) ),
  _patMets       ( iConfig.getParameter<edm::InputTag>( "patMets" ) ),
  _patMuons      ( iConfig.getParameter<edm::InputTag>( "patMuons" ) ),
  _patElectrons  ( iConfig.getParameter<edm::InputTag>( "patElectrons" ) ),
  _patPhotons    ( iConfig.getParameter<edm::InputTag>( "patPhotons" ) ),
  _patTaus       ( iConfig.getParameter<edm::InputTag>( "patTaus" ) ),
  _patCaloTowers ( iConfig.getParameter<edm::InputTag>( "patCaloTowers" ) ),
  _patTracks     ( iConfig.getParameter<edm::InputTag>( "patTracks" ) ),
  _patVertices   ( iConfig.getParameter<edm::InputTag>( "patVertices" ) ),
  L1JetCorrService_    ( iConfig.getParameter<std::string>( "L1JetCorrector" ) ),
  L2JetCorrService_    ( iConfig.getParameter<std::string>( "L2JetCorrector" ) ),
  L3JetCorrService_    ( iConfig.getParameter<std::string>( "L3JetCorrector" ) ),
  L4JetCorrService_    ( iConfig.getParameter<std::string>( "L4JetCorrector" ) ),
  L6JetCorrService_    ( iConfig.getParameter<std::string>( "L6JetCorrector" ) ),
  L5udsJetCorrService_ ( iConfig.getParameter<std::string>( "L5udsJetCorrector" ) ),
  L5gluJetCorrService_ ( iConfig.getParameter<std::string>( "L5gluonJetCorrector" ) ),
  L5cJetCorrService_   ( iConfig.getParameter<std::string>( "L5cJetCorrector" ) ),
  L5bJetCorrService_   ( iConfig.getParameter<std::string>( "L5bJetCorrector" ) ),
  L7udsJetCorrService_ ( iConfig.getParameter<std::string>( "L7udsJetCorrector" ) ),
  L7gluJetCorrService_ ( iConfig.getParameter<std::string>( "L7gluonJetCorrector" ) ),
  L7cJetCorrService_   ( iConfig.getParameter<std::string>( "L7cJetCorrector" ) ),
  L7bJetCorrService_   ( iConfig.getParameter<std::string>( "L7bJetCorrector" ) ),
  _doElectronJetCC ( iConfig.getParameter<bool>("doElectronJetCC" ) ),
  _EJselectionCfg(iConfig.getParameter<edm::ParameterSet>("ElectronJetCrossCleaning")),    
  _ElectronJetCC(reco::modules::make<ElectronJetCrossCleaner>(_EJselectionCfg)),

  _doPhotonJetCC   ( iConfig.getParameter<bool>("doPhotonJetCC" ) ),
  _PJselectionCfg(iConfig.getParameter<edm::ParameterSet>("PhotonJetCrossCleaning")),    
  _PhotonJetCC(reco::modules::make<PhotonJetCrossCleaner>(_PJselectionCfg)),

  _doMuonJetCC     ( iConfig.getParameter<bool>("doMuonJetCC" ) ),
  _MJselectionCfg(iConfig.getParameter<edm::ParameterSet>("MuonJetCrossCleaning")),
  _MuonJetCC(reco::modules::make<MuonJetCrossCleaner>(_MJselectionCfg)),
  
  _doElectronPhotonCC   ( iConfig.getParameter<bool>("doElectronPhotonCC" ) ),
  _ElectronPhotonCC()
  
{
  ccinstance_="cc";
  droppedinstance_="dropped";

  ///produces cross-cleaned collections of above objects
  //Alternative: produce cross-cleaning decision & MET correction per object
  produces<std::vector<pat::Jet> >(ccinstance_+"Jets");
  produces<std::vector<pat::MET> >(ccinstance_+"METs");
  produces<std::vector<pat::Muon> >(ccinstance_+"Muons");
  produces<std::vector<pat::Electron> >(ccinstance_+"Electrons");
  produces<std::vector<pat::Photon> >(ccinstance_+"Photons");
  produces<std::vector<pat::Tau> >(ccinstance_+"Taus");

  produces<std::vector<pat::Jet> >(droppedinstance_+"Jets");
  produces<std::vector<pat::MET> >(droppedinstance_+"METs");
  produces<std::vector<pat::Muon> >(droppedinstance_+"Muons");
  produces<std::vector<pat::Electron> >(droppedinstance_+"Electrons");
  produces<std::vector<pat::Photon> >(droppedinstance_+"Photons");
  produces<std::vector<pat::Tau> >(droppedinstance_+"Taus");

  //produces<std::vector<Track> >(); //there is nothing like this in PAT (yet?)
  //produces<std::vector<Tower> >(); //there is nothing like this in PAT (yet?)

  // booleans for jet corrections
  bl1_    = (L1JetCorrService_.compare("none")==0)    ? false : true;
  bl2_    = (L2JetCorrService_.compare("none")==0)    ? false : true;
  bl3_    = (L3JetCorrService_.compare("none")==0)    ? false : true;
  bl4_    = (L4JetCorrService_.compare("none")==0)    ? false : true;
  bl6_    = (L6JetCorrService_.compare("none")==0)    ? false : true;
  bl5uds_ = (L5udsJetCorrService_.compare("none")==0) ? false : true;
  bl5g_   = (L5gluJetCorrService_.compare("none")==0) ? false : true;
  bl5c_   = (L5cJetCorrService_.compare("none")==0)   ? false : true;
  bl5b_   = (L5bJetCorrService_.compare("none")==0)   ? false : true;
  bl7uds_ = (L7udsJetCorrService_.compare("none")==0) ? false : true;
  bl7g_   = (L7gluJetCorrService_.compare("none")==0) ? false : true;
  bl7c_   = (L7cJetCorrService_.compare("none")==0)   ? false : true;
  bl7b_   = (L7bJetCorrService_.compare("none")==0)   ? false : true;
}


PatCrossCleaner::~PatCrossCleaner()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PatCrossCleaner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    //Jets   
    Handle<edm::View<pat::Jet> > pJets;
    iEvent.getByLabel(_patJets,pJets);

    //MET
    Handle<edm::View<pat::MET> > pMets;
    iEvent.getByLabel(_patMets,pMets);
    edm::View<pat::MET> Mets = *pMets;

    //Muons   
    Handle<edm::View<pat::Muon> > pMuons;
    iEvent.getByLabel(_patMuons,pMuons);

    //Electrons   
    Handle<edm::View<pat::Electron> > pElectrons;
    iEvent.getByLabel(_patElectrons,pElectrons);

    //Photons   
    Handle<edm::View<pat::Photon> > pPhotons;
    iEvent.getByLabel(_patPhotons,pPhotons);

    //Taus   
    Handle<edm::View<pat::Tau> > pTaus;
    iEvent.getByLabel(_patTaus,pTaus);

    // CaloTowerConstituentsMap from EventSetup
    const IdealGeometryRecord& geoRecord = iSetup.get<IdealGeometryRecord>();
    ESHandle<CaloTowerConstituentsMap> constituentsMap;
    geoRecord.get(constituentsMap);

    //CaloTowers
    Handle<CaloTowerCollection> pTowers;
    iEvent.getByLabel(_patCaloTowers,pTowers);

    /*
    //Tracks
    Handle<reco::TrackCollection> pTracks;
    iEvent.getByLabel(_patTracks,pTracks);

    //Vertices
    Handle<reco::VertexCollection> pVertices;
    iEvent.getByLabel(_patVertices,pVertices);
    */

    ///The association map containing the objects that are to be dropped (key),
    ///and the objects that caused them to be dropped (value).
      CrossCleanerResult cleanerResult;
      CrossCleanerMap &  assMap = cleanerResult.map;

    ///1. Call the cross-cleaning algorithms:
    if (_doElectronJetCC)
      _ElectronJetCC.clean( *pElectrons, *pJets, *pTowers, assMap, *constituentsMap );
    if (_doPhotonJetCC)
      _PhotonJetCC.clean( *pPhotons, *pJets, assMap, *constituentsMap );
    if (_doMuonJetCC)
      _MuonJetCC.clean( *pMuons, *pJets, assMap );
    if (_doElectronPhotonCC)
      _ElectronPhotonCC.clean( *pElectrons, *pPhotons, assMap );
    //...


    ///2. Interference handling: 
    ///All object for which entries in the map "assMap" exist should be
    ///dropped. However dropping one object might solve conflicts of another 
    ///object, which therefore doesn't need to be dropped. This should be 
    ///handled, here.

      uint resolveIteration=0;
      const uint resolveIteration_max=10;
      //      bool resolved=false;
      CrossCleanerMap::iterator modifiedObjectIt = assMap.begin();
      while (resolveIteration<resolveIteration_max && modifiedObjectIt!=assMap.end()){//recalculate the end() at each test, since we are possibly changing the map on the fly
	//would this modifedObject be dropped?
	if (!modifiedObjectIt->second.keepFromModifiers()){
	  //do something, because it might be dropped because of something which needs to be dropped too.
	  std::vector<CrossCleanerModifier> ::iterator modifiersIt=modifiedObjectIt->second.modifiers.begin();
	  for (;modifiersIt!=modifiedObjectIt->second.modifiers.end();++modifiersIt){ //recalculate the end() at each test, since we are removing elements from the vector of modifiers on the fly
	    //do not bother if the modifier is not supposed to make the object dropped
	    if (modifiersIt->keepKeyObj) continue;
	    //is the modifier also to be modified?
	    CrossCleanerMap::iterator modifiedANDmodifier = assMap.find(modifiersIt->object);
	    if (modifiedANDmodifier!=assMap.end()){
	      //is the modifier also to be dropped?
	      if (!modifiedANDmodifier->second.keepFromModifiers()){

		//we have a conflict here
		//	      modifiedObjectIt->first is flagged to be dropped because of modifiersIt->object
		//	      and modifiersIt->object is flagged to be dropped too
		const edm::RefToBase<reco::Candidate> & toBeDropped = modifiedObjectIt->first;
		const edm::RefToBase<reco::Candidate> & toBeDroppedAndDropping = modifiersIt->object;
		const Provenance & prov1 = iEvent.getProvenance(toBeDropped.id());
		const Provenance & prov2 = iEvent.getProvenance(toBeDroppedAndDropping.id());
		
		LogDebug("PatCrossCleaner")<<" encounter a conflicts between:\n"
					   <<"a [1]: "<<toBeDropped->pdgId()<<" in"<< prov1.moduleLabel()<<" at index: "<<toBeDropped.key()
					   <<"and\n"
					   <<"a [2]: "<<toBeDroppedAndDropping->pdgId()<<" in"<< prov2.moduleLabel()<<" at index: "<<toBeDroppedAndDropping.key()
					   <<"\n forget about dropping [1] because of [2]";
		// need to void the dropping of toBeDropped due to toBeDroppedAndDropping. and rewind one step in the vector.
		modifiersIt=--(modifiedObjectIt->second.modifiers.erase(modifiersIt));
	      }//modifier will actually be dropped
	    }//modifier will actually be modified									
	  }//loop on modifiers
	}//modifiedObjectIt will not be kept

	//an object in the map has mo modifier anymore?
	if (modifiedObjectIt->second.modifiers.size()==0){
	  //remove the entry from the map
	  assMap.erase(modifiedObjectIt);
	  //start from the beginning again
	  //	  how to select the exact next entry, instead of rewinding it all?
	  modifiedObjectIt=assMap.begin();
	  resolveIteration++;//you cannot go back to the beginning more than N times.
	}
	else{
	  //	  otherwise, go to the next entry in the map.
	  ++modifiedObjectIt;
	}
      }//interference resolution main loop

      if (resolveIteration>=resolveIteration_max){
	edm::LogError("PatCrossCleaner")<<"could not resolve the cross cleaning interference. dropping ALL cross cleaning.";
	assMap.clear();
      }


    ///3. Produce clean and corrected collections
      //Electrons
      putObjects<pat::Electron>(cleanerResult, *pElectrons, iEvent, iSetup, "Electrons");

      //Jets      
      putObjects<pat::Jet>(cleanerResult, *pJets, iEvent, iSetup, "Jets");
      
      //Muons
      putObjects<pat::Muon>(cleanerResult, *pMuons, iEvent, iSetup, "Muons");
      
      //Photons
      putObjects<pat::Photon>(cleanerResult, *pPhotons, iEvent, iSetup, "Photons");
      
      //Taus
      putObjects<pat::Tau>(cleanerResult, *pTaus, iEvent, iSetup, "Taus");

      
      // 4. recalculate MET
      std::auto_ptr<std::vector<pat::MET> > corrMETs ( new std::vector<pat::MET> );
      corrMETs->push_back(cleanerResult.metCorrection.correct(Mets.front()));
      iEvent.put(corrMETs, ccinstance_+"METs");

      ///debugging output:
      //printConflicts(assMap);
      //printDropped(dropObj);
      //printCorrected(corrObj);

}

// --- debugging functions
void PatCrossCleaner::printConflicts(CrossCleanerMap& conflicts)
{
  using namespace std;
  for (CrossCleanerMap::const_iterator it=conflicts.begin(); 
       it!=conflicts.end(); ++it) {
     cout << "Object PDG-ID("<<setw(3)<<it->first->pdgId()
          << ") with phi="<<setw(8)<<it->first->phi()
          << ", eta="<<setw(9)<< it->first->eta()
          << ", pt="<<setw(9)<<it->first->pt()
	  << "; conflicting with the following " << it->second.modifiers.size()
	  << " object(s):" << endl;
     //if more detailed information on this object is needed, than it has to be casted back to an PATObject.	  
     for (std::vector<CrossCleanerModifier>::const_iterator
          cs=it->second.modifiers.begin(); cs!=it->second.modifiers.end(); ++cs) {
	cout << "     ->  object PDG-ID("<< setw(3)<<cs->object->pdgId() 
	     << ") with phi="<<setw(8)<<cs->object->phi()
             << ", eta="<<setw(9)<<cs->object->eta()
             << ", pT="<<setw(9)<<cs->object->pt()
             << ", energy="<<setw(9)<<cs->object->energy()
	     << endl;
     }	
  }     
}

void PatCrossCleaner::printCorrected(std::map<edm::RefToBase<reco::Candidate>, double>& corrected) const
{
  using namespace std;
  for(map<edm::RefToBase<reco::Candidate>, double>::const_iterator it = corrected.begin(); it!=corrected.end(); ++it)
  {
     cout << "Correcting object PDG-ID("<<setw(3)<<it->first->pdgId()
          << ") with phi="<<setw(8)<<it->first->phi()
          << ", eta="<<setw(9)<< it->first->eta()
          << ", pT="<<setw(9)<<it->first->pt()
          << ", energy="<<setw(9)<<it->first->energy()
	  << "; with an additional energy of " << it->second 
	  << endl;
  }
}
    
void PatCrossCleaner::printDropped(std::set<edm::RefToBase<reco::Candidate> > dropped) const
{
  using namespace std;
  for(set<edm::RefToBase<reco::Candidate> >::const_iterator it=dropped.begin(); it!=dropped.end(); ++it)
  {
     cout << "Dropping object PDG-ID("<<setw(3)<<(*it)->pdgId()
          << ") with phi="<<setw(8)<<(*it)->phi()
          << ", eta="<<setw(9)<< (*it)->eta()
          << ", pT="<<setw(9)<<(*it)->pt()
          << ", energy="<<setw(9)<<(*it)->energy()
	  << endl;
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
PatCrossCleaner::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PatCrossCleaner::endJob() {
  std::stringstream buffer;
  buffer<<" PatCrossCleaner statistics\n============================\n";
  for (std::map<int,std::pair<int,int> >::const_iterator it=_statistics.begin();
       it!=_statistics.end();++it){
    buffer << "PDG-ID "<<std::setw(3)<<it->first
          << ":  total = "<<std::setw(7)<<it->second.first+it->second.second
	  << ",  dropped = "<<std::setw(4)<<it->second.second;
     if (it->second.first>0)	  
       buffer << "  ("<<std::setw(1)<<(float)it->second.second/(it->second.first+it->second.second)*100.
               <<"%)";
     buffer<<std::endl;
  }     
  std::cout<<buffer.str();
  edm::LogInfo("PatCrossCleaner")<<buffer.str();
}

//define this as a plug-in
DEFINE_FWK_MODULE(PatCrossCleaner);
