//
// $Id: PATMuonProducer.cc,v 1.19.2.1 2008/11/25 15:39:40 gpetrucc Exp $
//

#include "PhysicsTools/PatAlgos/plugins/PATMuonProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/ParticleFlowCandidate/interface/IsolatedPFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/IsolatedPFCandidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/Association.h"

#include "TMath.h"

#include <vector>
#include <memory>


using namespace pat;
using namespace std;


PATMuonProducer::PATMuonProducer(const edm::ParameterSet & iConfig) :
  isolator_(iConfig.exists("isolation") ? iConfig.getParameter<edm::ParameterSet>("isolation") : edm::ParameterSet(), false),
  userDataHelper_ ( iConfig.getParameter<edm::ParameterSet>("userData") )

{

  
  // general configurables
  muonSrc_             = iConfig.getParameter<edm::InputTag>( "muonSource" );
  pfMuonSrc_           = iConfig.getParameter<edm::InputTag>( "pfMuonSource" );
  useParticleFlow_        = iConfig.getParameter<bool>( "useParticleFlow" );

  embedTrack_          = iConfig.getParameter<bool>         ( "embedTrack" );
  embedStandAloneMuon_ = iConfig.getParameter<bool>         ( "embedStandAloneMuon" );
  embedCombinedMuon_   = iConfig.getParameter<bool>         ( "embedCombinedMuon" );
  embedPFCandidate_   = iConfig.getParameter<bool>( "embedPFCandidate" );
  
  
  // MC matching configurables
  addGenMatch_   = iConfig.getParameter<bool>         ( "addGenMatch" );
  if (addGenMatch_) {
      embedGenMatch_ = iConfig.getParameter<bool>         ( "embedGenMatch" );
      if (iConfig.existsAs<edm::InputTag>("genParticleMatch")) {
          genMatchSrc_.push_back(iConfig.getParameter<edm::InputTag>( "genParticleMatch" ));
      } else {
          genMatchSrc_ = iConfig.getParameter<std::vector<edm::InputTag> >( "genParticleMatch" );
      }
  }
  
  // trigger matching configurables
  addTrigMatch_     = iConfig.getParameter<bool>            ( "addTrigMatch" );
  trigMatchSrc_     = iConfig.getParameter<std::vector<edm::InputTag> >( "trigPrimMatch" );
  
  // resolution configurables
  addResolutions_= iConfig.getParameter<bool>         ( "addResolutions" );
  
  // Efficiency configurables
  addEfficiencies_ = iConfig.getParameter<bool>("addEfficiencies");
  if (addEfficiencies_) {
     efficiencyLoader_ = pat::helper::EfficiencyLoader(iConfig.getParameter<edm::ParameterSet>("efficiencies"));
  }

  if (iConfig.exists("isoDeposits")) {
     edm::ParameterSet depconf = iConfig.getParameter<edm::ParameterSet>("isoDeposits");
     if (depconf.exists("tracker")) isoDepositLabels_.push_back(std::make_pair(TrackerIso, depconf.getParameter<edm::InputTag>("tracker")));
     if (depconf.exists("ecal"))    isoDepositLabels_.push_back(std::make_pair(ECalIso, depconf.getParameter<edm::InputTag>("ecal")));
     if (depconf.exists("hcal"))    isoDepositLabels_.push_back(std::make_pair(HCalIso, depconf.getParameter<edm::InputTag>("hcal")));
     if (depconf.exists("user")) {
        std::vector<edm::InputTag> userdeps = depconf.getParameter<std::vector<edm::InputTag> >("user");
        std::vector<edm::InputTag>::const_iterator it = userdeps.begin(), ed = userdeps.end();
        int key = UserBaseIso;
        for ( ; it != ed; ++it, ++key) {
            isoDepositLabels_.push_back(std::make_pair(IsolationKeys(key), *it));
        }
     }
  }

  // Check to see if the user wants to add user data
  useUserData_ = false;
  if ( iConfig.exists("userData") ) {
    useUserData_ = true;
  }

  // produces vector of muons
  produces<std::vector<Muon> >();

}


PATMuonProducer::~PATMuonProducer() {
}

void PATMuonProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
  
  if (isolator_.enabled()) isolator_.beginEvent(iEvent,iSetup);

  if (efficiencyLoader_.enabled()) efficiencyLoader_.newEvent(iEvent);

  std::vector<edm::Handle<edm::ValueMap<IsoDeposit> > > deposits(isoDepositLabels_.size());
  for (size_t j = 0, nd = deposits.size(); j < nd; ++j) {
    iEvent.getByLabel(isoDepositLabels_[j].second, deposits[j]);
  }

  // prepare the MC matching
  GenAssociations  genMatches(genMatchSrc_.size());
  if (addGenMatch_) {
    for (size_t j = 0, nd = genMatchSrc_.size(); j < nd; ++j) {
      iEvent.getByLabel(genMatchSrc_[j], genMatches[j]);
    }
  }

  // prepare the trigger matching
  TrigAssociations  trigMatches(trigMatchSrc_.size());
  if ( addTrigMatch_ ) {
    for ( size_t i = 0; i < trigMatchSrc_.size(); ++i ) {
      iEvent.getByLabel(trigMatchSrc_[i], trigMatches[i]);
    }
  }
  

  std::vector<Muon> * patMuons = new std::vector<Muon>();

  // loop over muons
  // Get the collection of muons from the event
  

  if( useParticleFlow_ ) {
    edm::Handle< reco::IsolatedPFCandidateCollection >  pfMuons;
    iEvent.getByLabel(pfMuonSrc_, pfMuons);
    unsigned index=0;
    for( reco::IsolatedPFCandidateConstIterator i = pfMuons->begin(); 
	 i != pfMuons->end(); ++i, ++index) {
      
      const reco::IsolatedPFCandidate& pfmu = *i;
      // std::cout<<pfmu<<std::endl;

      const reco::MuonRef& muonRef = pfmu.muonRef();
      assert( muonRef.isNonnull() );


      MuonBaseRef muonBaseRef(muonRef);
      Muon aMuon(muonBaseRef);

      reco::IsolatedPFCandidateRef pfRef( pfMuons, index );
      reco::CandidateBaseRef pfBaseRef( pfRef ); 
      
      fillMuon( aMuon, muonBaseRef, pfBaseRef, genMatches, trigMatches);
      
      aMuon.setPFCandidateRef( pfRef );
      if( embedPFCandidate_ ) aMuon.embedPFCandidate();
      
      patMuons->push_back(aMuon);
      
    } 
  }
  else {
    edm::Handle<edm::View<MuonType> > muons;
    iEvent.getByLabel(muonSrc_, muons);
    for (edm::View<MuonType>::const_iterator itMuon = muons->begin(); itMuon != muons->end(); ++itMuon) {
      
      
      // construct the Muon from the ref -> save ref to original object
      unsigned int idx = itMuon - muons->begin();
      MuonBaseRef muonRef = muons->refAt(idx);
      reco::CandidateBaseRef muonBaseRef( muonRef ); 
      
      Muon aMuon(muonRef);
      
      fillMuon( aMuon, muonRef, muonBaseRef, genMatches, trigMatches );
      
      // Isolation
      if (isolator_.enabled()) {
	isolator_.fill(*muons, idx, isolatorTmpStorage_);
	typedef pat::helper::MultiIsolator::IsolationValuePairs IsolationValuePairs;
	// better to loop backwards, so the vector is resized less times
	for (IsolationValuePairs::const_reverse_iterator it = isolatorTmpStorage_.rbegin(), ed = isolatorTmpStorage_.rend(); it != ed; ++it) {
	  aMuon.setIsolation(it->first, it->second);
	}
      }
      
      for (size_t j = 0, nd = deposits.size(); j < nd; ++j) {
	aMuon.setIsoDeposit(isoDepositLabels_[j].first, 
			    (*deposits[j])[muonRef]);
      }

      // add sel to selected
      edm::Ptr<MuonType> muonsPtr = muons->ptrAt(idx);
      if ( useUserData_ ) {
	userDataHelper_.add( aMuon, iEvent, iSetup );
      }

      patMuons->push_back(aMuon);
    }
    
  }

  // sort muons in pt
  std::sort(patMuons->begin(), patMuons->end(), pTComparator_);

  // put genEvt object in Event
  std::auto_ptr<std::vector<Muon> > ptr(patMuons);
  iEvent.put(ptr);

  if (isolator_.enabled()) isolator_.endEvent();
}

void PATMuonProducer::fillMuon( Muon& aMuon, 
				const MuonBaseRef& muonRef,
				const reco::CandidateBaseRef& baseRef,
				const GenAssociations& genMatches,
				const TrigAssociations& trigMatches ) const {
  

  if (embedTrack_) aMuon.embedTrack();
  if (embedStandAloneMuon_) aMuon.embedStandAloneMuon();
  if (embedCombinedMuon_) aMuon.embedCombinedMuon();
  
  // store the match to the generated final state muons
  if (addGenMatch_) {
    for(size_t i = 0, n = genMatches.size(); i < n; ++i) {      
      reco::GenParticleRef genMuon = (*genMatches[i])[baseRef];
      aMuon.addGenParticleRef(genMuon);
    }
    if (embedGenMatch_) aMuon.embedGenParticle();
  }

  // matches to trigger primitives
  if ( addTrigMatch_ ) {
    for ( size_t i = 0; i < trigMatches.size(); ++i ) {
      TriggerPrimitiveRef trigPrim = (*trigMatches[i])[baseRef];
      if ( trigPrim.isNonnull() && trigPrim.isAvailable() ) {
	aMuon.addTriggerMatch(*trigPrim);
      }
    }
  }
  
  if (efficiencyLoader_.enabled()) {
    efficiencyLoader_.setEfficiencies( aMuon, muonRef );
  }
  

}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATMuonProducer);
