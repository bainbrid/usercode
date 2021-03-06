//
// $Id: PATTauProducer.cc,v 1.18.2.1 2008/11/25 15:39:40 gpetrucc Exp $
//

#include "PhysicsTools/PatAlgos/plugins/PATTauProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/TauReco/interface/PFTau.h"
//#include "DataFormats/TauReco/interface/PFTauDiscriminatorByIsolation.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TauReco/interface/CaloTau.h"
#include "DataFormats/TauReco/interface/CaloTauDiscriminatorByIsolation.h"
#include "DataFormats/TauReco/interface/CaloTauDiscriminator.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

// #include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
// #include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
// #include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include <vector>
#include <memory>


using namespace pat;


PATTauProducer::PATTauProducer(const edm::ParameterSet & iConfig):
  userDataHelper_ ( iConfig.getParameter<edm::ParameterSet>("userData") )
{
  // initialize the configurables
  tauSrc_               = iConfig.getParameter<edm::InputTag>( "tauSource" );

  embedIsolationTracks_ = iConfig.getParameter<bool>         ( "embedIsolationTracks" );
  embedLeadTrack_       = iConfig.getParameter<bool>         ( "embedLeadTrack" );
  embedSignalTracks_    = iConfig.getParameter<bool>         ( "embedSignalTracks" );

  addGenMatch_    = iConfig.getParameter<bool>         ( "addGenMatch" );

  if (addGenMatch_) {
      embedGenMatch_ = iConfig.getParameter<bool>         ( "embedGenMatch" );
      if (iConfig.existsAs<edm::InputTag>("genParticleMatch")) {
          genMatchSrc_.push_back(iConfig.getParameter<edm::InputTag>( "genParticleMatch" ));
      } else {
          genMatchSrc_ = iConfig.getParameter<std::vector<edm::InputTag> >( "genParticleMatch" );
      }
  }

  addGenJetMatch_    = iConfig.getParameter<bool>         ( "addGenJetMatch" );
  if(addGenJetMatch_) {
    embedGenJetMatch_  = iConfig.getParameter<bool>         ( "embedGenJetMatch" );
    genJetMatchSrc_    = iConfig.getParameter<edm::InputTag>( "genJetMatch" );
  }

  addTrigMatch_   = iConfig.getParameter<bool>               ( "addTrigMatch" );
  trigMatchSrc_   = iConfig.getParameter<std::vector<edm::InputTag> >( "trigPrimMatch" );
  addResolutions_ = iConfig.getParameter<bool>         ( "addResolutions" );

  // tau ID configurables
  addTauID_       = iConfig.getParameter<bool>         ( "addTauID" );
  if ( addTauID_ ) {
    // it might be a single tau ID
    if (iConfig.existsAs<edm::InputTag>("tauIDSource")) {
      tauIDSrcs_.push_back(NameTag("", iConfig.getParameter<edm::InputTag>("tauIDSource")));
    }
    // or there might be many of them
    if (iConfig.existsAs<edm::ParameterSet>("tauIDSources")) {
      // please don't configure me twice
      if (!tauIDSrcs_.empty()) throw cms::Exception("Configuration") << 
	"PATTauProducer: you can't specify both 'tauIDSource' and 'tauIDSources'\n";
      // read the different tau ID names
      edm::ParameterSet idps = iConfig.getParameter<edm::ParameterSet>("tauIDSources");
      std::vector<std::string> names = idps.getParameterNamesForType<edm::InputTag>();
      for (std::vector<std::string>::const_iterator it = names.begin(), ed = names.end(); it != ed; ++it) {
	tauIDSrcs_.push_back(NameTag(*it, idps.getParameter<edm::InputTag>(*it)));
      }
    }
    // but in any case at least once
    if (tauIDSrcs_.empty()) throw cms::Exception("Configuration") <<
      "PATTauProducer: id addTauID is true, you must specify either:\n" <<
      "\tInputTag tauIDSource = <someTag>\n" << "or\n" <<
      "\tPSet tauIDSources = { \n" <<
      "\t\tInputTag <someName> = <someTag>   // as many as you want \n " <<
      "\t}\n";
  }

  // Efficiency configurables
  addEfficiencies_ = iConfig.getParameter<bool>("addEfficiencies");
  if (addEfficiencies_) {
     efficiencyLoader_ = pat::helper::EfficiencyLoader(iConfig.getParameter<edm::ParameterSet>("efficiencies"));
  }

  useUserData_ = false;
  if ( iConfig.exists("userData") ) {
    useUserData_ = true;
  }

  // produces vector of taus
  produces<std::vector<Tau> >();
}


PATTauProducer::~PATTauProducer() {
}


void PATTauProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {     

  if (efficiencyLoader_.enabled()) efficiencyLoader_.newEvent(iEvent);

  std::auto_ptr<std::vector<Tau> > patTaus(new std::vector<Tau>()); 

  edm::Handle<edm::View<TauType> > anyTaus;
  try {
    iEvent.getByLabel(tauSrc_, anyTaus);
  } catch (const edm::Exception &e) {
    edm::LogWarning("DataSource") << "WARNING! No Tau collection found. This missing input will not block the job. Instead, an empty tau collection is being be produced.";
    iEvent.put(patTaus);
    return;
  }
   
  // prepare the MC matching
  std::vector<edm::Handle<edm::Association<reco::GenParticleCollection> > > genMatches(genMatchSrc_.size());
  if (addGenMatch_) {
        for (size_t j = 0, nd = genMatchSrc_.size(); j < nd; ++j) {
            iEvent.getByLabel(genMatchSrc_[j], genMatches[j]);
        }
  }

  edm::Handle<edm::Association<reco::GenJetCollection> > genJetMatch;  
  if (addGenJetMatch_) iEvent.getByLabel(genJetMatchSrc_, genJetMatch); 

  for (size_t idx = 0, ntaus = anyTaus->size(); idx < ntaus; ++idx) {
    edm::RefToBase<TauType> tausRef = anyTaus->refAt(idx);
    edm::Ptr<TauType> tausPtr = anyTaus->ptrAt(idx);
    
    Tau aTau(tausRef);
    if (embedLeadTrack_)       aTau.embedLeadTrack();
    if (embedSignalTracks_)    aTau.embedSignalTracks();
    if (embedIsolationTracks_) aTau.embedIsolationTracks();

    // store the match to the generated final state muons
    if (addGenMatch_) {
      for(size_t i = 0, n = genMatches.size(); i < n; ++i) {
          reco::GenParticleRef genTau = (*genMatches[i])[tausRef];
          aTau.addGenParticleRef(genTau);
      }
      if (embedGenMatch_) aTau.embedGenParticle();
    }
    
    // store the match to the visible part of the generated tau
    if (addGenJetMatch_) {
      reco::GenJetRef genJetTau = (*genJetMatch)[tausRef];
      if (genJetTau.isNonnull() && genJetTau.isAvailable() ) {
        aTau.setGenJet( genJetTau );
      } // leave empty if no match found
    }
        
    // matches to trigger primitives
    if ( addTrigMatch_ ) {
      for ( size_t i = 0; i < trigMatchSrc_.size(); ++i ) {
        edm::Handle<edm::Association<TriggerPrimitiveCollection> > trigMatch;
        iEvent.getByLabel(trigMatchSrc_[i], trigMatch);
        TriggerPrimitiveRef trigPrim = (*trigMatch)[tausRef];
        if ( trigPrim.isNonnull() && trigPrim.isAvailable() ) {
          aTau.addTriggerMatch(*trigPrim);
        }
      }
    }

    // prepare ID extraction 
    if ( addTauID_ ) {
      std::vector<pat::Tau::IdPair> ids(tauIDSrcs_.size());
      for ( size_t i = 0; i < tauIDSrcs_.size(); ++i ) {
	edm::Handle<reco::CaloTauDiscriminator> caloTauIdDiscr;
	iEvent.getByLabel(tauIDSrcs_[i].second, caloTauIdDiscr);

	edm::Handle<reco::PFTauDiscriminator> pfTauIdDiscr;
	iEvent.getByLabel(tauIDSrcs_[i].second, pfTauIdDiscr);
	
	if ( typeid(*tausRef) == typeid(reco::PFTau) ) {
	  //std::cout << "filling PFTauDiscriminator '" << tauIDSrcs_[i].first << "' into pat::Tau object..." << std::endl;
	  edm::Handle<reco::PFTauCollection> pfTauCollection; 
	  iEvent.getByLabel(tauSrc_, pfTauCollection);

	  edm::Handle<reco::PFTauDiscriminator> pfTauIdDiscr;
	  iEvent.getByLabel(tauIDSrcs_[i].second, pfTauIdDiscr);

	  ids[i].first = tauIDSrcs_[i].first;
	  ids[i].second = getTauIdDiscriminator(pfTauCollection, idx, pfTauIdDiscr);
	} else if ( typeid(*tausRef) == typeid(reco::CaloTau) ) {
	  //std::cout << "filling CaloTauDiscriminator '" << tauIDSrcs_[i].first << "' into pat::Tau object..." << std::endl;
	  edm::Handle<reco::CaloTauCollection> caloTauCollection; 
	  iEvent.getByLabel(tauSrc_, caloTauCollection);

	  edm::Handle<reco::CaloTauDiscriminator> caloTauIdDiscr;
	  iEvent.getByLabel(tauIDSrcs_[i].second, caloTauIdDiscr);

	  ids[i].first = tauIDSrcs_[i].first;
	  ids[i].second = getTauIdDiscriminator(caloTauCollection, idx, caloTauIdDiscr);
	} else {
	  throw cms::Exception("Type Mismatch") <<
	    "PATTauProducer: unsupported datatype '" << typeid(*tausRef).name() << "' for tauSource\n";
	}
      }

      aTau.setTauIDs(ids);
    }

    if (efficiencyLoader_.enabled()) {
        efficiencyLoader_.setEfficiencies( aTau, tausRef );
    }

    if ( useUserData_ ) {
      userDataHelper_.add( aTau, iEvent, iSetup );
    }

    patTaus->push_back(aTau);
  }

  // sort taus in pT
  std::sort(patTaus->begin(), patTaus->end(), pTTauComparator_);

  // put genEvt object in Event
  iEvent.put(patTaus);
}

template <typename TauCollectionType, typename TauDiscrType>
float PATTauProducer::getTauIdDiscriminator(const edm::Handle<TauCollectionType>& tauCollection, size_t tauIdx, const edm::Handle<TauDiscrType>& tauIdDiscr)
{
  edm::Ref<TauCollectionType> tauRef(tauCollection, tauIdx);
  return (*tauIdDiscr)[tauRef];
}     

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATTauProducer);


