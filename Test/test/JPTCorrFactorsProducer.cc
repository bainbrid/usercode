#include "bainbrid/Test/test/JPTCorrFactorsProducer.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <vector>

// -----------------------------------------------------------------------------
//
JPTCorrFactorsProducer::JPTCorrFactorsProducer( const edm::ParameterSet& pset ) 
  : uncorrected_( pset.getParameter<edm::InputTag>("UncorrectedJets") ),
    corrected_( pset.getParameter<edm::InputTag>("CorrectedJets") )
{
  produces< edm::ValueMap<pat::JetCorrFactors> >();
}

// -----------------------------------------------------------------------------
//
JPTCorrFactorsProducer::~JPTCorrFactorsProducer() {;}

// -----------------------------------------------------------------------------
//
bool JPTCorrFactorsProducer::matchJetsByCaloTowers( const pat::Jet& jet1,
						    const pat::Jet& jet2 ) {
  
  std::vector< edm::Ptr<CaloTower> > towers1 = jet1.getCaloConstituents();
  std::vector< edm::Ptr<CaloTower> > towers2 = jet2.getCaloConstituents();
  
  if ( towers1.empty() || 
       towers2.empty() || 
       towers1.size() != towers2.size() ) { return false; }
  
  std::vector< edm::Ptr<CaloTower> >::const_iterator ii = towers1.begin();
  std::vector< edm::Ptr<CaloTower> >::const_iterator jj = towers1.end();
  for ( ; ii != jj; ++ii ) {
    std::vector< edm::Ptr<CaloTower> >::const_iterator iii = towers2.begin();
    std::vector< edm::Ptr<CaloTower> >::const_iterator jjj = towers2.end();
    for ( ; iii != jjj; ++iii ) { if ( *iii == *ii ) { break; } }
    if ( iii == towers2.end() ) { return false; }
  }
  
  return true;
  
}

// -----------------------------------------------------------------------------
//
void JPTCorrFactorsProducer::produce( edm::Event& event, const edm::EventSetup& setup ) {

  std::vector<pat::JetCorrFactors> corrs;

  float l1 = -1;
  float l2 = -1;
  float l3 = -1;
  float l4 = -1;
  pat::JetCorrFactors::FlavourCorrections l5;
  pat::JetCorrFactors::FlavourCorrections l6;
  pat::JetCorrFactors::FlavourCorrections l7;
  
  // Retrieve uncorrected pat::Jets
  edm::Handle< edm::View<pat::Jet> > uncorrected;
  event.getByLabel( uncorrected_, uncorrected );

  // Retrieve corrected pat::Jets
  edm::Handle< edm::View<pat::Jet> > corrected;
  event.getByLabel( corrected_, corrected );

  std::vector<bool> matched;
  matched.clear();
  matched.resize( corrected->size(), false );
  
  // Match jets  
  edm::View<pat::Jet>::const_iterator ii = uncorrected->begin();
  edm::View<pat::Jet>::const_iterator jj = uncorrected->end();
  for ( ; ii != jj; ++ii ) {
    edm::View<pat::Jet>::const_iterator iii = corrected->begin();
    edm::View<pat::Jet>::const_iterator jjj = corrected->end();
    for ( ; iii != jjj; ++iii ) {
      uint32_t index = static_cast<uint32_t>( iii - corrected->begin() );
      if ( !matched[index] ) {
	if ( matchJetsByCaloTowers( *ii, *iii ) ) {
	  if ( ii->pt() > 0. ) { l3 = iii->pt() / ii->pt(); }
	  std::string name = "JPT";
	  corrs.push_back( pat::JetCorrFactors( name,
						l1,
						l2,
						l3, //@@
						l4,
						l5,
						l6,
						l7 ) );
	  matched[index] = true;
	  continue;
	}
      }
    }
  }
  
  // Construct container
  std::auto_ptr< edm::ValueMap<pat::JetCorrFactors> > corrections( new edm::ValueMap<pat::JetCorrFactors>() );
  edm::ValueMap<pat::JetCorrFactors>::Filler filler( *corrections );

  // Populate container 
  filler.insert( uncorrected, corrs.begin(), corrs.end() );
  filler.fill(); 
  
  // Put container in Event
  event.put(corrections);

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JPTCorrFactorsProducer);
