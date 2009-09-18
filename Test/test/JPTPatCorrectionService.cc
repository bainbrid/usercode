#include "bainbrid/Test/test/JPTPatCorrectionService.h"
#include "bainbrid/Test/test/JPTPatCorrector.h"

// -----------------------------------------------------------------------------
//
JPTPatCorrectionService::JPTPatCorrectionService( const edm::ParameterSet& pset ) 
  : mCorrector( new JPTPatCorrector(pset) )
{
  std::string label = pset.getParameter<std::string>("label");
  setWhatProduced(this,label);
  findingRecord<JetCorrectionsRecord>();
}

// -----------------------------------------------------------------------------
//
JPTPatCorrectionService::~JPTPatCorrectionService() {;}

// -----------------------------------------------------------------------------
//
boost::shared_ptr<JetCorrector> JPTPatCorrectionService::produce( const JetCorrectionsRecord& ) {
  return mCorrector;
}
  
// -----------------------------------------------------------------------------
//
void JPTPatCorrectionService::setIntervalFor( const edm::eventsetup::EventSetupRecordKey&, 
					      const edm::IOVSyncValue&, 
					      edm::ValidityInterval& fIOV ) {
  fIOV = edm::ValidityInterval( edm::IOVSyncValue::beginOfTime(),
				edm::IOVSyncValue::endOfTime() );
}

// -----------------------------------------------------------------------------
//
#include "FWCore/Framework/interface/SourceFactory.h"
DEFINE_FWK_EVENTSETUP_SOURCE(JPTPatCorrectionService);
