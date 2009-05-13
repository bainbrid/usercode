#include "bainbrid/Test/test/JPTCorrectionService.h"
#include "bainbrid/Test/test/JPTCorrector.h"

// -----------------------------------------------------------------------------
//
JPTCorrectionService::JPTCorrectionService( const edm::ParameterSet& pset ) 
  : mCorrector( new JPTCorrector(pset) )
{
  std::string label = pset.getParameter<std::string>("label");
  setWhatProduced(this,label);
  findingRecord<JetCorrectionsRecord>();
}

// -----------------------------------------------------------------------------
//
JPTCorrectionService::~JPTCorrectionService() {;}

// -----------------------------------------------------------------------------
//
boost::shared_ptr<JetCorrector> JPTCorrectionService::produce( const JetCorrectionsRecord& ) {
  return mCorrector;
}
  
// -----------------------------------------------------------------------------
//
void JPTCorrectionService::setIntervalFor( const edm::eventsetup::EventSetupRecordKey&, 
					   const edm::IOVSyncValue&, 
					   edm::ValidityInterval& fIOV ) {
  fIOV = edm::ValidityInterval( edm::IOVSyncValue::beginOfTime(),
				edm::IOVSyncValue::endOfTime() );
}

// -----------------------------------------------------------------------------
//
#include "FWCore/Framework/interface/SourceFactory.h"
DEFINE_FWK_EVENTSETUP_SOURCE(JPTCorrectionService);
