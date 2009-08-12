#ifndef bainbrid_Test_JPTCorrectionService_h
#define bainbrid_Test_JPTCorrectionService_h

#include "bainbrid/Test/test/JPTCorrector.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "boost/shared_ptr.hpp"

class JPTCorrectionService : public edm::ESProducer, public edm::EventSetupRecordIntervalFinder {

public:
  
  JPTCorrectionService( const edm::ParameterSet& );
  ~JPTCorrectionService();
  
  boost::shared_ptr<JetCorrector> produce( const JetCorrectionsRecord& );
  
  void setIntervalFor( const edm::eventsetup::EventSetupRecordKey&, 
		       const edm::IOVSyncValue&, 
		       edm::ValidityInterval& );
  
 private:

  boost::shared_ptr<JetCorrector> mCorrector;

};

#endif // bainbrid_Test_JPTCorrectionService_h
