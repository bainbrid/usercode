#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_SEAL_MODULE();

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "bainbrid/Test/test/JPTCorrectionProducer.cc"
typedef cms::JPTCorrectionProducer<pat::Jet> PatJetCorrectionProducer;
DEFINE_ANOTHER_FWK_MODULE(PatJetCorrectionProducer);

#include "FWCore/Framework/interface/SourceFactory.h"
#include "bainbrid/Test/test/JPTCorrectionService.h"
DEFINE_FWK_EVENTSETUP_SOURCE(JPTCorrectionService);
