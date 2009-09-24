#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_SEAL_MODULE();

using namespace cms;

#include "JetMETCorrections/Modules/interface/JetCorrectionProducer.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
typedef JetCorrectionProducer<pat::Jet> PatJetCorrectionProducer;
DEFINE_ANOTHER_FWK_MODULE(PatJetCorrectionProducer);

#include "FWCore/Framework/interface/SourceFactory.h"
#include "JetMETCorrections/Modules/interface/JetCorrectionService.h"
#include "bainbrid/Test/test/JPTCorrector.h"
DEFINE_JET_CORRECTION_SERVICE(JPTCorrector,JPTCorrectionService);

#include "FWCore/Framework/interface/SourceFactory.h"
#include "JetMETCorrections/Modules/interface/JetCorrectionService.h"
#include "bainbrid/Test/test/JPTPatCorrector.h"
DEFINE_JET_CORRECTION_SERVICE(JPTPatCorrector,JPTPatCorrectionService);
