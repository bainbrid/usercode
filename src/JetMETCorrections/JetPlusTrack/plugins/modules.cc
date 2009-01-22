#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
DEFINE_SEAL_MODULE();

#include "JetMETCorrections/JetPlusTrack/interface/EnergyScaleHistogrammer.h"
DEFINE_ANOTHER_FWK_SERVICE(EnergyScaleHistogrammer);

#include "JetMETCorrections/JetPlusTrack/plugins/EnergyScaleAnalyzer.h"
DEFINE_ANOTHER_FWK_MODULE(EnergyScaleAnalyzer);

