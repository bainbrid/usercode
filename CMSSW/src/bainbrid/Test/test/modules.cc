#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
DEFINE_SEAL_MODULE();

#include "bainbrid/Test/test/EnergyScaleHistogrammer.h"
DEFINE_ANOTHER_FWK_SERVICE(EnergyScaleHistogrammer);

#include "bainbrid/Test/test/EnergyScaleAnalyzer.h"
DEFINE_ANOTHER_FWK_MODULE(EnergyScaleAnalyzer);

