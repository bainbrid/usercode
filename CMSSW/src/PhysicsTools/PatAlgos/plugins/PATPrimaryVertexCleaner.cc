#include "PhysicsTools/PatAlgos/interface/PATPrimaryVertexSelector.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"

typedef ObjectSelector<PATPrimaryVertexSelector> PATPrimaryVertexCleaner;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATPrimaryVertexCleaner);
