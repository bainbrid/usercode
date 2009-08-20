#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Common/interface/View.h"
typedef edm::View<reco::LeafCandidate> Candidates;

#include "CommonTools/UtilAlgos/interface/NtpProducer.h"
typedef NtpProducer<Candidates> LorentzVectorNtuplizer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LorentzVectorNtuplizer);
