#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Common/interface/View.h"
typedef edm::View<reco::LeafCandidate> Candidates;

#include "CommonTools/UtilAlgos/interface/NtpProducer.h"
typedef NtpProducer<Candidates> CandidateNtuplizer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CandidateNtuplizer);
