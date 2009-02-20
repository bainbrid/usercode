#include "SusyAnalysis/PatCrossCleaner/interface/CrossCleanerResult.h"

namespace edm{
  bool operator< (const edm::RefToBase<reco::Candidate> & r1,
		  const edm::RefToBase<reco::Candidate> & r2){
    bool result=(r1.id()<r2.id());
    if (r1.id()==r2.id()) result = (r1.key()<r2.key());
    return result; 
  }
}
