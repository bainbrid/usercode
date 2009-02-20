#ifndef CrossCleanerResult_h
#define CrossCleanerResult_h

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/PatCandidates/interface/MET.h"

namespace edm{
  bool operator< (const edm::RefToBase<reco::Candidate> & r1,
		  const edm::RefToBase<reco::Candidate> & r2);
}

struct CrossCleanerModifier{
  //  CrossCleanerModifier(const bool &k, const edm::RefToBase<reco::Candidate>&o,const double& m=0)
  //    : keepKeyObj(k), object(o), modEnergy(m) {}
  CrossCleanerModifier(const edm::RefToBase<reco::Candidate>&o,const double& m)
    : keepKeyObj(true), object(o), modEnergy(m) {}
  CrossCleanerModifier(const edm::RefToBase<reco::Candidate>&o)
    : keepKeyObj(false), object(o), modEnergy(0) {}

    bool keepKeyObj; // decision to drop or modify, based on this object
    edm::RefToBase<reco::Candidate> object;
    double modEnergy;
};

class CrossCleanerValue { 
 public:
  CrossCleanerValue() : keep(true){}
  bool keep; // overall decision to drop or to modify
  std::vector< CrossCleanerModifier > modifiers;
  bool keepFromModifiers(){
    if (modifiers.size()==0) return true;
    for (uint m=0;m!=modifiers.size();++m)
      if (modifiers[m].keepKeyObj) return true;
    return false;
  }
};

typedef std::map<edm::RefToBase<reco::Candidate>, CrossCleanerValue>       CrossCleanerMap; 
     
class CrossCleanerMETCorrection {
 public:
  CrossCleanerMETCorrection(double px, double py) : dMETx(px),dMETy(py){}
  CrossCleanerMETCorrection() : dMETx(0),dMETy(0){}
  double dMETx;
  double dMETy;
  CrossCleanerMETCorrection & operator+=(const CrossCleanerMETCorrection & add){
    dMETx+=add.dMETx;
    dMETy+=add.dMETy;
    return *this;
  }
  pat::MET correct(const pat::MET & original){
    pat::MET newMET(original);
    double corrMETx = original.px() - dMETx;
    double corrMETy = original.py() - dMETy;
    double corrMETpt = sqrt(corrMETx*corrMETx + corrMETy*corrMETy);
    const pat::MET::LorentzVector corrMetp4( corrMETx, corrMETy, 0., corrMETpt );
    newMET.setP4( corrMetp4 );
    return newMET;
  }
};

class CrossCleanerResult {
 public:
  typedef CrossCleanerMap MapType;
  CrossCleanerResult(){}
  //implicit operator
    operator MapType&() { return map;}

  //map 
  std::map<edm::RefToBase<reco::Candidate>, CrossCleanerValue> map;
  //corrections to MET
  CrossCleanerMETCorrection metCorrection;
};




#endif
