#include "JetMETCorrections/Algorithms/interface/L3AbsoluteCorrector.h"
#include "CondFormats/JetMETObjects/interface/SimpleL3AbsoluteCorrector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

using namespace std;


L3AbsoluteCorrector::L3AbsoluteCorrector (const edm::ParameterSet& fConfig) {
 std::string file="CondFormats/JetMETObjects/data/"+fConfig.getParameter <std::string> ("tagName")+".txt";
 edm::FileInPath f1(file);
 mSimpleCorrector = new SimpleL3AbsoluteCorrector (f1.fullPath());
}

L3AbsoluteCorrector::~L3AbsoluteCorrector () {
  delete mSimpleCorrector;
} 

double L3AbsoluteCorrector::correction (const LorentzVector& fJet) const {
  return mSimpleCorrector->correctionPtEta(fJet.pt(), fJet.eta());
}

double L3AbsoluteCorrector::correction (const reco::Jet& fJet) const {
  return correction (fJet.p4());
}
