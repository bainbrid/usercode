#ifndef bainbrid_Test_MHT_h
#define bainbrid_Test_MHT_h

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include <vector>
#include <string>
#include <map>

class MHT : public edm::EDProducer {

 public:
  
  explicit MHT( const edm::ParameterSet& );
  ~MHT() {;}
  
 private:
  
  typedef reco::LeafCandidate Candidate;
  typedef Candidate::LorentzVector LorentzVector;
  typedef std::vector<Candidate> Candidates;

  void beginJob( const edm::EventSetup& );
  void produce( edm::Event&, const edm::EventSetup& );
  double ht( const std::vector<Candidate>& );
  double mht( const std::vector<Candidate>& );
  
  bool getPhotons( const edm::Event&, Candidates& );
  bool getJets( const edm::Event&, Candidates&, Candidates& );
  bool getMuons( const edm::Event&, Candidates& );
  bool getElectrons( const edm::Event&, Candidates& );
  
  double Et( const LorentzVector& );
  LorentzVector t4Mom( const LorentzVector& );

 private:

  edm::InputTag jets_;
  edm::InputTag muons_;
  edm::InputTag electrons_;
  edm::InputTag photons_;

  edm::InputTag caloMet_;
  edm::InputTag genMet_;

  double jetEt_;
  double jetEta_;
  double jetEMfrac_;

  double muonPt_;
  double muonEta_;
  double muonTrkIso_;

  double electronPt_;
  double electronEta_;
  double electronTrkIso_;

  double photonEt_;
  double photonEta_;

  double totalHt_;
  
  uint32_t minObjects_;
  uint32_t minJets_;
  uint32_t minMuons_;
  uint32_t minElectrons_;
  uint32_t minPhotons_;
  
};

inline double MHT::Et( const LorentzVector& input ) {
  double et = input.E() * input.E() - input.Pz() * input.Pz();
  et = et < 0. ? -1.*sqrt(-1.*et) : sqrt(et);
  return et;
}

inline MHT::LorentzVector MHT::t4Mom( const LorentzVector& input ) {
  return LorentzVector( input.Px(), input.Py(), 0., Et(input) );
}

#endif // bainbrid_Test_MHT_h




