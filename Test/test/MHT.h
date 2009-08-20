#ifndef bainbrid_Test_MHT_h
#define bainbrid_Test_MHT_h

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include <vector>
#include <string>
#include <map>

class TH1D;
class TH2D;

class MHT : public edm::EDAnalyzer {

 public:
  
  explicit MHT( const edm::ParameterSet& );
  ~MHT() {;}
  
 private:
  
  typedef reco::LeafCandidate Candidate;
  typedef Candidate::LorentzVector LorentzVector;

  void beginJob( const edm::EventSetup& );
  void analyze( const edm::Event&, const edm::EventSetup& );
  double ht( const std::vector<Candidate>& );
  double mht( const std::vector<Candidate>& );
  
  bool getPhotons( const edm::Event&, std::vector<Candidate>& );
  bool getJets( const edm::Event&, std::vector<Candidate>& );
  bool getMuons( const edm::Event&, std::vector<Candidate>& );
  bool getElectrons( const edm::Event&, std::vector<Candidate>& );
  
  TH1D* histo( const std::string& histogram_name );
  TH2D* histo2d( const std::string& histogram_name );

  double Et( const LorentzVector& );
  LorentzVector t4Mom( const LorentzVector& );

 private:

  edm::InputTag jets_;
  edm::InputTag muons_;
  edm::InputTag electrons_;
  edm::InputTag photons_;

  edm::InputTag met_;  
  edm::InputTag ccMet_;  
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
  
  std::map<std::string,TH1D*> histos_; 
  std::map<std::string,TH2D*> histos2d_; 
  
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




