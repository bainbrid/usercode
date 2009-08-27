#ifndef bainbrid_Test_Ntuplizer_h
#define bainbrid_Test_Ntuplizer_h

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include <vector>

class TTree;

/**
   
 */
class Ntuplizer : public edm::EDAnalyzer {

 public:
  
  explicit Ntuplizer( const edm::ParameterSet& );
  ~Ntuplizer() {;}
  
 private:

  void beginJob( const edm::EventSetup& );
  void analyze( const edm::Event&, const edm::EventSetup& );

  void init( uint32_t size );
  void reset( uint32_t size );
  void add( uint32_t index, const reco::LeafCandidate& );
  const char* name( const edm::InputTag& tag, 
		    const std::string& name, 
		    const std::string& type = "",
		    const std::string& index = "" );
  
  static const uint32_t SIZE;
  
  TTree* tree_;
  std::vector<edm::InputTag> tags_;
  std::vector<std::string> names_;

  typedef double Double;
  typedef std::vector<Double> VDouble;
  typedef std::vector<VDouble> VVDouble;
  typedef uint32_t Int;
  typedef std::vector<uint32_t> VInt;
  typedef std::vector<VInt> VVInt;
  
  VInt     n_;
  VVDouble e_;
  VVDouble et_;
  VVDouble p_;
  VVDouble pt_;
  VVDouble px_;
  VVDouble py_;
  VVDouble pz_;
  VVDouble eta_;
  VVDouble phi_;
  VVDouble mass_;
  VVInt    charge_;
  VVInt    pdgId_;
  VVDouble vx_;
  VVDouble vy_;
  VVDouble vz_;  
  
  bool k_; // kinematics
  bool t_; // transverse
  bool d_; // direction
  bool i_; // particle id
  bool v_; // vertex
  
};


#endif // bainbrid_Test_Ntuplizer_h
