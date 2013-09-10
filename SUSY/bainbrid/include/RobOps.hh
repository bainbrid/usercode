#ifndef hadronic_include_RobOps_hh
#define hadronic_include_RobOps_hh

#include "EventData.hh"
#include "Jet.hh"
#include "Operation.hh"
#include "Types.hh"
#include "Math/VectorUtil.h"
#include "TH1F.h"


namespace Operation {
 
  class RobAlphaT : public Operation::_Base {
  public:
    RobAlphaT( float );
    ~RobAlphaT() {;} 
    bool Process( Event::Data& );
    std::ostream& Description( std::ostream& );
  private:
    float cut_;
  }; 

  /**
   */
  class RobOps : public Operation::_Base {
    
  public:

    RobOps( const Utils::ParameterSet& );

    ~RobOps();

    void Start( Event::Data& );

    bool Process( Event::Data& );

    std::ostream& Description( std::ostream& );

    // Utility methods

    static std::pair<double,double> dalitz( const Event::Data&, 
					    LorentzV& mht,
					    LorentzV& lv1,
					    LorentzV& lv2,
					    bool thrust = false );
    
    static std::vector<LorentzV> genJets( const Event::Data& );
    
    static double ht( const std::vector<LorentzV>& );
    static LorentzV mht( const std::vector<LorentzV>& );
    static double meff( const std::vector<LorentzV>& );
    
  private:
    
    std::string algo_;
    double cut_;
    std::vector<double> jets_;

  }; 
  
} 

#endif // hadronic_include_RobOps_hh
