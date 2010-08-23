#ifndef hadronic_include_ToyMC_hh
#define hadronic_include_ToyMC_hh

#include "PlottingBase.hh"
#include "Utils.hh"

class TFile;

namespace Operation {

  class ToyMC : public PlottingBase {

  public:

    ToyMC( const Utils::ParameterSet& );
    ~ToyMC();
    
    bool constrain( Double_t x1, 
		    Double_t x2, 
		    Double_t x3 );
    
    void integrate( Int_t depth, 
		    Int_t ndepth, 
		    Int_t* nbins, 
		    Double_t x3,
		    Double_t alpha_t,
		    Double_t& numerator,
		    Double_t& denominator, 
		    Double_t xlow, 
		    Double_t xhigh, 
		    Double_t ylow, 
		    Double_t yhigh );

  private:

    bool Process( Event::Data& ) { return true; }
    void BookHistos() {;}
    std::ostream& Description( std::ostream& os ) { return os; }

    TFile* output_;
    
  }; 

}

#endif // hadronic_include_ToyMC_hh
