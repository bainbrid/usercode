#ifndef hadronic_include_RobPlottingOps_hh
#define hadronic_include_RobPlottingOps_hh

#include "PlottingBase.hh"
#include "Utils.hh"
#include "Types.hh"

class TH1D;
class TH2D;

namespace Operation {

  class RobPlottingOps : public PlottingBase {

  public:

    RobPlottingOps( const Utils::ParameterSet& );
    ~RobPlottingOps();

    void Start( Event::Data& );
    bool Process( Event::Data& );

    std::ostream& Description( std::ostream& );

  private:

    void BookHistos();
    LorentzV RecoilMET( Event::Data& );

    std::string dirName_;
    UInt_t nMin_;
    UInt_t nMax_;
    bool verbose_;

    // Dalitz plots
    bool dalitz_;
    std::vector<TH2D*> hDalitzRhoARhoB_;
    std::vector<TH2D*> hDalitzRhoARhoB2_;
    std::vector<TH2D*> hDalitzRhoARhoB3_;
    std::vector<TH2D*> hDalitzRhoARhoBn_;
    std::vector<TH2D*> hDalitzMTj1Mjj_;
    std::vector<TH2D*> hDalitzMTj2Mjj_;
    std::vector<TH2D*> hDalitzMTj1MTjj_;
    std::vector<TH2D*> hDalitzMTj2MTjj_;
    std::vector<TH2D*> hDalitzMT2j1M2jj_;
    std::vector<TH2D*> hDalitzMT2j2M2jj_;
    std::vector<TH2D*> hDalitzMT2j1MT2jj_;
    std::vector<TH2D*> hDalitzMT2j2MT2jj_;

    void dalitz();
    bool dalitz( Event::Data& );
    
    // AlphaT plots
    bool alphaT_;
    std::vector<TH1D*> hAlphaT_;
    std::vector<TH2D*> hAlphaTVal_;
    std::vector<TH2D*> hAlphaTnJets_;
    void alphaT();
    bool alphaT( Event::Data& );

    // Gen PtHat plots
    bool ptHat_;
    std::vector<TH1D*> hPtHat_;
    void ptHat();
    bool ptHat( Event::Data& );

    // MET plots
    bool met_;
    std::vector<TH1D*> hMetDiff1_;
    std::vector<TH1D*> hMetDiff2_;
    std::vector<TH1D*> hMetCC_;
    std::vector<TH1D*> hMetRaw_;
    std::vector<TH1D*> hMetCalo_;
    void met();
    bool met( Event::Data& );

    // Event summary
    bool summary_;
    bool summary(  Event::Data& );

    // Cross-cleaning output
    bool cc_;
    bool cc( Event::Data& );
    
    // Kinematic plots
    bool kine_;
    std::vector<TH1D*> hEta_;
    std::vector<TH1D*> hDEta_;
    std::vector<TH1D*> hDR_;
    void kine();
    bool kine( Event::Data& );
    
    std::vector<TH1D*> hHt_;
    bool ht_;
    void ht();
    bool ht( Event::Data& );

    std::vector<TH1D*> hGenHt_;
    bool genHt_;
    void genHt();
    bool genHt( Event::Data& );

    std::vector<TH1D*> hMht_;
    bool mht_;
    void mht();
    bool mht( Event::Data& );

    std::vector<TH1D*> hMeff_;
    bool meff_;
    void meff();
    bool meff( Event::Data& );

    std::vector<TH1D*> hGenMeff_;
    bool genMeff_;
    void genMeff();
    bool genMeff( Event::Data& );

    std::vector<TH1D*> hMt2_;
    bool mt2_;
    void mt2();
    bool mt2( Event::Data& );

    // Jet response plots
    bool response_;
    void response();
    bool response( Event::Data& );
    std::vector<TH2D*> hRespCorr_;
    std::vector<TH2D*> hRespEtaCorr_;
    std::vector<TH2D*> hRespRaw_;
    std::vector<TH2D*> hRespEtaRaw_;
    std::vector<TH2D*> hCorr_;
    std::vector<TH2D*> hJetID_;

    // AlphaT ratio
    bool ratio_;
    double alphaTcut_;
    void ratio();
    bool ratio( Event::Data& );
    std::vector<TH1D*> hHtEqPre_; 
    std::vector<TH1D*> hHtEqPost_; 
    std::vector<TH1D*> hGenHtEqPre_; 
    std::vector<TH1D*> hGenHtEqPost_; 
    std::vector<TH1D*> hHtGtPre_; 
    std::vector<TH1D*> hHtGtPost_; 
    std::vector<TH1D*> hGenHtGtPre_; 
    std::vector<TH1D*> hGenHtGtPost_; 
    std::vector<TH1D*> hHtLtPre_;
    std::vector<TH1D*> hHtLtPost_;
    std::vector<TH1D*> hGenHtLtPre_;
    std::vector<TH1D*> hGenHtLtPost_;
    
  }; 

}

#endif // hadronic_include_RobPlottingOps_hh
