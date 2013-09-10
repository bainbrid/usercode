#ifndef hadronic_include_DeadRegionsOps_hh
#define hadronic_include_DeadRegionsOps_hh

#include "Jet.hh"
#include "PlottingBase.hh"
#include "Types.hh"
#include "Utils.hh"
#include <vector>

class TH1D;
class TH2D;

namespace Operation {

  class DeadRegionsOps : public PlottingBase {

  public:

    DeadRegionsOps( const Utils::ParameterSet& );
    ~DeadRegionsOps();
    void Start( Event::Data& );
    bool Process( Event::Data& );
    std::ostream& Description( std::ostream& );
    
  private:

    void BookHistos();
    
    void badlyMeasuredGenJets( Event::Data& ev );

    typedef std::vector<Event::Jet const*> Jets;

    Jets::const_iterator minBiasedDeltaPhi( Jets& jets,
					    double& min_biased_dphi, 
					    double& response,
					    double& resp );
    
    LorentzV matchedGenJet( Event::Data& ev,
			    Jets::const_iterator jet );
    
    double closestDeadRegion( Jets::const_iterator jet,
			      double& deta,
			      double& dphi );
    
  private:

    typedef std::pair<double,double> Region;
    typedef std::vector<Region> Regions;
    Regions dead_;

    double minDPhi_;
    double respMax_;
    bool veto_;
    bool problem_;
    
    bool verbose_;
    bool histos_;
    std::string dirName_;
    UInt_t nMin_;
    UInt_t nMax_;
    
    std::vector<TH2D*> hDPhiNewVsOld_;

    std::vector<TH1D*> hRespGen_;
    std::vector<TH2D*> hRespProjVsGen_;
    std::vector<TH2D*> hRespRecoilVsGen_;

    std::vector<TH1D*> hRespProj_;
    std::vector<TH1D*> hRespRecoil_;
    std::vector<TH2D*> hRespRecoilVsProj_;

    std::vector<TH2D*> hRespProjVsDPhi_;
    std::vector<TH2D*> hRespRecoilVsDPhi_;

    std::vector<TH1D*> hDEta_;
    std::vector<TH1D*> hDPhi_;
    std::vector<TH1D*> hDR_;

    std::vector<TH2D*> hGenEtaVsGenPhi_;
    std::vector<TH2D*> hGenEtaVsGenPhiAdj_;

  }; 

}

#endif // hadronic_include_DeadRegionsOps_hh
