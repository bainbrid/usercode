#ifndef PhotonJetCrossCleaner_h
#define PhotonJetCrossCleaner_h

/**
    \class pat::PhotonJetCrossCleaner PhotonJetCrossCleaner.h "SusyAnalysis/PatUtils/PhotonJetCrossCleaner.h"
    \brief cross cleans objects
**/

#include "SusyAnalysis/PatCrossCleaner/interface/CrossCleanerResult.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include <DataFormats/Common/interface/Handle.h>
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"

namespace pat {

  /// Structure defining the electron - jet cross cleaning
  struct PhotonJetCleaning {
    double deltaR_min; ///min. distance in dR between photon & jet
    double IsoValueCut;
    std::string IsolationKey; // key for isolation method as defined in DataFormats/PatCandidates/interface/Isolation.h
    std::string PhotonID;
  };


  class PhotonJetCrossCleaner {
  public:
    PhotonJetCrossCleaner( const PhotonJetCleaning& cfg ) : config_( cfg ) {
	this->setIsolationMethod( cfg.IsolationKey );
    }
    ~PhotonJetCrossCleaner() {}

    void clean(
           const edm::View< pat::Photon>& Photons,
           const edm::View< pat::Jet>& Jets,
	   CrossCleanerMap & assMap,
	   const CaloTowerConstituentsMap& constituentsMap
         ) const;

  private:
    ///calculates shared energy between electron & jet: based on SusyAnalyzer
    double SharedEnergy_( const pat::Photon& photon,
                          const pat::Jet& jet,
			  const CaloTowerConstituentsMap& constituentsMap ) const;
    bool isIsolated_(const pat::Photon& photon) const;
    void setIsolationMethod(const std::string isoMethod);

    IsolationKeys isolationMethod_;
    PhotonJetCleaning config_;

  }; // class
} // namespace

#endif
