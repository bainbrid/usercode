#ifndef ElectronJetCrossCleaner_h
#define ElectronJetCrossCleaner_h

/**
    \class pat::ElectronJetCrossCleaner ElectronJetCrossCleaner.h "SusyAnalysis/PatUtils/ElectronJetCrossCleaner.h"
    \brief cross cleans objets

    \version $Id: ElectronJetCrossCleaner.h,v 1.4 2009/02/10 13:26:14 georgia Exp $
**/

#include "SusyAnalysis/PatCrossCleaner/interface/CrossCleanerResult.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include <DataFormats/Common/interface/Handle.h>
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"

#include <cmath> // georgia

namespace pat {

  /// Structure defining the electron - jet cross cleaning
  struct ElectronJetCleaning {
    bool   SusyAnalyzerCleaning; //to choose SusyAnalyzer-like cleaning
    double deltaR_min; ///min. distance in dR between electron & jet
    double SharedEtoJetE; // for SusyAnalyzer style cleaning
    double IsoValueCut; // Cut on isolation (caloIso a.t.m.)
    double SharedEForNIsoEle; // Shared Energey treshold for not-isolated electrons
    std::string IsolationKey; // key for isolation method as defined in DataFormats/PatCandidates/interface/Isolation.h
    std::string ElectronID; // Electron ID to check before x-cleaning (only cut based so far)
  };


  class ElectronJetCrossCleaner {
  public:
    ElectronJetCrossCleaner( const ElectronJetCleaning& cfg ) : config_( cfg ) {
	this->setIsolationMethod( cfg.IsolationKey );
    }
    ~ElectronJetCrossCleaner() {}

    void clean(
           const edm::View< pat::Electron>& Electrons,
           const edm::View< pat::Jet>& Jets,
	   const CaloTowerCollection& towers,
	   CrossCleanerMap & assMap,
	   const CaloTowerConstituentsMap& constituentsMap
         ) const;

  private:
    void runSusyAnalyzerCleaning(
           const edm::View< pat::Electron>& Electrons,
           const edm::View< pat::Jet>& Jets,
	   const CaloTowerCollection& towers,
	   CrossCleanerMap & assMap,
	   const CaloTowerConstituentsMap& constituentsMap
         ) const;
    void runCleaning( // placeholder for other cleaning algorithm
           const edm::View< pat::Electron>& Electrons,
           const edm::View< pat::Jet>& Jets,
	   const CaloTowerCollection& towers,
	   CrossCleanerMap & assMap,
	   const CaloTowerConstituentsMap& constituentsMap
         ) const;
    ///calculates shared energy between electron & jet
    void SharedEnergy_( const pat::Electron& electron,
                          const pat::Jet& jet,
			  const CaloTowerConstituentsMap& constituentsMap, 
			math::XYZVector* sharedMomentum ) const; //georgia
    bool isIsolated_(const pat::Electron&) const;
    void setIsolationMethod(const std::string isoMethod);

    IsolationKeys isolationMethod_;
    ElectronJetCleaning config_;

  }; // class
} // namespace

#endif
