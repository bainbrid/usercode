#ifndef ElectronPhotonCrossCleaner_h
#define ElectronPhotonCrossCleaner_h

/**
    \class pat::ElectronPhotonCrossCleaner ElectronPhotonCrossCleaner.h "SusyAnalysis/PatUtils/ElectronPhotonCrossCleaner.h"
    \brief cross cleans objets

    \version $Id: ElectronPhotonCrossCleaner.h,v 1.1 2008/07/16 17:39:22 bmura Exp $
**/

#include "SusyAnalysis/PatCrossCleaner/interface/CrossCleanerResult.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include <DataFormats/Common/interface/Handle.h>

namespace pat {
  class ElectronPhotonCrossCleaner {
  public:
    ElectronPhotonCrossCleaner() {}
    ~ElectronPhotonCrossCleaner() {}

    void clean(
           const edm::View< pat::Electron>& Electrons,
           const edm::View< pat::Photon>& Photons,
	   CrossCleanerMap & assMap
         ) const;

  private:

  }; // class
} // namespace

#endif
