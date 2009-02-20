#ifndef MuonJetCrossCleaner_h
#define MuonJetCrossCleaner_h

/**
    \class pat::MuonJetCrossCleaner MuonJetCrossCleaner.h "SusyAnalysis/PatUtils/MuonJetCrossCleaner.h"
    \brief cross cleans objets

    \version $Id: MuonJetCrossCleaner.h,v 1.3 2009/02/04 17:00:46 bmura Exp $
**/

#include "SusyAnalysis/PatCrossCleaner/interface/CrossCleanerResult.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <DataFormats/Common/interface/Handle.h>

namespace pat {

  /// Structure defining the muon - jet cross cleaning
  struct MuonJetCleaning {
    double deltaR_min; ///min. distance in dR between muon & jet
    double caloIso_max; ///maximum calorimeter isolation energy in order to consider muon for removal
    double trackIso_max; ///maximum tracker isolation energy in order to consider muon for removal
    std::string muonID; ///any considered muon must fulfill ID requirements
    bool modifyJetEnergy;
  };


  class MuonJetCrossCleaner {
  public:
    MuonJetCrossCleaner( const MuonJetCleaning& cfg ) : config_( cfg ) {
	this->setMuonID(config_.muonID);
    }
    ~MuonJetCrossCleaner() {}

    void clean(
           const edm::View< pat::Muon>& Muons,
           const edm::View< pat::Jet>& Jets,
	   CrossCleanerMap & assMap
         ) const;

  private:
    void setMuonID(std::string muonID);
    MuonJetCleaning config_;
    reco::Muon::SelectionType muonID_;

  }; // class
} // namespace

#endif
