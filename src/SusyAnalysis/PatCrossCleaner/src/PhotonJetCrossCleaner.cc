
#include "SusyAnalysis/PatCrossCleaner/interface/PhotonJetCrossCleaner.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace pat;
using namespace std;


void PhotonJetCrossCleaner::clean( 
             const edm::View<Photon>& Photons,
             const edm::View<Jet>& Jets,
	     CrossCleanerMap & assMap,
	     const CaloTowerConstituentsMap& constituentsMap
	     ) const
{
    for (unsigned int iPhoton=0; iPhoton!=Photons.size(); ++iPhoton)
    {
	// check photon ID first
	std::string photonID = config_.PhotonID;
	if( photonID == "LooseEM" )
	{
	    if( !Photons[iPhoton].isLooseEM() ) continue;
	}
	else if( photonID == "LoosePhoton" )
	{
	    if( !Photons[iPhoton].isLoosePhoton() ) continue;
	}
	else if( photonID == "TightPhoton" )
	{
	    if( !Photons[iPhoton].isTightPhoton() ) continue;
	}
	else
	{
	    LogDebug("PhotonJetCrossCleaner") << "No matching photon ID method. Using entire collection.";
	}

		
	edm::RefToBase<reco::Candidate> photonRef( Photons.refAt(iPhoton) );
	double dR_min = 1000;
	unsigned int iClosestJet=Jets.size();
	for (unsigned int iJet=0; iJet!=Jets.size(); ++iJet)
	{
	    double dR = ::deltaR(Photons[iPhoton], Jets[iJet]);
	    if ( dR < dR_min )
	    {
		dR_min = dR;
		iClosestJet = iJet;
	    }
	}

        if ( dR_min > config_.deltaR_min ) continue;

	if (iClosestJet==Jets.size()) continue;

	bool isolated = isIsolated_(Photons[iPhoton]);
	double sharedE = SharedEnergy_(Photons[iPhoton],Jets[iClosestJet],constituentsMap);
	edm::RefToBase<reco::Candidate> jetRef( Jets.refAt(iClosestJet) );
	if( sharedE > 0 && isolated )
	{
            if( sharedE > Jets[iClosestJet].correctedJet("RAW").energy()-0.01 )
	    {
	       	assMap[jetRef].modifiers.push_back(CrossCleanerModifier(photonRef));
		LogDebug("PhotonJetCrossCleaner") << "photon/jet overlap, dropping the jet";
	    }
            else
	    {
		assMap[jetRef].modifiers.push_back(CrossCleanerModifier(photonRef, -sharedE));
		LogDebug("PhotonJetCrossCleaner") << "photon/jet overlap, modifying the jet energy by: "<<-sharedE;
	    }
	}
	else if( sharedE > 0 && !isolated )
	{
	    double energyCorr=Photons[iPhoton].energy()-sharedE;
	    assMap[photonRef].modifiers.push_back(CrossCleanerModifier(jetRef));
	    assMap[jetRef].modifiers.push_back(CrossCleanerModifier(photonRef, energyCorr));
	    LogDebug("PhotonJetCrossCleaner") << "photon/jet overlap. dropping the photon, modifying the jet energy by: "<<energyCorr;
	}
    }
}


double PhotonJetCrossCleaner::SharedEnergy_( 
             const pat::Photon& emobject,
             const pat::Jet& jet,
	     const CaloTowerConstituentsMap& constituentsMap ) const
{
  double result = 0.0;

  std::map<CaloTowerDetId, double> etowers;
  
  // Photon supercluster and its crystals
  reco::SuperClusterRef sc = emobject.superCluster();
  std::vector<DetId>  scXtals = sc->getHitsByDetId();
  
  // Jet calo towers
  //  std::vector<CaloTowerRef> jt = jet.originalObjectRef()->getConstituents();
  std::vector<CaloTowerPtr> jt = jet.getCaloConstituents();
  
  for(unsigned int xi=0; xi<scXtals.size(); ++xi)
  {
      CaloTowerDetId towerDetId = constituentsMap.towerOf(scXtals[xi]);

      for (std::vector <CaloTowerPtr>::const_iterator tow = jt.begin(); tow!=jt.end(); ++tow)
      {
	  if( (*tow)->id()==towerDetId)
	  {
	      etowers[(*tow)->id()] =  (*tow)->emEnergy();
	  }
      }

  }

  // Or search the other way round
  // (as soon as CaloTowerConstituentsMap::constituetsOf() function works)

  for(std::map< CaloTowerDetId, double >::const_iterator itMap=etowers.begin(); itMap!=etowers.end(); ++itMap )
  {
      result+=itMap->second;
  }

  // Make sure that the shared energy is not larger than the energy of the photon
  // This can happen if there is much EM energy in a tower but the supercluster
  // does actually only contain very few crystal of this tower  
  if( result > emobject.energy() )
      return emobject.energy();
  else
    return result;
}

bool PhotonJetCrossCleaner::isIsolated_(const pat::Photon& photon) const
{
    float isoValue = photon.isolation(this->isolationMethod_);
    return ( isoValue < config_.IsoValueCut );
}

void PhotonJetCrossCleaner::setIsolationMethod(const std::string isoMethod)
{
    IsolationKeys isoKey;
    if( isoMethod == "CaloIso" )
	isoKey = CaloIso;
    else if( isoMethod == "ECalIso" )
	isoKey = ECalIso;
    else if( isoMethod == "HCalIso" )
	isoKey = HCalIso;
    else if( isoMethod == "TrackerIso" )
	isoKey = TrackerIso;
    else if( isoMethod == "User1Iso" )
	isoKey = User1Iso;
    else if( isoMethod == "User2Iso" )
	isoKey = User2Iso;
    else if( isoMethod == "User3Iso" )
	isoKey = User3Iso;
    else if( isoMethod == "User4Iso" )
	isoKey = User4Iso;
    else if( isoMethod == "User5Iso" )
	isoKey = User5Iso;
    else
    {
    	LogDebug("PhotonJetCrossCleaner") << "No match for isolation method. Using CaloIso";
	isoKey = CaloIso;
    }

    isolationMethod_ =  isoKey;
}
