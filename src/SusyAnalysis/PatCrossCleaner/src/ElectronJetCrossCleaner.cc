
#include "SusyAnalysis/PatCrossCleaner/interface/ElectronJetCrossCleaner.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace pat;
using namespace std;


void ElectronJetCrossCleaner::clean( 
             const edm::View<Electron>& Electrons,
             const edm::View<Jet>& Jets,
	     const CaloTowerCollection& towers,
	     CrossCleanerMap & assMap,
	     const CaloTowerConstituentsMap& constituentsMap
	     ) const
{
    if (Jets.size() > 0 && Electrons.size() > 0 )
    {
	//if(config_.SusyAnalyzerCleaning)
	    this->runSusyAnalyzerCleaning(Electrons, Jets, towers, assMap, constituentsMap);
        // placeholder for other cleaning algorithm
	//else
	//    this->runCleaning(Electrons, Jets, towers, assMap, constituentsMap);
    }
}

void ElectronJetCrossCleaner::runSusyAnalyzerCleaning( 
             const edm::View<Electron>& Electrons,
             const edm::View<Jet>& Jets,
	     const CaloTowerCollection& towers,
	     CrossCleanerMap & assMap,
	     const CaloTowerConstituentsMap& constituentsMap
	     ) const
{
    for (unsigned int iElectron=0; iElectron!=Electrons.size(); ++iElectron)
    {
	// Check ElectronID first
	if( !Electrons[iElectron].electronID(config_.ElectronID) ) continue;
	 
	edm::RefToBase<reco::Candidate> electronRef( Electrons.refAt(iElectron) );
	bool isolated = isIsolated_(Electrons[iElectron]);

	double dR_min = 1000.;
	unsigned int iClosestJet=Jets.size();
	for (unsigned int iJet=0; iJet!=Jets.size(); ++iJet)
	{
	    double dR = ::deltaR(Electrons[iElectron], Jets[iJet] );
	    if( dR < dR_min )
	    {
		dR_min = dR;
		iClosestJet = iJet;
	    } 
	}
        if ( dR_min > config_.deltaR_min ) continue;

	if (iClosestJet==Jets.size()) continue;

	double sharedE = SharedEnergy_(Electrons[iElectron],Jets[iClosestJet],constituentsMap);
	edm::RefToBase<reco::Candidate> jetRef( Jets.refAt(iClosestJet) );
	if( sharedE > 0. && isolated )
	{
	    if( sharedE/Jets[iClosestJet].correctedJet("RAW").energy() > config_.SharedEtoJetE )
	    {
		assMap[jetRef].modifiers.push_back( CrossCleanerModifier(electronRef));
		LogDebug("ElectronJetCrossCleaner") << "electron/jet overlap. dropping the jet";
	    }
            else
            {
	        assMap[jetRef].modifiers.push_back(CrossCleanerModifier(electronRef, -sharedE));
	        LogDebug("ElectronJetCrossCleaner") << "electron/jet overlap.  modifying the jet energy by: "<<-sharedE;
            }
	}
	else if( sharedE > config_.SharedEForNIsoEle && !isolated )
	{
	    double energyCorr=Electrons[iElectron].energy()-sharedE;
	    assMap[electronRef].modifiers.push_back(CrossCleanerModifier(jetRef));
	    assMap[jetRef].modifiers.push_back(CrossCleanerModifier(electronRef, energyCorr));
	    LogDebug("ElectronJetCrossCleaner") << "electron/jet overlap. dropping the electron, modifying the jet energy by: "<<energyCorr;
	}
	
        //if (sharedE>0)        
        //cout << "Electron_"<<iElectron<<" ("<< Electrons.size()<<")"
        //     << ", Jet_"<<iClosestJet<<";  E_e="<<Electrons[iElectron].energy()
        //     << ", E_jet="<<Jets[iClosestJet].correctedJet("RAW").energy()
        //     << ", shared_E="<< sharedE
        //     <<endl;
    }
}



double ElectronJetCrossCleaner::SharedEnergy_( 
             const pat::Electron& emobject,
             const pat::Jet& jet,
	     const CaloTowerConstituentsMap& constituentsMap ) const
{
  double result = 0.0;

  std::map<CaloTowerDetId, double> etowers;
  
  // Electron supercluster and its crystals
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
  // (once CaloTowerConstituentsMap::constituetsOf() function works)

  for(std::map< CaloTowerDetId, double >::const_iterator itMap=etowers.begin(); itMap!=etowers.end(); ++itMap )
  {
      result+=itMap->second;
  }

  // Make sure that the shared energy is not larger than the energy of the electron
  // This can happen if there is much EM energy in a tower but the supercluster
  // does actually only contain very few crystal of this tower  
  if( result > emobject.energy() )
      return emobject.energy();
  else
  return result;
}



bool ElectronJetCrossCleaner::isIsolated_(const pat::Electron& electron) const
{
    float isoValue = electron.isolation(this->isolationMethod_);
    return ( isoValue < config_.IsoValueCut );
}

void ElectronJetCrossCleaner::setIsolationMethod(const std::string isoMethod)
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
    	LogDebug("ElectronJetCrossCleaner") << "No match for isolation method. Using CaloIso";
	isoKey = CaloIso;
    }

    isolationMethod_ =  isoKey;
}
