
#include "SusyAnalysis/PatCrossCleaner/interface/ElectronPhotonCrossCleaner.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Flags.h"

using namespace pat;
using namespace std;


void ElectronPhotonCrossCleaner::clean( 
             const edm::View<pat::Electron>& Electrons,
             const edm::View<pat::Photon>& Photons,
	     CrossCleanerMap & assMap
	     ) const
{
    for (unsigned int iPhoton=0; iPhoton!=Photons.size(); ++iPhoton)
    {
	const pat::Photon& photon = Photons[iPhoton];
	reco::SuperClusterRef scPhoton = photon.superCluster();

	for (unsigned int iElectron=0; iElectron!=Electrons.size(); ++iElectron)
	{
	    const pat::Electron electron = Electrons[iElectron];
	    reco::SuperClusterRef scElectron = electron.superCluster();

	    // Check overlap: Could in principle use the pat flag
	    // if ( !pat::Flags::test(photon, pat::Flags::Overlap::Electrons) )
	    // but then I don't know the modifier.
	    // Maybe there is a better way to check if the superclusters are the same?
	    if ( scElectron==scPhoton )	
	    {
		edm::RefToBase<reco::Candidate> photonRef( Photons.refAt(iPhoton) );
		edm::RefToBase<reco::Candidate> electronRef( Electrons.refAt(iElectron) );
		//remove photon
		assMap[photonRef].modifiers.push_back( CrossCleanerModifier(electronRef));
		LogDebug("ElectronPhotonCrossCleaner") << "electron/photon overlap. dropping the photon";
	    }
	}
    }
}    

