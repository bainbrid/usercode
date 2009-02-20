
#include "SusyAnalysis/PatCrossCleaner/interface/MuonJetCrossCleaner.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/MuonReco/interface/Muon.h"

using namespace pat;
using namespace std;


void MuonJetCrossCleaner::clean( 
             const edm::View<Muon>& Muons,
             const edm::View<Jet>& Jets,
	     CrossCleanerMap & assMap
	     ) const
{
    for (unsigned int iMuon=0; iMuon!=Muons.size(); ++iMuon)
    {
      // first check muon ID
      if (!Muons[iMuon].isGood(muonID_)) continue;
	
      // do not touch muon that are isolated
      if (Muons[iMuon].caloIso() < config_.caloIso_max) continue;
      if (Muons[iMuon].trackIso() < config_.trackIso_max) continue;

	edm::RefToBase<reco::Candidate> muonRef( Muons.refAt(iMuon) );
	for (unsigned int iJet=0; iJet!=Jets.size(); ++iJet)
       	{
	  double dR = ::deltaR(Muons[iMuon], Jets[iJet] );
	    if ( dR > config_.deltaR_min ) continue;
	    
	    edm::RefToBase<reco::Candidate> jetRef( Jets.refAt(iJet) );	    
	    double sharedEnergy = muonRef->energy();
	    LogDebug("MuonJetCrossCleaner")<<" a muon within dR :"<< dR <<" of a Jet has energy: "<<sharedEnergy
					   <<"\n modifying the jet energy AND dropping the muon.";
	    
	    assMap[muonRef].modifiers.push_back( CrossCleanerModifier(jetRef));
	    if( config_.modifyJetEnergy )
		assMap[jetRef].modifiers.push_back( CrossCleanerModifier(muonRef, sharedEnergy));
	}
    }
}    

void MuonJetCrossCleaner::setMuonID(std::string muonID)
{
    // All options from DataFormats/MuonReco/interface/Muon.h
    reco::Muon::SelectionType idType;
    if( muonID == "All" )
	idType = reco::Muon::All;                              // dummy options - always true
    else if( muonID == "AllGlobalMuons" )	
	idType = reco::Muon::AllGlobalMuons;                   // checks isGlobalMuon flag
    else if( muonID == "AllStandAloneMuons" )	
	idType = reco::Muon::AllStandAloneMuons;               // checks isStandAloneMuon flag
    else if( muonID == "AllTrackerMuons" )	
	idType = reco::Muon::AllTrackerMuons;                  // checks isTrackerMuon flag
    else if( muonID == "TrackerMuonArbitrated" )	
	idType = reco::Muon::TrackerMuonArbitrated;            // resolve ambiguity of sharing segments
    else if( muonID == "AllArbitrated" )	
	idType = reco::Muon::AllArbitrated;                    // all muons with the tracker muon arbitrated
    else if( muonID == "GlobalMuonPromptTight" )	
	idType = reco::Muon::GlobalMuonPromptTight;            // global muons with tighter fit requirements
    else if( muonID == "TMLastStationLoose" )	
	idType = reco::Muon::TMLastStationLoose;               // penetration depth loose selector
    else if( muonID == "TMLastStationTight" )	
	idType = reco::Muon::TMLastStationTight;               // penetration depth tight selector
    else if( muonID == "TM2DCompatibilityLoose" )	
	idType = reco::Muon::TM2DCompatibilityLoose;           // likelihood based loose selector
    else if( muonID == "TM2DCompatibilityTight" )	
	idType = reco::Muon::TM2DCompatibilityTight;           // likelihood based tight selector
    else if( muonID == "TMOneStationLoose" )	
	idType = reco::Muon::TMOneStationLoose;                // require one well matched segment
    else if( muonID == "TMOneStationTight" )	
	idType = reco::Muon::TMOneStationTight;                // require one well matched segment
    else if( muonID == "TMLastStationOptimizedLowPtLoose" )	
	idType = reco::Muon::TMLastStationOptimizedLowPtLoose; // combination of TMLastStation and TMOneStation
    else if( muonID == "TMLastStationOptimizedLowPtTight" )	
	idType = reco::Muon::TMLastStationOptimizedLowPtTight; // combination of TMLastStation and TMOneStation
    else
    {
    	LogDebug("MuonJetCrossCleaner") << "No match for given ID. Using 'All' muons";
	idType = reco::Muon::AllGlobalMuons;
    }

    muonID_ = idType;
}
