// Associate jets with tracks by simple "dR" criteria
// Fedor Ratnikov (UMd), Aug. 28, 2007
// $Id: JetTracksAssociationDRVertex.cc,v 1.2 2007/09/19 18:02:40 fedor Exp $

#include "RecoJets/JetAssociationAlgorithms/interface/JetTracksAssociationDRVertex.h"
#include "DataFormats/Math/interface/deltaR.h"

// -----------------------------------------------------------------------------
//
JetTracksAssociationDRVertex::JetTracksAssociationDRVertex( double fDr ) 
  : JetTracksAssociationDR(fDr),
    propagatedTracks_()
{;}

// -----------------------------------------------------------------------------
//
void JetTracksAssociationDRVertex::produce( Association* fAssociation, 
					    const Jets& fJets,
					    const Tracks& fTracks ) 
{
  JetRefs jets;
  createJetRefs( fJets, jets );
  TrackRefs tracks;
  createTrackRefs( fTracks, tracks );
  produce( fAssociation, jets, tracks );
}

// -----------------------------------------------------------------------------
//
void JetTracksAssociationDRVertex::produce( Association* fAssociation, 
					    const JetRefs& fJets,
					    const TrackRefs& fTracks ) 
{
  //clear();
  propagateTracks( fTracks ); 
  associateTracksToJets( fAssociation, fJets, fTracks ); 
}

// -----------------------------------------------------------------------------
//
void JetTracksAssociationDRVertex::associateTracksToJet( reco::TrackRefVector& associated,
							 const reco::Jet& fJet,
							 const TrackRefs& fTracks ) 
{
  associated.clear();
  std::vector<math::RhoEtaPhiVector>::const_iterator ii = propagatedTracks_.begin();
  std::vector<math::RhoEtaPhiVector>::const_iterator jj = propagatedTracks_.end();
  for ( ; ii != jj; ++ii ) {
    uint32_t index = ii - propagatedTracks_.begin();
    double dR2 = deltaR2( fJet.eta(), fJet.phi(), ii->eta(), ii->phi() );
    if ( dR2 < mDeltaR2Threshold ) { associated.push_back( fTracks[index] ); }
  }
}

// -----------------------------------------------------------------------------
//
void JetTracksAssociationDRVertex::propagateTracks( const TrackRefs& fTracks ) 
{
  propagatedTracks_.clear();
  propagatedTracks_.reserve( fTracks.size() );
  TrackRefs::const_iterator ii = fTracks.begin();
  TrackRefs::const_iterator jj = fTracks.end();
  for ( ; ii != jj; ++ii ) {
    propagatedTracks_.push_back( math::RhoEtaPhiVector( (**ii).p(), (**ii).eta(), (**ii).phi() ) ); 
  }
}  


