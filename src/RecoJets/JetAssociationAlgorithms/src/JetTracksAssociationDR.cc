// Associate jets with tracks by simple "dR" criteria
// Fedor Ratnikov (UMd), Aug. 28, 2007
// $Id: $

#include "RecoJets/JetAssociationAlgorithms/interface/JetTracksAssociationDR.h"

// -----------------------------------------------------------------------------
//
JetTracksAssociationDR::JetTracksAssociationDR( double fDr ) 
  : mDeltaR2Threshold(fDr*fDr)
{;}

// -----------------------------------------------------------------------------
//
void JetTracksAssociationDR::associateTracksToJets( Association* fAssociation, 
						    const JetRefs& fJets,
						    const TrackRefs& fTracks ) 
{
  JetRefs::const_iterator ii = fJets.begin();
  JetRefs::const_iterator jj = fJets.end();
  for ( ; ii != jj; ++ii ) {
    reco::TrackRefVector associated;
    associateTracksToJet( associated, **ii, fTracks );
    reco::JetTracksAssociation::setValue( fAssociation, *ii, associated );
  }
}

// -----------------------------------------------------------------------------
//
void JetTracksAssociationDR::createJetRefs( const Jets& input, JetRefs& output ) {
  output.clear();
  output.reserve( input->size() );
  for ( unsigned ii = 0; ii < input->size(); ++ii ) { 
    output.push_back( input->refAt(ii) );
  }
}

// -----------------------------------------------------------------------------
//
void JetTracksAssociationDR::createTrackRefs( const Tracks& input, TrackRefs& output ) {
  output.clear();
  output.reserve( input->size() );
  for ( unsigned ii = 0; ii < input->size(); ++ii ) { 
    output.push_back( reco::TrackRef( input, ii ) );
  }
}
