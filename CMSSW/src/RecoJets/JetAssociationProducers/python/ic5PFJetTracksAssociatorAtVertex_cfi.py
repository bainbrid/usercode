import FWCore.ParameterSet.Config as cms

from RecoJets.JetAssociationProducers.j2tParametersVX_cfi import *
ic5PFJetTracksAssociatorAtVertex = cms.EDFilter("JetTracksAssociatorAtVertex",
    j2tParametersVX,
    jets = cms.InputTag("iterativeCone5PFJets")
)


