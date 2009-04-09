import FWCore.ParameterSet.Config as cms

simplePhotonIDAnalysis = cms.EDFilter(
    "SimplePhotonIDAnalysis",
    Photons      = cms.untracked.InputTag("patPhotonIDProducer"),
    OtherPhotons = cms.untracked.InputTag("selectedLayer1Photons"),
    )

