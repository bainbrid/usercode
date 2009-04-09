import FWCore.ParameterSet.Config as cms

simplePhotonIDAnalysis = cms.EDFilter(
    "SimplePhotonIDAnalysis",
    PhotonsWithOldID = cms.untracked.InputTag("selectedLayer1Photons"),
    PhotonsWithNewID = cms.untracked.InputTag("patPhotonIDProducer"),
    Jets             = cms.untracked.InputTag("selectedLayer1Jets"),
    )

