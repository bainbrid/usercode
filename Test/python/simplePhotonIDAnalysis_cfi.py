import FWCore.ParameterSet.Config as cms

simplePhotonIDAnalysis = cms.EDFilter(
    "SimplePhotonIDAnalysis",
    Photons      = cms.InputTag("patPhotonIDProducer"),
    OtherPhotons = cms.InputTag("selectedLayer1Photons"),
    GenParticles = cms.InputTag("genParticles"),
    PhotonsPdgIds = cms.vint32(),
    LeptonsPdgIds = cms.vint32(),
    ParticlesPdgIds = cms.vint32(),
    )
