import FWCore.ParameterSet.Config as cms

process = cms.Process("VerySimplePATAnalysis")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:myPATfile.root')
)

process.MessageLogger = cms.Service("MessageLogger")

process.verySimplePATAnalysis = cms.EDFilter("VerySimplePATAnalysis",
    electronTag = cms.untracked.InputTag("selectedLayer1Electrons"),
    tauTag      = cms.untracked.InputTag("selectedLayer1Taus"),
    muonTag     = cms.untracked.InputTag("selectedLayer1Muons"),
    jetTag      = cms.untracked.InputTag("selectedLayer1Jets"),
    photonTag   = cms.untracked.InputTag("selectedLayer1Photons"),
    metTag      = cms.untracked.InputTag("selectedLayer1METs")
)

process.TFileService = cms.Service("TFileService", fileName = cms.string('histo.root') )

process.p = cms.Path(process.verySimplePATAnalysis)

