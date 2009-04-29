import FWCore.ParameterSet.Config as cms

process = cms.Process("COMBINATION")

process.load("DQM.SiStripCommon.MessageLogger_cfi")

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("file:test.root")
    )

process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring('file:pat.root')
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.test = cms.EDAnalyzer(
    "TestCombination",
    JetCollection    = cms.InputTag("selectedLayer1Jets"),
    PhotonCollection = cms.InputTag("selectedLayer1Photons"),
    JetPt            = cms.double(50.),
    JetEta           = cms.double(2.4),
    PhotonPt         = cms.double(30.),
    PhotonEta        = cms.double(2.4),
    TotalPt          = cms.double(300.),
    )

process.p1 = cms.Path( process.test )

