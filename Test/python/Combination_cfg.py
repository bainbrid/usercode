import FWCore.ParameterSet.Config as cms

process = cms.Process("COMBINATION")

process.load("DQM.SiStripCommon.MessageLogger_cfi")
process.MessageLogger.debugModules = cms.untracked.vstring('')

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("file:combinations.root")
    )

process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring('file:pat.root')
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    )

process.anal = cms.EDAnalyzer("EventContentAnalyzer")

process.test = cms.EDAnalyzer(
    "TestCombination",
    # Misc
    MaximumObjects = cms.int32(30),
    TestObjects    = cms.int32(-1),
    # Collections
    Photons        = cms.InputTag("selectedLayer1Photons"),
    Jets           = cms.InputTag("selectedLayer1Jets"),
    Muons          = cms.InputTag("selectedLayer1Muons"),
    Electrons      = cms.InputTag("selectedLayer1Electrons"),
    # Object selection
    PhotonEt       = cms.double(50.),
    PhotonEta      = cms.double(2.4),
    JetEt          = cms.double(50.),
    JetEta         = cms.double(2.4),
    JetEMfraction  = cms.double(0.9),
    MuonPt         = cms.double(10.),
    MuonEta        = cms.double(2.4),
    MuonTrkIso     = cms.double(10.0),
    ElectronPt     = cms.double(10.),
    ElectronEta    = cms.double(2.4),
    ElectronTrkIso = cms.double(1.0),
    # Event selection
    TotalEt        = cms.double(350.),
    )

process.p1 = cms.Path( process.test )

