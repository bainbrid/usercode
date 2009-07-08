import FWCore.ParameterSet.Config as cms

process = cms.Process("EFF")

process.Timing =cms.Service("Timing")

process.Tracer = cms.Service(
    "Tracer",
    sourceSeed = cms.untracked.string("$$")
    )

process.load("DQM.SiStripCommon.MessageLogger_cfi")
process.MessageLogger.debugModules = cms.untracked.vstring('')

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.test = cms.EDAnalyzer(
    "TestEfficiency",
    OutputFile  = cms.untracked.string("result.root"),
    SignalFiles = cms.untracked.vstring(
    "GM1c_Jets0+.root"
    ),
    BkgdFiles   = cms.untracked.vstring(
    "PhotonJetsMadgraph200-Inf_Jets0+.root",
    "QCDMadgraph500-1000_Jets0+.root",
    ),
    Histograms  = cms.untracked.vstring(
    "test/AlphaT/AlphaT",
    "test/AlphaT/BetaT",
    "test/BiasedAlphaT/BiasedAlphaT",
    "test/BiasedAlphaT/BiasedBetaT",
    "test/Recoil/NewAlphaT",
    "test/Recoil/NewBetaT",
    #"test/Recoil/DPHI_BetweenJetSystemAndThrustAxis",
    ),
    )

process.p = cms.Path( 
    process.test
    )

