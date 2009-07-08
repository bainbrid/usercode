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

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("eff.root")
    )

process.test = cms.EDAnalyzer(
    "TestMerge",
    OutputFile = cms.untracked.string("QCDMadgraph500-1000_Jets0+.root"),
    InputFiles = cms.untracked.VPSet(
    cms.PSet(
    FileName = cms.untracked.string("/home/bainbrid/susy/crab/alphaT/090606/QCDMadgraph500-1000/Jets0+/QCDMadgraph500-1000_Jets0+.root"),
    XSection = cms.untracked.double(4200),
    ),
    ),
    NormalisedLumi = cms.untracked.double(100.),
    )

process.p = cms.Path( 
    process.test
    )

