import FWCore.ParameterSet.Config as cms

process = cms.Process("DUMP")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring('file:./INPUT.root')
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
    )

process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.p1 = cms.Path( process.content )

