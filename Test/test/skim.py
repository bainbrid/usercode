import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring('file:INPUT.root')
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.p1 = cms.Path( process.content )

process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('skim.root'),
    outputCommands = cms.untracked.vstring(
    'drop *',
    'keep recoGenJets_iterativeCone5GenJets_*_*',
    'keep recoGenParticles_genParticles_*_*',
    'keep recoCaloJets_iterativeCone5CaloJets_*_*',
    'keep recoTracks_globalMuons_*_*'
    )
    )

process.e1 = cms.EndPath( process.output )
