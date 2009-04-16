import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    PATSummaryTables = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
    )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTBar-2_2_X_2008-11-03-STARTUP_V7-AODSIM.100.root')
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V11::All')

from PhysicsTools.PatAlgos.tools.jetTools import *

process.p = cms.Path(
    process.patDefaultSequence  
    )

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('JptCorrectedJets.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('drop *')
    )

process.outpath = cms.EndPath(process.out)

from PhysicsTools.PatAlgos.patEventContent_cff import *
process.out.outputCommands += patEventContent
process.out.outputCommands += ["keep *_selectedLayer1Jets*_*_*", "keep *_layer1METs*_*_*"]
