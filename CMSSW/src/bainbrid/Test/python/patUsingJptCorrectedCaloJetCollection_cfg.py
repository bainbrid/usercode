import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATLayer0Summary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    PATLayer0Summary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
    )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTBar-2_2_X_2008-11-03-STARTUP_V7-AODSIM.100.root')
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V9::All')

process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08Redigi_cff")
process.load("JetMETCorrections.Configuration.JetPlusTrackCorrections_cff")
process.load("JetMETCorrections.Configuration.ZSPJetCorrections219_cff")
process.prefer("L2L3JetCorrectorIC5JPT") 
process.s = cms.Sequence( process.ZSPJetCorrections * process.JetPlusTrackCorrections * process.L2L3CorJetIC5JPT )

from PhysicsTools.PatAlgos.tools.jetTools import *

addJetCollection(process,
                 'JetPlusTrackZSPCorJetIcone5',
                 'JPT', 
                 runCleaner='CaloJet',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=None,
                 doType1MET=False,
                 doL1Counters=False
                 )

process.p = cms.Path(
    process.s *
    process.patLayer0 *
    process.patLayer1
    )

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('JptCorrectedJets.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('drop *')
    )

process.outpath = cms.EndPath(process.out)

process.load("PhysicsTools.PatAlgos.patLayer1_EventContent_cff")
process.out.outputCommands.extend(process.patLayer1EventContent.outputCommands)
process.out.outputCommands.extend(["keep *_selectedLayer1Jets*_*_*"])
