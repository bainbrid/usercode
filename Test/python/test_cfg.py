# TW modifications to analyse my Layer 1 PAT file...
import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    #'file:/home/bainbrid/data/susy/photonid/samples22X/CMSSW_2_2_5/src/production/signal/pat.root'
    #'file:/home/bainbrid/data/susy/photonid/samples22X/CMSSW_2_2_5/src/production/signal/loose/pat.root'
    'file:/home/bainbrid/data/susy/photonid/samples22X/CMSSW_2_2_5/src/production/signal/tight/pat.root'
    #'file:/home/bainbrid/data/susy/photonid/samples22X/CMSSW_2_2_5/src/production/background/pat.root'
    #'file:/home/bainbrid/data/susy/photonid/samples21X/CMSSW_2_2_5/src/production/background/pat.root'
    )
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.load("DQM.SiStripCommon.MessageLogger_cfi")
process.MessageLogger.debugModules = cms.untracked.vstring('')

process.load("PhysicsTools.PatAlgos.producersLayer1.photonIDProducer_cfi")

process.simplePhotonIDAnalysis = cms.EDFilter(
    "SimplePhotonIDAnalysis",
    electronTag = cms.untracked.InputTag("selectedLayer1Electrons"),
    tauTag      = cms.untracked.InputTag("selectedLayer1Taus"),
    muonTag     = cms.untracked.InputTag("selectedLayer1Muons"),
    jetTag      = cms.untracked.InputTag("selectedLayer1Jets"),
    photonTag   = cms.untracked.InputTag("selectedLayer1Photons"),
    photonIDTag = cms.untracked.InputTag("selectedLayer1PhotonIDs"),
    metTag      = cms.untracked.InputTag("selectedLayer1METs")
    )

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string('SimplePhotonIDAnalysis.root')
    )

process.p = cms.Path(
    process.selectedLayer1PhotonIDs *
    process.simplePhotonIDAnalysis
    )

