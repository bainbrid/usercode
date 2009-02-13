import FWCore.ParameterSet.Config as cms

from JetMETCorrections.JetPlusTrack.L1L3Corrections_JPT_cff import *

correctedCaloJetSeq = cms.Sequence(
    L1CorJetIC5ZSP *
    L3CorJetIC5JPT
    )

keepCorrectedCaloJets = cms.PSet(
    outputCommands = cms.untracked.vstring(
    'keep recoCaloJets_L1CorJetIC5ZSP_*_*',
    'keep recoCaloJets_L3CorJetIC5JPT_*_*'
    )
    )

from JetMETCorrections.JetPlusTrack.EnergyScaleHistogrammer_cfi import * 
energyScaleHistogrammer.RootFileName = 'CaloJetEnergyScale.root'

from JetMETCorrections.JetPlusTrack.EnergyScaleAnalyzer_cfi import * 

rawCaloHistos = energyScaleAnalyzer.clone(
    GenObjectType  = cms.string('GenJet'),
    RecoObjectType = cms.string('CaloJet'),
    GenObjectTag   = cms.InputTag('iterativeCone5GenJetsNoNuBSM'),
    RecoObjectTag  = cms.InputTag('iterativeCone5CaloJets'),
    verbose = cms.untracked.bool(False)
    )

zspCaloHistos = energyScaleAnalyzer.clone(
    GenObjectType  = cms.string('GenJet'),
    RecoObjectType = cms.string('CaloJet'),
    GenObjectTag   = cms.InputTag('iterativeCone5GenJetsNoNuBSM'),
    RecoObjectTag  = cms.InputTag('L1CorJetIC5ZSP'),
    verbose = cms.untracked.bool(False)
    )

jptCaloHistos = energyScaleAnalyzer.clone(
    GenObjectType  = cms.string('GenJet'),
    RecoObjectType = cms.string('CaloJet'),
    GenObjectTag   = cms.InputTag('iterativeCone5GenJetsNoNuBSM'),
    RecoObjectTag  = cms.InputTag('L3CorJetIC5JPT'),
    verbose = cms.untracked.bool(False)
    )

correctedCaloJetHistos = cms.Sequence(
    rawCaloHistos +
    zspCaloHistos +
    jptCaloHistos
    )
