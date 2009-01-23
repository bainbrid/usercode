import FWCore.ParameterSet.Config as cms

from JetMETCorrections.JetPlusTrack.L1L3Corrections_JPT_cff import *

from RecoJets.JetAssociationProducers.iterativeCone5JTA_cff import *
L1CorJetIC5ZSPJetTracksAssociatorAtVertex = iterativeCone5JetTracksAssociatorAtVertex.clone() 
L1CorJetIC5ZSPJetTracksAssociatorAtVertex.jets = cms.InputTag("L1CorJetIC5ZSP")
L1CorJetIC5ZSPJetTracksAssociatorAtCaloFace = iterativeCone5JetTracksAssociatorAtCaloFace.clone()
L1CorJetIC5ZSPJetTracksAssociatorAtCaloFace.jets = cms.InputTag("L1CorJetIC5ZSP")
L1CorJetIC5ZSPJetExtender = iterativeCone5JetExtender.clone() 
L1CorJetIC5ZSPJetExtender.jets = cms.InputTag("L1CorJetIC5ZSP")
L1CorJetIC5ZSPJetExtender.jet2TracksAtVX = cms.InputTag("L1CorJetIC5ZSPJetTracksAssociatorAtVertex")
L1CorJetIC5ZSPJetExtender.jet2TracksAtCALO = cms.InputTag("L1CorJetIC5ZSPJetTracksAssociatorAtCaloFace")
L1CorJetIC5ZSPJTA = cms.Sequence(
    L1CorJetIC5ZSPJetTracksAssociatorAtVertex *
    L1CorJetIC5ZSPJetTracksAssociatorAtCaloFace *
    L1CorJetIC5ZSPJetExtender
    )

# this works
correctedCaloJetSeq = cms.Sequence( L1CorJetIC5ZSP * L1CorJetIC5ZSPJTA * L3CorJetIC5JPT )

# Useful products
keepCorrectedCaloJets = cms.PSet(
    outputCommands = cms.untracked.vstring(
    'keep recoGenJets_iterativeCone5GenJets_*_*',
    'keep recoCaloJets_iterativeCone5CaloJets_*_*',
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
    GenObjectTag   = cms.InputTag('iterativeCone5GenJets'),
    RecoObjectTag  = cms.InputTag('iterativeCone5CaloJets')
    )

zspCaloHistos = energyScaleAnalyzer.clone(
    GenObjectType  = cms.string('GenJet'),
    RecoObjectType = cms.string('CaloJet'),
    GenObjectTag   = cms.InputTag('iterativeCone5GenJets'),
    RecoObjectTag  = cms.InputTag('L1CorJetIC5ZSP')
    )

jptCaloHistos = energyScaleAnalyzer.clone(
    GenObjectType  = cms.string('GenJet'),
    RecoObjectType = cms.string('CaloJet'),
    GenObjectTag   = cms.InputTag('iterativeCone5GenJets'),
    RecoObjectTag  = cms.InputTag('L3CorJetIC5JPT')
    )

correctedCaloJetHistos = cms.Sequence(
    rawCaloHistos +
    zspCaloHistos +
    jptCaloHistos
    )
