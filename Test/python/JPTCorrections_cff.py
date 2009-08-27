import FWCore.ParameterSet.Config as cms

# -------------------- Electron ID --------------------

from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *

eleIdRobustLoose = eidCutBasedExt.clone()
eleIdRobustLoose.electronIDType = cms.string('robust')
eleIdRobustLoose.electronQuality = cms.string('loose')

eleIdRobustTight = eidCutBasedExt.clone()
eleIdRobustTight.electronIDType = cms.string('robust')
eleIdRobustTight.electronQuality = cms.string('tight')

eleIdRobustHighEnergy = eidCutBasedExt.clone()
eleIdRobustHighEnergy.electronIDType = cms.string('robust')
eleIdRobustHighEnergy.electronQuality = cms.string('highenergy')

eleIdLoose = eidCutBasedExt.clone()
eleIdLoose.electronIDType = cms.string('classbased')
eleIdLoose.electronQuality = cms.string('loose')

eleIdTight = eidCutBasedExt.clone()
eleIdTight.electronIDType = cms.string('classbased')
eleIdTight.electronQuality = cms.string('tight')

eleIdSequence = cms.Sequence(
    eleIdRobustLoose +
    eleIdRobustTight +
    eleIdRobustHighEnergy +
    eleIdLoose +
    eleIdTight
    )

# -------------------- JetTrackAssociation --------------------

from JetMETCorrections.Configuration.JetCorrectionsRecord_cfi import *
from RecoJets.Configuration.RecoJetAssociations_cff import *
from RecoJets.JetAssociationProducers.iterativeCone5JTA_cff import*

JetTracksAssociationAtVertex = iterativeCone5JetTracksAssociatorAtVertex.clone() 
JetTracksAssociationAtVertex.jets = cms.InputTag("ZSPJetCorJetIcone5")

JetTracksAssociationAtCaloFace = iterativeCone5JetTracksAssociatorAtCaloFace.clone()
JetTracksAssociationAtCaloFace.jets = cms.InputTag("ZSPJetCorJetIcone5")

JetTracksAssociationExtender = iterativeCone5JetExtender.clone() 
JetTracksAssociationExtender.jets = cms.InputTag("ZSPJetCorJetIcone5")
JetTracksAssociationExtender.jet2TracksAtCALO = cms.InputTag("JetTracksAssociationAtCaloFace")
JetTracksAssociationExtender.jet2TracksAtVX = cms.InputTag("JetTracksAssociationAtVertex")

JetTrackAssociations = cms.Sequence(
    JetTracksAssociationAtVertex *
    JetTracksAssociationAtCaloFace *
    JetTracksAssociationExtender
    )

# -------------------- JPT Corrections --------------------

from bainbrid.Test.JPTCorrections_cfi import *

JPTCorrectionIC5 = cms.ESSource(
    "JPTCorrectionService",
    JPTCorrection,
    label = cms.string('JPTCorrectionIC5'),
    )

JPTCorrectorIC5 = cms.EDProducer("CaloJetCorrectionProducer",
    src = cms.InputTag("ZSPJetCorJetIcone5"),
    correctors = cms.vstring('JPTCorrectionIC5'),
    alias = cms.untracked.string('JPTCorrectorIC5')
)

# -------------------- Sequence --------------------

JPTCorrections = cms.Sequence(
    eleIdSequence *
    JetTrackAssociations *
    JPTCorrectorIC5
    )

