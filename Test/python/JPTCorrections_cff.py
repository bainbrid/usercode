import FWCore.ParameterSet.Config as cms

# Jet-Tracks Association

from JetMETCorrections.Configuration.JetCorrectionsRecord_cfi import *
from RecoJets.Configuration.RecoJetAssociations_cff import *
from RecoJets.JetAssociationProducers.iterativeCone5JTA_cff import*

# -------------------- Electron ID --------------------

from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *

electronIdRobustLoose = eidCutBasedExt.clone()
electronIdRobustLoose.electronQuality = 'robust'

electronIdRobustTight = eidCutBasedExt.clone()
electronIdRobustTight.electronQuality = 'robust'
electronIdRobustTight.robustEleIDCuts.barrel = [0.015, 0.0092, 0.020, 0.0025]
electronIdRobustTight.robustEleIDCuts.endcap = [0.018, 0.025, 0.020, 0.0040]

electronIdRobustHighEnergy = eidCutBasedExt.clone()
electronIdRobustHighEnergy.electronQuality = 'robust'
electronIdRobustHighEnergy.robustEleIDCuts.barrel = [0.050, 0.011, 0.090, 0.005]
electronIdRobustHighEnergy.robustEleIDCuts.endcap = [0.100, 0.0275, 0.090, 0.007]

electronIdLoose = eidCutBasedExt.clone()
electronIdLoose.electronQuality = 'loose'

electronIdTight = eidCutBasedExt.clone()
electronIdTight.electronQuality = 'loose'

electronIdSequence = cms.Sequence(
    electronIdRobustLoose +
    electronIdRobustTight +
    electronIdRobustHighEnergy +
    electronIdLoose +
    electronIdTight
    )

# -------------------- JetTrackAssociation --------------------

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

# Sequence

JPTCorrections = cms.Sequence(
    electronIdSequence *
    JetTrackAssociations *
    JPTCorrectorIC5
    )

