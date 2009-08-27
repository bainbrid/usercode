import FWCore.ParameterSet.Config as cms

# -------------------- Electron ID --------------------

from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *

eidRobustLoose = eidCutBasedExt.clone()
eidRobustLoose.electronQuality = 'robust'

eidRobustTight = eidCutBasedExt.clone()
eidRobustTight.electronQuality = 'robust'
eidRobustTight.robusttightEleIDCuts.barrel = [0.015, 0.0092, 0.020, 0.0025]
eidRobustTight.robusttightEleIDCuts.endcap = [0.018, 0.025, 0.020, 0.0040]

eidRobustHighEnergy = eidCutBasedExt.clone()
eidRobustHighEnergy.electronQuality = 'robust'
eidRobustHighEnergy.robusthighenergyEleIDCuts.barrel = [0.050, 0.011, 0.090, 0.005]
eidRobustHighEnergy.robusthighenergyEleIDCuts.endcap = [0.100, 0.0275, 0.090, 0.007]

eidLoose = eidCutBasedExt.clone()
eidLoose.electronQuality = 'loose'

eidTight = eidCutBasedExt.clone()
eidTight.electronQuality = 'loose'

eidSequence = cms.Sequence(
    eidRobustLoose +
    eidRobustTight +
    eidRobustHighEnergy +
    eidLoose +
    eidTight
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
    eidSequence *
    JetTrackAssociations *
    JPTCorrectorIC5
    )

