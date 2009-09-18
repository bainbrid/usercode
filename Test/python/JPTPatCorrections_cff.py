import FWCore.ParameterSet.Config as cms

from bainbrid.Test.JPTCorrections_cff import *

# -------------------- ZSP Corrections on PAT --------------------

ZSPPatJetCorrectorIC5 = cms.ESSource(
    "ZSPJetCorrectionService",
    tagName = cms.string('ZSP_CMSSW219_Iterative_Cone_05'),
    label = cms.string('ZSPJetCorrectorIcone5')
    )

ZSPPatJetCorrectionIC5 = cms.EDProducer(
    "PatJetCorrectionProducer",
    src = cms.InputTag("cleanLayer1Jets"),
    correctors = cms.vstring('ZSPPatJetCorrectorIcone5'),
    alias = cms.untracked.string('ZSPPatJetCorrectionIC5')
    )

# -------------------- JTA on PAT --------------------

PatJetTracksAssociationAtVertex = iterativeCone5JetTracksAssociatorAtVertex.clone() 
PatJetTracksAssociationAtVertex.jets = cms.InputTag("ZSPPatJetCorrectionIC5")

PatJetTracksAssociationAtCaloFace = iterativeCone5JetTracksAssociatorAtCaloFace.clone()
PatJetTracksAssociationAtCaloFace.jets = cms.InputTag("ZSPPatJetCorrectionIC5")

PatJetTracksAssociationExtender = iterativeCone5JetExtender.clone() 
PatJetTracksAssociationExtender.jets = cms.InputTag("ZSPPatJetCorrectionIC5")
PatJetTracksAssociationExtender.jet2TracksAtCALO = cms.InputTag("PatJetTracksAssociationAtCaloFace")
PatJetTracksAssociationExtender.jet2TracksAtVX = cms.InputTag("PatJetTracksAssociationAtVertex")

PatJetTrackAssociations = cms.Sequence(
    PatJetTracksAssociationAtVertex *
    PatJetTracksAssociationAtCaloFace *
    PatJetTracksAssociationExtender
    )

# -------------------- JPT Corrections on PAT --------------------

JPTPatCorrectorIC5 = cms.ESSource(
    "JPTPatCorrectionService",
    # Import JPT configuration
    JPTCorrection,
    # JTA on-the-fly
    AllowOnTheFly     = cms.bool(True),
    Tracks            = cms.InputTag("generalTracks"),
    Propagator        = cms.string('SteppingHelixPropagatorAlong'),
    ConeSize          = cms.double(0.5),
    # Misc
    UsePatCollections = cms.bool(True),
    label = cms.string('JPTPatCorrectorIC5'),
    )

JPTPatCorrectionIC5 = cms.EDProducer("CaloJetCorrectionProducer",
    src = cms.InputTag("ZSPJetCorJetIcone5"),
    correctors = cms.vstring('JPTPatCorrectorIC5'),
    alias = cms.untracked.string('JPTPatCorrectionIC5')
)

# -------------------- Sequences --------------------

JPTPatCorrections = cms.Sequence(
    eleIdSequence *
    ZSPPatJetCorrectionIC5 *
    PatJetTrackAssociations *
    JPTPatCorrectionIC5
    )

