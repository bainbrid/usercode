import FWCore.ParameterSet.Config as cms

# Correctors and producers for ZSP correction

ak5ZSPCorrector = cms.ESSource(
    "ZSPJetCorrectionService",
    tagName = cms.vstring('ZSP_CMSSW332_Iterative_Cone_05_PU0'),
    tagNameOffset = cms.vstring(),
    label = cms.string('ak5ZSPCorrector'),
    PU = cms.int32(-1),
    FixedPU = cms.int32(0),
    )
 
ak5ZSPProducer = cms.EDProducer(
    "CaloJetCorrectionProducer",
    src = cms.InputTag("ak5CaloJets"),
    correctors = cms.vstring('ak5ZSPCorrector'),
    alias = cms.untracked.string('ak5ZSPProducer'),
    )

# Jet-track association using AK5 ZSP-corrected jet collection as input

from RecoJets.JetAssociationProducers.ak5JTA_cff import*

ak5ZSPAtVertex = ak5JetTracksAssociatorAtVertex.clone()
ak5ZSPAtVertex.jets = cms.InputTag("ak5ZSPProducer")

ak5ZSPAtCaloFace = ak5JetTracksAssociatorAtCaloFace.clone()
ak5ZSPAtCaloFace.jets = cms.InputTag("ak5ZSPProducer")

ak5ZSPExtender = ak5JetExtender.clone()
ak5ZSPExtender.jets = cms.InputTag("ak5ZSPProducer")
ak5ZSPExtender.jet2TracksAtCALO = cms.InputTag("ak5ZSPAtVertex")
ak5ZSPExtender.jet2TracksAtVX = cms.InputTag("ak5ZSPAtCaloFace")

# Correctors and Producers for scalar JPT correction

from JetMETCorrections.Configuration.JetPlusTrackCorrections_cfi import *

ic5sJPTCorrector = cms.ESSource(
    "JetPlusTrackCorrectionService",
    cms.PSet(JPTZSPCorrectorICone5),
    label = cms.string('ic5sJPTCorrector'),
    )
ic5sJPTCorrector.JetTracksAssociationAtVertex = cms.InputTag("ak5ZSPAtVertex")
ic5sJPTCorrector.JetTracksAssociationAtCaloFace = cms.InputTag("ak5ZSPAtCaloFace")
ic5sJPTCorrector.JetSplitMerge = cms.int32(2)
ic5sJPTCorrector.ElectronIds = cms.InputTag("eidTight")
ic5sJPTCorrector.VectorialCorrection = cms.bool(False)

ak5sJPTProducer = cms.EDProducer(
    "CaloJetCorrectionProducer",
    src = cms.InputTag("ak5ZSPProducer"),
    correctors = cms.vstring('ic5sJPTCorrector'),
    alias = cms.untracked.string('ak5sJPTProducer'),
    )

# Correctors and Producers for vectorial JPT correction

ic5vJPTCorrector = cms.ESSource(
    "JetPlusTrackCorrectionService",
    cms.PSet(JPTZSPCorrectorICone5),
    label = cms.string('ic5vJPTCorrector'),
    )
ic5vJPTCorrector.JetTracksAssociationAtVertex = cms.InputTag("ak5ZSPAtVertex")
ic5vJPTCorrector.JetTracksAssociationAtCaloFace = cms.InputTag("ak5ZSPAtCaloFace")
ic5vJPTCorrector.JetSplitMerge = cms.int32(2)
ic5vJPTCorrector.ElectronIds = cms.InputTag("eidTight")
ic5vJPTCorrector.VectorialCorrection = cms.bool(True)
ic5vJPTCorrector.UseResponseInVecCorr = cms.bool(False)

ak5vJPTProducer = cms.EDProducer(
    "CaloJetCorrectionProducer",
    src = cms.InputTag("ak5ZSPProducer"),
    correctors = cms.vstring('ic5vJPTCorrector'),
    alias = cms.untracked.string('ak5vJPTProducer'),
    )

# Sequences for ZSP, JTA and JPT

ak5ZSP = cms.Sequence( 
    ak5ZSPProducer *
    ak5ZSPAtVertex *
    ak5ZSPAtCaloFace *
    ak5ZSPExtender
    )

ak5JPT = cms.Sequence( ak5sJPTProducer + ak5vJPTProducer )

JPT = cms.Sequence( ak5ZSP * ak5JPT )

