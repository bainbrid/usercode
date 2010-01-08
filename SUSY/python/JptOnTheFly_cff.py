import FWCore.ParameterSet.Config as cms

# On-the-fly ZSP and JPT corrections

from JetMETCorrections.Configuration.ZSPJetCorrections332_cff import *
from JetMETCorrections.Configuration.JetPlusTrackCorrections_cff import *

# Services and Producers for scalar JPT correction

ic5sJPTService = JetPlusTrackZSPCorrectorIcone5.clone()
sc5sJPTService = JetPlusTrackZSPCorrectorSiscone5.clone()
ak5sJPTService = JetPlusTrackZSPCorrectorAntiKt5.clone()

ic5sJPTService.VectorialCorrection = cms.bool(False)
sc5sJPTService.VectorialCorrection = cms.bool(False)
ak5sJPTService.VectorialCorrection = cms.bool(False)

ic5sJPTProducer = JetPlusTrackZSPCorJetIcone5.clone()
sc5sJPTProducer = JetPlusTrackZSPCorJetSiscone5.clone()
ak5sJPTProducer = JetPlusTrackZSPCorJetAntiKt5.clone()

ic5sJPTProducer.correctors = cms.vstring('ic5sJPTProducer')
sc5sJPTProducer.correctors = cms.vstring('sc5sJPTProducer')
ak5sJPTProducer.correctors = cms.vstring('ak5sJPTProducer')

# Services and Producers for vectorial JPT correction

ic5vJPTService = JetPlusTrackZSPCorrectorIcone5.clone()
sc5vJPTService = JetPlusTrackZSPCorrectorSiscone5.clone()
ak5vJPTService = JetPlusTrackZSPCorrectorAntiKt5.clone()

ic5vJPTService.VectorialCorrection = cms.bool(True)
sc5vJPTService.VectorialCorrection = cms.bool(True)
ak5vJPTService.VectorialCorrection = cms.bool(True)

ic5vJPTProducer = JetPlusTrackZSPCorJetIcone5.clone()
sc5vJPTProducer = JetPlusTrackZSPCorJetSiscone5.clone()
ak5vJPTProducer = JetPlusTrackZSPCorJetAntiKt5.clone()

ic5vJPTProducer.correctors = cms.vstring('ic5vJPTProducer')
sc5vJPTProducer.correctors = cms.vstring('sc5vJPTProducer')
ak5vJPTProducer.correctors = cms.vstring('ak5vJPTProducer')

# Sequences for ZSP corrections and jet-track association

ic5ZSP = cms.Sequence( 
    ZSPJetCorJetIcone5 *
    JPTeidTight *
    ZSPiterativeCone5JetTracksAssociatorAtVertex *
    ZSPiterativeCone5JetTracksAssociatorAtCaloFace *
    ZSPiterativeCone5JetExtender 
    )

sc5ZSP = cms.Sequence( 
    ZSPJetCorJetSiscone5 *
    JPTeidTight *
    ZSPSisCone5JetTracksAssociatorAtVertex *
    ZSPSisCone5JetTracksAssociatorAtCaloFace *
    ZSPSisCone5JetExtender 
    )

ak5ZSP = cms.Sequence( 
    ZSPJetCorJetAntiKt5 *
    JPTeidTight *
    ZSPAntiKt5JetTracksAssociatorAtVertex *
    ZSPAntiKt5JetTracksAssociatorAtCaloFace *
    ZSPAntiKt5JetExtender 
    )

# Sequences scalar and vectorial JPT corrections

ic5sJPT = cms.Sequence( ic5ZSP * ic5sJPTProducer )
sc5sJPT = cms.Sequence( ic5ZSP * sc5sJPTProducer )
ak5sJPT = cms.Sequence( ic5ZSP * ak5sJPTProducer )

ic5vJPT = cms.Sequence( ic5ZSP * ic5vJPTProducer )
sc5vJPT = cms.Sequence( ic5ZSP * sc5vJPTProducer )
ak5vJPT = cms.Sequence( ic5ZSP * ak5vJPTProducer )

sJPT = cms.Sequence( ic5sJPT + sc5sJPT + ak5sJPT )
vJPT = cms.Sequence( ic5vJPT + sc5vJPT + ak5vJPT )
JPT  = cms.Sequence( sJPT + vJPT )
