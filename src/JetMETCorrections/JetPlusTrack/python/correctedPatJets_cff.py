import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.L2L3Corrections_Summer08_cff import *
from JetMETCorrections.Configuration.L7PartonCorrections_cff import *
from JetMETCorrections.JetPlusTrack.L1L3Corrections_JPT_cff import *

from RecoJets.JetAssociationProducers.iterativeCone5JTA_cff import *

from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *
patJetCorrFactors = jetCorrFactors.clone(
    L1JetCorrector      = cms.string('L2RelativeJetCorrectorIC5ZSP'),
    L2JetCorrector      = cms.string('none'),
    L3JetCorrector      = cms.string('L3AbsoluteJetCorrectorIC5JPT'),
    L7udsJetCorrector   = cms.string('none'),
    L7gluonJetCorrector = cms.string('none'),
    L7cJetCorrector     = cms.string('none'),
    L7bJetCorrector     = cms.string('none')
    )

from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import * 
correctedPatJets = allLayer1Jets.clone(
    jetSource            = cms.InputTag("iterativeCone5CaloJets"),
    addJetCorrFactors    = cms.bool(False),
    jetCorrFactorsSource = cms.InputTag("patJetCorrFactors"),
    embedCaloTowers      = cms.bool(False),
    addResolutions       = cms.bool(False),
    useNNResolutions     = cms.bool(False),
    addBTagInfo          = cms.bool(False),
    addDiscriminators    = cms.bool(False),
    addTagInfoRefs       = cms.bool(False),
    addAssociatedTracks  = cms.bool(False),
    addJetCharge         = cms.bool(False),
    addTrigMatch         = cms.bool(False),
    addGenPartonMatch    = cms.bool(False),
    embedGenPartonMatch  = cms.bool(False),
    addGenJetMatch       = cms.bool(False),
    addPartonJetMatch    = cms.bool(False),
    getJetMCFlavour      = cms.bool(False),
    addEfficiencies      = cms.bool(False)
    )

correctedPatJetSeq = cms.Sequence(
    iterativeCone5JTA *
    patJetCorrFactors *
    correctedPatJets
    )

keepCorrectedPatJets = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoGenJets_iterativeCone5GenJets_*_*',
                                           'keep recoCaloJets_iterativeCone5CaloJets_*_*',
                                           'keep patJets_correctedPatJets_*_*',
                                           'keep patJetCorrFactorsedmValueMap_patJetCorrFactors_*_*'
                                           )
    )


