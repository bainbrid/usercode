import FWCore.ParameterSet.Config as cms

JPTCorrFactors = cms.EDProducer(
    "JPTCorrFactorsProducer",
    UncorrectedJets = cms.InputTag("allLayer1JetsIC5"), 
    CorrectedJets   = cms.InputTag("allLayer1JetsIC5JPT"),
    )
