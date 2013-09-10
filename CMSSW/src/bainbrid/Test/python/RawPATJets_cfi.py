import FWCore.ParameterSet.Config as cms

RawPATJets = cms.EDProducer(
    "RawPATJetProducer",
    JetCollection = cms.InputTag("allLayer1JetsIC5"), 
    )
