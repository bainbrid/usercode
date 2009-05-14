import FWCore.ParameterSet.Config as cms

energyScaleAnalyzer = cms.EDAnalyzer(
    "EnergyScaleAnalyzer",
    GenObjectType  = cms.string('GenJet'),
    RecoObjectType = cms.string('CaloJet'),
    GenObjectTag   = cms.InputTag('iterativeCone5GenJets'),
    RecoObjectTaag = cms.InputTag('iterativeCone5CaloJets')
    )

#friendlyClassName, moduleLabel, productInstanceName and processName
