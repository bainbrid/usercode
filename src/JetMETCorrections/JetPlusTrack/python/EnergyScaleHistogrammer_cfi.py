import FWCore.ParameterSet.Config as cms

energyScaleHistogrammer = cms.Service(
    "EnergyScaleHistogrammer",
    RootFileName = cms.string('EnergyScale.root')
    )

