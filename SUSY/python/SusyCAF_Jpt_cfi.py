import FWCore.ParameterSet.Config as cms

# Scalar JPT corrections to jet collections

susycafak5sjptjetreco = cms.EDProducer(
    "SusyCAF_CaloJet",
    InputTag = cms.InputTag('ak5sJPTProducer'),
    PrimaryVertexTag = cms.InputTag('offlinePrimaryVertices'),
    MaxD0Trk = cms.double(0.02),
    PtErrFracTrk = cms.double(0.2),
    Prefix = cms.string('ak5sJpt'),
    Suffix = cms.string('Calo')
    )

# Vectorial JPT corrections to jet collections

susycafak5vjptjetreco = cms.EDProducer(
    "SusyCAF_CaloJet",
    InputTag = cms.InputTag('ak5vJPTProducer'),
    PrimaryVertexTag = cms.InputTag('offlinePrimaryVertices'),
    MaxD0Trk = cms.double(0.02),
    PtErrFracTrk = cms.double(0.2),
    Prefix = cms.string('ak5vJpt'),
    Suffix = cms.string('Calo')
    )

