import FWCore.ParameterSet.Config as cms

# Ntuple entries for scalar JPT corrections to jet collections

susycafic5sjptjet = cms.EDProducer(
    "SusyCAF_PatJet",
    InputTag = cms.InputTag('ic5sJPTProducer'),
    PrimaryVertexTag = cms.InputTag('offlinePrimaryVertices'),
    MaxD0Trk = cms.double(0.02),
    PtErrFracTrk = cms.double(0.2),
    Prefix = cms.string('ic5sJpt'),
    Suffix = cms.string('Pat')
    )

susycafsc5sjptjet = cms.EDProducer(
    "SusyCAF_PatJet",
    InputTag = cms.InputTag('sc5sJPTProducer'),
    PrimaryVertexTag = cms.InputTag('offlinePrimaryVertices'),
    MaxD0Trk = cms.double(0.02),
    PtErrFracTrk = cms.double(0.2),
    Prefix = cms.string('sc5sJpt'),
    Suffix = cms.string('Pat')
    )

susycafak5sjptjet = cms.EDProducer(
    "SusyCAF_PatJet",
    InputTag = cms.InputTag('ak5sJPTProducer'),
    PrimaryVertexTag = cms.InputTag('offlinePrimaryVertices'),
    MaxD0Trk = cms.double(0.02),
    PtErrFracTrk = cms.double(0.2),
    Prefix = cms.string('ak5sJpt'),
    Suffix = cms.string('Pat')
    )

# Ntuple entries for vectorial JPT corrections to jet collections

susycafic5vjptjet = cms.EDProducer(
    "SusyCAF_PatJet",
    InputTag = cms.InputTag('ic5vJPTProducer'),
    PrimaryVertexTag = cms.InputTag('offlinePrimaryVertices'),
    MaxD0Trk = cms.double(0.02),
    PtErrFracTrk = cms.double(0.2),
    Prefix = cms.string('ic5vJpt'),
    Suffix = cms.string('Pat')
    )

susycafsc5vjptjet = cms.EDProducer(
    "SusyCAF_PatJet",
    InputTag = cms.InputTag('sc5vJPTProducer'),
    PrimaryVertexTag = cms.InputTag('offlinePrimaryVertices'),
    MaxD0Trk = cms.double(0.02),
    PtErrFracTrk = cms.double(0.2),
    Prefix = cms.string('sc5vJpt'),
    Suffix = cms.string('Pat')
    )

susycafak5vjptjet = cms.EDProducer(
    "SusyCAF_PatJet",
    InputTag = cms.InputTag('ak5vJPTProducer'),
    PrimaryVertexTag = cms.InputTag('offlinePrimaryVertices'),
    MaxD0Trk = cms.double(0.02),
    PtErrFracTrk = cms.double(0.2),
    Prefix = cms.string('ak5vJpt'),
    Suffix = cms.string('Pat')
    )

