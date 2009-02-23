import FWCore.ParameterSet.Config as cms

#################### Level-1 "offset" corrections ##################################

L1JetCorrectorZSP = cms.ESSource(
    "ZSPJetCorrectionService",
    label = cms.string('L1OffsetJetCorrectorZSP'),
    tagName = cms.string('ZSP_CMSSW219_Iterative_Cone_05')
    )

L1CorJetIC5ZSP = cms.EDProducer(
    "CaloJetCorrectionProducer",
    src = cms.InputTag("iterativeCone5CaloJets"),
    correctors = cms.vstring('L1OffsetJetCorrectorZSP'),
    verbose = cms.untracked.bool(False)
    )

#################### Level-3 "absolute" corrections ##################################

# Required by JPT algorithm
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi import * 

L3JetCorrectorJPT = cms.ESSource(
    "JetPlusTrackCorrectionService",
    label = cms.string('L3AbsoluteJetCorrectorJPT'),
    # Look-up tables
    NonEfficiencyFile     = cms.string('CMSSW_167_TrackNonEff'),
    NonEfficiencyFileResp = cms.string('CMSSW_167_TrackLeakage'),
    ResponseFile          = cms.string('CMSSW_167_response'),
    # Access to tracks and muons
    trackSrc     = cms.InputTag("generalTracks"),
    UseQuality   = cms.bool(True),
    TrackQuality = cms.string('highPurity'),
    muonSrc      = cms.InputTag("globalMuons"),
    # Jet-tracks association (null InputTags mean use "on-the-fly")
    JetTrackCollectionAtVertex = cms.InputTag(""), 
    JetTrackCollectionAtCalo   = cms.InputTag(""),
    Propagator = cms.string('SteppingHelixPropagatorAlong'),
    coneSize   = cms.double(0.5),
    # Other misc algo configuration
    respalgo           = cms.int32(5),
    AddOutOfConeTracks = cms.bool(True)
    )

L3CorJetIC5JPT = cms.EDProducer(
    "CaloJetCorrectionProducer",
    src = cms.InputTag("L1CorJetIC5ZSP"),
    correctors = cms.vstring('L3AbsoluteJetCorrectorJPT'),
    verbose = cms.untracked.bool(False)
    )


#################### L1+L3 corrections ##################################

L1L3CorJetIC5JPT = cms.EDProducer(
    "CaloJetCorrectionProducer",
    src = cms.InputTag("iterativeCone5CaloJets"),
    correctors = cms.vstring('L1OffsetJetCorrectorZSP',
                             'L3AbsoluteJetCorrectorJPT'),
    verbose = cms.untracked.bool(False)
    )

