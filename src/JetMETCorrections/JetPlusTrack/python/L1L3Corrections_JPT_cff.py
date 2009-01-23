import FWCore.ParameterSet.Config as cms

#################### Level-1 "offset" corrections ##################################

L1JetCorrectorIC5ZSP = cms.ESSource(
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

L3JetCorrectorIC5JPT = cms.ESSource(
    "JetPlusTrackCorrectionService",
    label = cms.string('L3AbsoluteJetCorrectorJPT'),
    # Input files
    NonEfficiencyFile     = cms.string('CMSSW_167_TrackNonEff'),
    NonEfficiencyFileResp = cms.string('CMSSW_167_TrackLeakage'),
    ResponseFile          = cms.string('CMSSW_167_response'),
    # JetTrackAssociation objects built from ZSP-corrected ("L1CorJetIC5ZSP") jets 
    JetTrackCollectionAtVertex = cms.InputTag("L1CorJetIC5ZSPJetTracksAssociatorAtVertex"),
    JetTrackCollectionAtCalo   = cms.InputTag("L1CorJetIC5ZSPJetTracksAssociatorAtCaloFace"),
    # Muon collection
    muonSrc = cms.InputTag("globalMuons"),
    # Track quality definition
    TrackQuality = cms.string('highPurity'),
    UseQuality   = cms.bool(True),
    # Configuration of algorithm
    respalgo           = cms.int32(5),
    AddOutOfConeTracks = cms.bool(True)
    )

L3CorJetIC5JPT = cms.EDProducer(
    "CaloJetCorrectionProducer",
    src = cms.InputTag("L1CorJetIC5ZSP"),
    correctors = cms.vstring('L3AbsoluteJetCorrectorJPT'),
    verbose = cms.untracked.bool(False)
    )
