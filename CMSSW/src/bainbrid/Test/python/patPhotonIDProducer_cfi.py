import FWCore.ParameterSet.Config as cms

patPhotonIDProducer = cms.EDProducer(
    "PATPhotonIDProducer",
    Photons = cms.InputTag("selectedLayer1Photons"),

    photonProducer = cms.string('photons'),                              
    photonLabel = cms.string(''),
    photonIDAssociationLabel = cms.string('PhotonAssociatedID'),
    photonIDLabel = cms.string('PhotonIDCutBasedProducer'),
    barrelEcalRecHitProducer = cms.string('ecalRecHit'),
    barrelEcalRecHitCollection = cms.string('EcalRecHitsEB'),
    endcapEcalRecHitProducer = cms.string('ecalRecHit'),
    endcapEcalRecHitCollection = cms.string('EcalRecHitsEE'),
    HcalRecHitProducer = cms.string('hbhereco'),
    HcalRecHitCollection = cms.string(''),
    GsfRecoCollection = cms.InputTag("pixelMatchGsfElectrons"),
    modulePhiBoundary =   cms.double(0.0087),
    moduleEtaBoundary = cms.vdouble(0.0, 0.05, 0.4, 0.5, 0.75, 0.85, 1.1, 1.2, 1.4, 1.6),
    trackProducer = cms.InputTag("generalTracks"),
    doCutBased = cms.bool(True),

    RequireNotElectron = cms.bool(False),                    
    RequireFiducial = cms.bool(False),
    DoHollowConeTrackIsolationCut = cms.bool(True),
    DoSolidConeTrackIsolationCut = cms.bool(False),
    DoHollowConeNTrkCut = cms.bool(False),
    DoSolidConeNTrkCut = cms.bool(False),
    DoHadOverEMCut = cms.bool(False),
    DoEtaWidthCut = cms.bool(False),
    DoHcalRecHitIsolationCut = cms.bool(True),
    DoEcalRecHitIsolationCut = cms.bool(True),
    DoR9Cut = cms.bool(True),                               

    #configuration
    isolationtrackThreshold = cms.double(0.0),
    TrackConeOuterRadius = cms.double(0.4),
    TrackConeInnerRadius = cms.double(0.04),
    EcalRecHitInnerRadius = cms.double(0.0),
    EcalRecHitOuterRadius = cms.double(0.4),
    EcalRecHitEtaSlice = cms.double(0.0),
    EcalRecThresh = cms.double(0.0),
    HcalRecHitInnerRadius = cms.double(0.1),
    HcalRecHitOuterRadius = cms.double(0.4),
    HcalRecHitEtaSlice = cms.double(0.),
    HcalRecHitThresh = cms.double(0.0),

    #LooseEM cuts EB
    LooseEMHollowTrkEB = cms.double(999.9),
    LooseEMSolidTrkEB  = cms.double(999.9),
    LooseEMSolidNTrkEB = cms.int32(999),
    LooseEMHollowNTrkEB = cms.int32(999),
    LooseEMEtaWidthEB = cms.double(999.9),
    LooseEMHadOverEMEB = cms.double(999.9),
    LooseEMEcalRecHitIsoEB = cms.double(20.0),
    LooseEMHcalRecHitIsoEB = cms.double(10.0),
    LooseEMR9CutEB = cms.double(0.0),

    #LoosePhoton cuts EB  
    LoosePhotonHollowTrkEB = cms.double(30.0),
    LoosePhotonSolidTrkEB  = cms.double(999.9),
    LoosePhotonSolidNTrkEB = cms.int32(999),
    LoosePhotonHollowNTrkEB = cms.int32(999),
    LoosePhotonEtaWidthEB = cms.double(999.9),
    LoosePhotonHadOverEMEB = cms.double(999.9),
    LoosePhotonEcalRecHitIsoEB = cms.double(20.0),
    LoosePhotonHcalRecHitIsoEB = cms.double(10.0),
    LoosePhotonR9CutEB = cms.double(0.0),

    #TightPhoton cuts EB
    TightPhotonHollowTrkEB = cms.double(30.0),
    TightPhotonSolidTrkEB  = cms.double(999.9),
    TightPhotonSolidNTrkEB = cms.int32(999),
    TightPhotonHollowNTrkEB = cms.int32(999),
    TightPhotonEtaWidthEB = cms.double(999.9),
    TightPhotonHadOverEMEB = cms.double(999.9),
    TightPhotonEcalRecHitIsoEB = cms.double(20.0),
    TightPhotonHcalRecHitIsoEB = cms.double(10.0),
    TightPhotonR9CutEB = cms.double(0.8),

    #LooseEM cuts EB
    LooseEMHollowTrkEE = cms.double(999.9),
    LooseEMSolidTrkEE  = cms.double(999.9),
    LooseEMSolidNTrkEE = cms.int32(999),
    LooseEMHollowNTrkEE = cms.int32(999),
    LooseEMEtaWidthEE = cms.double(999.9),
    LooseEMHadOverEMEE = cms.double(999.9),
    LooseEMEcalRecHitIsoEE = cms.double(20.0),
    LooseEMHcalRecHitIsoEE = cms.double(10.0),
    LooseEMR9CutEE = cms.double(0.0),

    #LoosePhoton cuts EB  
    LoosePhotonHollowTrkEE = cms.double(30.0),
    LoosePhotonSolidTrkEE  = cms.double(999.9),
    LoosePhotonSolidNTrkEE = cms.int32(999),
    LoosePhotonHollowNTrkEE = cms.int32(999),
    LoosePhotonEtaWidthEE = cms.double(999.9),
    LoosePhotonHadOverEMEE = cms.double(999.9),
    LoosePhotonEcalRecHitIsoEE = cms.double(20.0),
    LoosePhotonHcalRecHitIsoEE = cms.double(10.0),
    LoosePhotonR9CutEE = cms.double(0.0),

    #TightPhoton cuts EB
    TightPhotonHollowTrkEE = cms.double(30.0),
    TightPhotonSolidTrkEE  = cms.double(999.9),
    TightPhotonSolidNTrkEE = cms.int32(999),
    TightPhotonHollowNTrkEE = cms.int32(999),
    TightPhotonEtaWidthEE = cms.double(999.9),
    TightPhotonHadOverEMEE = cms.double(999.9),
    TightPhotonEcalRecHitIsoEE = cms.double(20.0),
    TightPhotonHcalRecHitIsoEE = cms.double(10.0),
    TightPhotonR9CutEE = cms.double(0.8)

    )

