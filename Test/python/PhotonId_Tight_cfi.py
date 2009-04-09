import FWCore.ParameterSet.Config as cms

# Based on: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/RecoEgamma/PhotonIdentification/python/photonId_cfi.py?revision=1.7.2.1&view=markup

from RecoEgamma.PhotonIdentification.photonId_cfi import *

RegeneratedPhotonID = PhotonIDProd.clone(
    
    # General
    
    photonProducer = cms.string('photons'),
    photonLabel = cms.string(''),
    
    photonIDAssociationLabel = cms.string('RegeneratedPhotonAssociatedID'),
    photonIDLabel = cms.string('PhotonIDCutBasedProducer'),

    barrelEcalRecHitProducer = cms.string('ecalRecHit'),
    barrelEcalRecHitCollection = cms.string('EcalRecHitsEB'),
    endcapEcalRecHitProducer = cms.string('ecalRecHit'),
    endcapEcalRecHitCollection = cms.string('EcalRecHitsEE'),

    HcalRecHitProducer = cms.string('hbhereco'),
    HcalRecHitCollection = cms.string(''),
    
    GsfRecoCollection = cms.InputTag("pixelMatchGsfElectrons"),
    modulePhiBoundary =   cms.double( 0.0087 ),
    moduleEtaBoundary = cms.vdouble( 0.0, 0.05, 0.4, 0.5, 0.75, 0.85, 1.1, 1.2, 1.4, 1.6 ),
    
    trackProducer = cms.InputTag("generalTracks"),
    
    # Cut-based analysis and switches
    doCutBased                    = cms.bool(True),
    DoHollowConeTrackIsolationCut = cms.bool(True),
    DoEcalRecHitIsolationCut      = cms.bool(True),
    DoHcalRecHitIsolationCut      = cms.bool(True),
    DoR9Cut                       = cms.bool(True), 
    RequireNotElectron            = cms.bool(False),                    
    RequireFiducial               = cms.bool(False),
    DoSolidConeTrackIsolationCut  = cms.bool(False),
    DoHollowConeNTrkCut           = cms.bool(False),
    DoSolidConeNTrkCut            = cms.bool(False),
    DoHadOverEMCut                = cms.bool(False),
    DoEtaWidthCut                 = cms.bool(False),
    
    # "Jurassic" footprint for ECAL, "hollow cone" for tracks and HCAL
    
    TrackConeOuterRadius  = cms.double(0.4),
    TrackConeInnerRadius  = cms.double(0.04),
    EcalRecHitInnerRadius = cms.double(0.06),
    EcalRecHitOuterRadius = cms.double(0.4),
    EcalRecHitEtaSlice    = cms.double(0.04),
    HcalRecHitInnerRadius = cms.double(0.1),
    HcalRecHitOuterRadius = cms.double(0.4),

    # Loose EM (EB and EE)

    LooseEMEcalRecHitIsoEB = cms.double(10.0),
    LooseEMHcalRecHitIsoEB = cms.double(5.0),

    LooseEMEcalRecHitIsoEE = cms.double(10.0),
    LooseEMHcalRecHitIsoEE = cms.double(10.0),

    # Loose Photon (EB and EE)
    
    LoosePhotonEcalRecHitIsoEB = cms.double(10.0),
    LoosePhotonHcalRecHitIsoEB = cms.double(5.0),
    LoosePhotonHollowTrkEB     = cms.double(30.0),
    
    LoosePhotonEcalRecHitIsoEE = cms.double(10.0),
    LoosePhotonHcalRecHitIsoEE = cms.double(10.0),
    LoosePhotonHollowTrkEE     = cms.double(30.0),
    
    # Tight Photon (EB and EE)

    TightPhotonEcalRecHitIsoEB = cms.double(10.0),
    TightPhotonHcalRecHitIsoEB = cms.double(5.0),
    TightPhotonHollowTrkEB     = cms.double(30.0),
    TightPhotonR9CutEB         = cms.double(0.8),

    TightPhotonEcalRecHitIsoEE = cms.double(10.0),
    TightPhotonHcalRecHitIsoEE = cms.double(10.0),
    TightPhotonHollowTrkEE     = cms.double(30.0),
    TightPhotonR9CutEE         = cms.double(0.8),

    # Other miscellaneous

    )
