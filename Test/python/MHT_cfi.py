import FWCore.ParameterSet.Config as cms

MHT = cms.EDAnalyzer(
    "MHT",
    
    Jets           = cms.InputTag("allLayer1Jets"),
    Muons          = cms.InputTag("allLayer1Muons"),
    Electrons      = cms.InputTag("allLayer1Electrons"),
    Photons        = cms.InputTag("allLayer1Photons"),
    
    CaloMET        = cms.InputTag("allLayer1METs"),
    GenMET         = cms.InputTag("genMet"),
    
    JetEt          = cms.double(50.),
    JetEta         = cms.double(2.4),
    JetEMfraction  = cms.double(0.9),
    
    MuonPt         = cms.double(10.),
    MuonEta        = cms.double(2.4),
    MuonTrkIso     = cms.double(10.0),
    
    ElectronPt     = cms.double(10.),
    ElectronEta    = cms.double(2.4),
    ElectronTrkIso = cms.double(1.0),
    
    PhotonEt       = cms.double(25.),
    PhotonEta      = cms.double(2.4),
    
    TotalEt        = cms.double(350.),
    
    MinObjects     = cms.int32(0),
    MinJets        = cms.int32(0),
    MinMuons       = cms.int32(0),
    MinElectrons   = cms.int32(0),
    MinPhotons     = cms.int32(0),
    
    )
