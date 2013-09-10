import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLE")

process.Timing =cms.Service("Timing")

process.Tracer = cms.Service(
    "Tracer",
    sourceSeed = cms.untracked.string("$$")
    )

process.load("DQM.SiStripCommon.MessageLogger_cfi")
process.MessageLogger.debugModules = cms.untracked.vstring('')

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('IDEAL_V12::All')

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring('file:/data2/bainbrid/data/QCD1000toInf-madgraph_229_SUSYPAT_V5_v1_6e6ee89347d8aaa6f7b14dff1a0707fb_patLayer1_1.root')
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
)

process.load("SusyAnalysis.PatCrossCleaner.patCrossCleaner_cfi")
process.patcrosscleaner.patMuons      = cms.InputTag("allLayer1Muons")
process.patcrosscleaner.patElectrons  = cms.InputTag("allLayer1Electrons")
process.patcrosscleaner.patPhotons    = cms.InputTag("allLayer1Photons")
process.patcrosscleaner.patTaus       = cms.InputTag("allLayer1Taus")
process.patcrosscleaner.patJets       = cms.InputTag("allLayer1JetsIC5")
process.patcrosscleaner.patMets       = cms.InputTag("allLayer1METsIC5")

process.patcrosscleaner.doMuonJetCC        = True
process.patcrosscleaner.doElectronJetCC    = True
process.patcrosscleaner.doPhotonJetCC      = True
process.patcrosscleaner.doElectronPhotonCC = True

process.patcrosscleaner.L1JetCorrector      = 'none'
process.patcrosscleaner.L2JetCorrector      = 'L2RelativeJetCorrectorIC5Calo'
process.patcrosscleaner.L3JetCorrector      = 'L3AbsoluteJetCorrectorIC5Calo'
process.patcrosscleaner.L4JetCorrector      = 'none'
process.patcrosscleaner.L5udsJetCorrector   = 'none'
process.patcrosscleaner.L5gluonJetCorrector = 'none'
process.patcrosscleaner.L5cJetCorrector     = 'none'
process.patcrosscleaner.L5bJetCorrector     = 'none'
process.patcrosscleaner.L6JetCorrector      = 'none'
process.patcrosscleaner.L7udsJetCorrector   = 'L7PartonJetCorrectorIC5qJet'
process.patcrosscleaner.L7gluonJetCorrector = 'L7PartonJetCorrectorIC5gJet'
process.patcrosscleaner.L7cJetCorrector     = 'L7PartonJetCorrectorIC5cJet'
process.patcrosscleaner.L7bJetCorrector     = 'L7PartonJetCorrectorIC5bJet'

process.patcrosscleaner.MuonJetCrossCleaning.deltaR_min   = 0.2
process.patcrosscleaner.MuonJetCrossCleaning.caloIso_max  = 10.0
process.patcrosscleaner.MuonJetCrossCleaning.trackIso_max = 10.0
process.patcrosscleaner.MuonJetCrossCleaning.MuonID = 'AllGlobalMuons'

process.patcrosscleaner.ElectronJetCrossCleaning.SusyAnalyzerCleaning = True
process.patcrosscleaner.ElectronJetCrossCleaning.deltaR_min        = 0.5
process.patcrosscleaner.ElectronJetCrossCleaning.SharedEtoJetE     = 0.7
process.patcrosscleaner.ElectronJetCrossCleaning.SharedEForNIsoEle = -1.
process.patcrosscleaner.ElectronJetCrossCleaning.IsolationKey  = 'TrackerIso'
process.patcrosscleaner.ElectronJetCrossCleaning.IsoValueCut   = 1.
process.patcrosscleaner.ElectronJetCrossCleaning.ElectronID   = 'eidRobustLoose'

process.patcrosscleaner.patPhotons = cms.InputTag("patPhotonIDProducer")
process.patcrosscleaner.PhotonJetCrossCleaning.deltaR_min   = 0.4
process.patcrosscleaner.PhotonJetCrossCleaning.IsoValueCut  = 1000000.
process.patcrosscleaner.PhotonJetCrossCleaning.IsolationKey = 'CaloIso'
process.patcrosscleaner.PhotonJetCrossCleaning.SharedEtoJetE = 0.7
process.patcrosscleaner.PhotonJetCrossCleaning.PhotonID = 'TightPhoton'

process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")

process.load("PF.Susy.patFromPF2PAT_cff")
process.load("PF.Susy.PF2PAT_cff")

process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08Redigi_cff")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("bkgd.root")
    )

process.load("bainbrid.Test.patPhotonIDProducer_cff")
process.patPhotonIDProducer.Photons = cms.InputTag("allLayer1Photons")

process.test = cms.EDAnalyzer(
    "TestCombination",
    # Misc
    MaximumObjects = cms.int32(30),
    TestObjects    = cms.int32(-1),
    # Collections
    Photons        = cms.InputTag("patcrosscleaner:ccPhotons"),
    Jets           = cms.InputTag("patcrosscleaner:ccJets"),
    Muons          = cms.InputTag("patcrosscleaner:ccMuons"),
    Electrons      = cms.InputTag("patcrosscleaner:ccElectrons"),
    # Object selection
    PhotonEt       = cms.double(50.),
    PhotonEta      = cms.double(2.4),
    JetEt          = cms.double(50.),
    JetEta         = cms.double(2.4),
    JetEMfraction  = cms.double(0.9),
    MuonPt         = cms.double(10.),
    MuonEta        = cms.double(2.4),
    MuonTrkIso     = cms.double(10.0),
    ElectronPt     = cms.double(10.),
    ElectronEta    = cms.double(2.4),
    ElectronTrkIso = cms.double(1.0),
    # Event selection
    TotalEt        = cms.double(350.),
    MinJets        = cms.untracked.int32(0),
    MaxJets        = cms.untracked.int32(30),
    )

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
    process.dump *
    process.patPhotonIDProducer *
    process.patcrosscleaner *
    process.test
    )

