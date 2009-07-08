import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

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
process.GlobalTag.globaltag = cms.string('IDEAL_V11::All')

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring('file:/data2/bainbrid/data/gmsb/pat/GMSB_GM1g_1.root')
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
)

process.load("SusyAnalysis.PatCrossCleaner.patCrossCleaner_cfi")
process.patcrosscleaner.doMuonJetCC        = False
process.patcrosscleaner.doElectronJetCC    = False
process.patcrosscleaner.doElectronPhotonCC = False
process.patcrosscleaner.doPhotonJetCC      = True

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
    fileName = cms.string("output.root")
    )

process.load("bainbrid.Test.patPhotonIDProducer_cff")

process.jetMatchedCCPhotons = cms.EDFilter(
    "MCMatcherByPt",
    src = cms.InputTag("patcrosscleaner:ccJets"),
    matched = cms.InputTag("genParticles"),
    mcPdgId = cms.vint32(22),
    mcStatus = cms.vint32(1),
    maxDeltaR = cms.double(0.4),
    maxDPtRel = cms.double(1.0),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    )

process.jetMatchedCCElectrons = cms.EDFilter(
    "MCMatcherByPt",
    src = cms.InputTag("patcrosscleaner:ccJets"),
    matched = cms.InputTag("genParticles"),
    mcPdgId = cms.vint32(11),
    mcStatus = cms.vint32(1),
    maxDeltaR = cms.double(0.4),
    maxDPtRel = cms.double(1.0),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    )

process.jetMatchedCCPartons = cms.EDFilter(
    "MCMatcherByPt",
    src = cms.InputTag("patcrosscleaner:ccJets"),
    matched = cms.InputTag("genParticles"),
    mcPdgId = cms.vint32(1,2,3,4,5,21),
    mcStatus = cms.vint32(3),
    maxDeltaR = cms.double(0.4),
    maxDPtRel = cms.double(1.0),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    )

process.jetMatchedDroppedPhotons = cms.EDFilter(
    "MCMatcherByPt",
    src = cms.InputTag("patcrosscleaner:droppedJets"),
    matched = cms.InputTag("genParticles"),
    mcPdgId = cms.vint32(22),
    mcStatus = cms.vint32(1),
    maxDeltaR = cms.double(0.4),
    maxDPtRel = cms.double(1.0),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    )

process.jetMatchedDroppedElectrons = cms.EDFilter(
    "MCMatcherByPt",
    src = cms.InputTag("patcrosscleaner:droppedJets"),
    matched = cms.InputTag("genParticles"),
    mcPdgId = cms.vint32(11),
    mcStatus = cms.vint32(1),
    maxDeltaR = cms.double(0.4),
    maxDPtRel = cms.double(1.0),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    )

process.jetMatchedDroppedPartons = cms.EDFilter(
    "MCMatcherByPt",
    src = cms.InputTag("patcrosscleaner:droppedJets"),
    matched = cms.InputTag("genParticles"),
    mcPdgId = cms.vint32(1,2,3,4,5,21),
    mcStatus = cms.vint32(3),
    maxDeltaR = cms.double(0.4),
    maxDPtRel = cms.double(1.0),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    )

process.photonMatchedCCPhotons = cms.EDFilter(
    "MCMatcherByPt",
    src = cms.InputTag("patcrosscleaner:ccPhotons"),
    matched = cms.InputTag("genParticles"),
    mcPdgId = cms.vint32(22),
    mcStatus = cms.vint32(1),
    maxDeltaR = cms.double(0.4),
    maxDPtRel = cms.double(1.0),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    )

process.photonMatchedDroppedPhotons = cms.EDFilter(
    "MCMatcherByPt",
    src = cms.InputTag("patcrosscleaner:droppedPhotons"),
    matched = cms.InputTag("genParticles"),
    mcPdgId = cms.vint32(22),
    mcStatus = cms.vint32(1),
    maxDeltaR = cms.double(0.4),
    maxDPtRel = cms.double(1.0),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    )

process.test = cms.EDAnalyzer(
    "TestPatCrossCleaner",
    #
    RawPhotons       = cms.InputTag("patPhotonIDProducer"),
    RawJets          = cms.InputTag("selectedLayer1Jets"),
    RawElectrons     = cms.InputTag("selectedLayer1Electrons"),
    RawMuons         = cms.InputTag("selectedLayer1Muons"),
    #
    CCPhotons        = cms.InputTag("patcrosscleaner:ccPhotons"),
    CCJets           = cms.InputTag("patcrosscleaner:ccJets"),
    CCElectrons      = cms.InputTag("patcrosscleaner:ccElectrons"),
    CCMuons          = cms.InputTag("patcrosscleaner:ccMuons"),
    #
    DroppedPhotons   = cms.InputTag("patcrosscleaner:droppedPhotons"),
    DroppedJets      = cms.InputTag("patcrosscleaner:droppedJets"),
    DroppedElectrons = cms.InputTag("patcrosscleaner:droppedElectrons"),
    DroppedMuons     = cms.InputTag("patcrosscleaner:droppedMuons"),
    #
    JetMatchedCCPhotons   = cms.InputTag("jetMatchedCCPhotons"),
    JetMatchedCCElectrons = cms.InputTag("jetMatchedCCElectrons"),
    JetMatchedCCPartons   = cms.InputTag("jetMatchedCCPartons"),
    #
    JetMatchedDroppedPhotons   = cms.InputTag("jetMatchedDroppedPhotons"),
    JetMatchedDroppedElectrons = cms.InputTag("jetMatchedDroppedElectrons"),
    JetMatchedDroppedPartons   = cms.InputTag("jetMatchedDroppedPartons"),
    #
    PhotonMatchedCCPhotons      = cms.InputTag("photonMatchedCCPhotons"),
    PhotonMatchedDroppedPhotons = cms.InputTag("photonMatchedDroppedPhotons"),
    )

process.p = cms.Path(
    process.patPhotonIDProducer *
    process.patcrosscleaner *
    process.jetMatchedCCPhotons *
    process.jetMatchedCCElectrons *
    process.jetMatchedCCPartons *
    process.jetMatchedDroppedPhotons *
    process.jetMatchedDroppedElectrons *
    process.jetMatchedDroppedPartons *
    process.photonMatchedCCPhotons *
    process.photonMatchedDroppedPhotons *
    process.test
    )
