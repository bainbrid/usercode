import FWCore.ParameterSet.Config as cms

# 1: GMSB,22X,Summer08Redigi 
# 2: PhotonJets,21X,Summer08Redigi 
# 3: PhotonJets,22X,Summer08Redigi(Fall08) 
Vers = str("1")

Location = str("/home/bainbrid/susy/MISC/")

if Vers == str("1") :
    Source  = str("GMSB") # Options: "GMSB", "PhotonJets"
    Release = str("22X") # Options: "22X", "21X"
    Corrs   = str("Summer08Redigi") # Options: "Summer08Redigi", "Summer08", "Winter09"
elif Vers == str("2") :
    Source  = str("PhotonJets") # Options: "GMSB", "PhotonJets"
    Release = str("21X") # Options: "22X", "21X"
    Corrs   = str("Summer08Redigi") # Options: "Summer08Redigi", "Summer08", "Winter09"
elif Vers == str("3") :
    Source  = str("PhotonJets") # Options: "GMSB", "PhotonJets"
    Release = str("22X") # Options: "22X", "21X"
    Corrs   = str("Summer08Redigi") # Options: "Summer08Redigi", "Summer08", "Winter09"
else :
    print "UNKNOWN VERSION!"

Name = str(Source + "_" + Release + "_" + Corrs)

process = cms.Process("ANALYSIS")

process.Timing =cms.Service("Timing")
process.Tracer = cms.Service(
    "Tracer",
    sourceSeed = cms.untracked.string("$$")
    )

process.load("DQM.SiStripCommon.MessageLogger_cfi")
process.MessageLogger.debugModules = cms.untracked.vstring('')
process.MessageLogger.destinations = cms.untracked.vstring(
    "cerr", "susyPhotonID_" + Name, "info", "warning", "error"
    )

inputFiles = cms.untracked.vstring("file:" + Location + "susyPatLayer1_" + Name + ".root")
process.source = cms.Source(
    "PoolSource",
    fileNames = inputFiles
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    )

photonPdgIds = cms.vint32(22)
process.matchedPhotons = cms.EDFilter(
    "MCMatcher",
    src = cms.InputTag("allLayer0Photons"),
    matched = cms.InputTag("genParticles"),
    mcPdgId = photonPdgIds,
    mcStatus = cms.vint32(1),
    checkCharge = cms.bool(False),
    maxDeltaR = cms.double(0.4),
    maxDPtRel = cms.double(1.0),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    )

electronPdgIds = cms.vint32(11)
process.matchedElectrons = cms.EDFilter(
    "MCMatcher",
    src = cms.InputTag("allLayer0Photons"),
    matched = cms.InputTag("genParticles"),
    mcPdgId = electronPdgIds,
    mcStatus = cms.vint32(1),
    checkCharge = cms.bool(False),
    maxDeltaR = cms.double(0.4),
    maxDPtRel = cms.double(1.0),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    )

particlePdgIds = cms.vint32()
process.matchedParticles = cms.EDFilter(
    "MCMatcher",
    src = cms.InputTag("allLayer0Photons"),
    matched = cms.InputTag("genParticles"),
    mcPdgId = particlePdgIds,
    mcStatus = cms.vint32(1),
    checkCharge = cms.bool(False),
    maxDeltaR = cms.double(0.4),
    maxDPtRel = cms.double(5.0),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    )

# Rebuild pat::Photons with rebuilt MC matching
from PhysicsTools.PatAlgos.producersLayer1.photonProducer_cfi import allLayer1Photons
process.matchedLayer1Photons = allLayer1Photons.clone()
process.matchedLayer1Photons.photonSource = cms.InputTag("allLayer0Photons")
process.matchedLayer1Photons.genParticleMatch = cms.VInputTag(
    cms.InputTag("matchedPhotons"),
    cms.InputTag("matchedElectrons"),
    cms.InputTag("matchedParticles"),
    )

# Take rebuilt pat::Photons as input
process.load("bainbrid.Test.patPhotonIDProducer_cff")
process.patPhotonIDProducer.PatPhotons = cms.InputTag("matchedLayer1Photons")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
)

process.load("SusyAnalysis.PatCrossCleaner.patCrossCleaner_cfi")
process.patcrosscleaner.doMuonJetCC         = True
process.patcrosscleaner.doElectronJetCC     = True
process.patcrosscleaner.doPhotonJetCC       = True
process.patcrosscleaner.doElectronPhotonCC  = True
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
process.patcrosscleaner.MuonJetCrossCleaning.deltaR_min   = 0.2
process.patcrosscleaner.MuonJetCrossCleaning.caloIso_max  = 10.0
process.patcrosscleaner.MuonJetCrossCleaning.trackIso_max = 10.0
process.patcrosscleaner.MuonJetCrossCleaning.MuonID = 'AllGlobalMuons'

process.load("JetMETCorrections.Configuration.L7PartonCorrections_cff")
if Corrs == str("Summer08Redigi") :
    process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08Redigi_cff")
elif Corrs == str("Summer08") :
    process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08_cff")
elif Corrs == str("Winter09") :
    process.load("JetMETCorrections.Configuration.L2L3Corrections_Winter09_cff")
else :
    print "UNKNOWN CORRECTIONS!"

process.load("bainbrid.Test.simplePhotonIDAnalysis_cfi")
process.simplePhotonIDAnalysis.Photons        = cms.InputTag("patcrosscleaner:ccPhotons")
process.simplePhotonIDAnalysis.OtherPhotons   = cms.InputTag("")
process.simplePhotonIDAnalysis.PhotonPdgIds   = photonPdgIds
process.simplePhotonIDAnalysis.ElectronPdgIds = electronPdgIds 
process.simplePhotonIDAnalysis.ParticlePdgIds = particlePdgIds

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("file:" + Location + "susyPhotonID_" + Name + ".root")
    )

process.p = cms.Path(
    process.matchedPhotons *
    process.matchedElectrons *
    process.matchedParticles *
    process.matchedLayer1Photons *
    process.patPhotonIDProducer *
    process.patcrosscleaner *
    process.simplePhotonIDAnalysis
    )

