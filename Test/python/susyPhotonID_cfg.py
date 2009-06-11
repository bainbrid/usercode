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
    maxDeltaR = cms.double(0.2),
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
    maxDeltaR = cms.double(0.2),
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
    maxDeltaR = cms.double(0.2),
    maxDPtRel = cms.double(1.0),
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

process.load("bainbrid.Test.simplePhotonIDAnalysis_cfi")
process.simplePhotonIDAnalysis.PhotonPdgIds = photonPdgIds
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
    process.simplePhotonIDAnalysis
    )

