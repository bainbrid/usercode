import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLE")

process.Timing =cms.Service("Timing")

process.MessageLogger = cms.Service(
    "MessageLogger",
    ntuple = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO')
    ),
    destinations = cms.untracked.vstring('ntuple')
    )

# Include PAT Layer 0 & 1 if not running on pattified data
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('IDEAL_V11::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

# CaloTowerConstituentsMap needed for Electron/Photon-Jet cleaning
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.CaloTowerConstituentsMapBuilder = cms.ESProducer(
    "CaloTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
    )

# Cross-cleaner setup
process.load("SusyAnalysis.PatCrossCleaner.patCrossCleaner_cfi")
# Switch on/off some components
process.patcrosscleaner.doMuonJetCC        = True
#process.patcrosscleaner.doElectronJetCC    = True
process.patcrosscleaner.doElectronJetCC    = True
process.patcrosscleaner.doPhotonJetCC      = True
process.patcrosscleaner.doElectronPhotonCC = True
# Change the jet energy corrections
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

# Parameters for electron-jet cross-cleaning
process.patcrosscleaner.ElectronJetCrossCleaning.SusyAnalyzerCleaning = True
process.patcrosscleaner.ElectronJetCrossCleaning.deltaR_min        = 0.5
process.patcrosscleaner.ElectronJetCrossCleaning.SharedEtoJetE     = 0.7
process.patcrosscleaner.ElectronJetCrossCleaning.SharedEForNIsoEle = -1.
process.patcrosscleaner.ElectronJetCrossCleaning.IsolationKey  = 'TrackerIso'
process.patcrosscleaner.ElectronJetCrossCleaning.IsoValueCut   = 1.
process.patcrosscleaner.ElectronJetCrossCleaning.ElectronID   = 'eidRobustLoose'
# Parameters for photon-jet cross-cleaning
process.patcrosscleaner.PhotonJetCrossCleaning.deltaR_min   = 0.5
process.patcrosscleaner.PhotonJetCrossCleaning.IsoValueCut  = 0.3
process.patcrosscleaner.PhotonJetCrossCleaning.IsolationKey = 'CaloIso'
process.patcrosscleaner.PhotonJetCrossCleaning.PhotonID = 'TightPhoton'
# Parameters for muon-jet cross-cleaning
process.patcrosscleaner.MuonJetCrossCleaning.deltaR_min   = 0.2
process.patcrosscleaner.MuonJetCrossCleaning.caloIso_max  = 10.0
process.patcrosscleaner.MuonJetCrossCleaning.trackIso_max = 10.0
process.patcrosscleaner.MuonJetCrossCleaning.MuonID = 'TMLastStationTight'

# clone for JPT corrections
process.load("JetMETCorrections.JetPlusTrack.L1L3Corrections_JPT_cff")
process.patcrosscleanerJPT = process.patcrosscleaner.clone()
process.patcrosscleanerJPT.L1JetCorrector      = 'L1OffsetJetCorrectorZSP'
process.patcrosscleanerJPT.L2JetCorrector      = 'none'
process.patcrosscleanerJPT.L3JetCorrector      = 'L3AbsoluteJetCorrectorJPT'
process.patcrosscleanerJPT.L7udsJetCorrector   = 'none'
process.patcrosscleanerJPT.L7gluonJetCorrector = 'none'
process.patcrosscleanerJPT.L7cJetCorrector     = 'none'
process.patcrosscleanerJPT.L7bJetCorrector     = 'none'
process.patcrosscleanerJPT.MuonJetCrossCleaning.modifyJetEnergy = False

process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")
process.load("PF.Susy.patFromPF2PAT_cff")
process.load("PF.Susy.PF2PAT_cff")
process.load("PF.Susy.patLayer1_EventContent_cff")

#process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08_cff")
#process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08Redigi_cff")
#process.load("JetMETCorrections.Configuration.L2L3Corrections_Winter09_cff")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string('ntuple.root')
    )

process.Tracer = cms.Service(
    "Tracer",
    sourceSeed = cms.untracked.string("$$")
    )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    'file:pat.root'
    )
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

#process.genParticles.abortOnUnknownPDGCode = False

process.load("PF.Susy.hemisphere_cfi")
process.load("PF.Susy.pfhemisphere_cfi")
process.load("SusyAnalysis.AnalysisSkeleton.dijet_cfi")
process.load("SusyAnalysis.AnalysisSkeleton.pfdijet_cfi")

process.p = cms.Path(
    process.patcrosscleaner *
    process.patcrosscleanerJPT * 
    process.selectedLayer2Hemispheres *
    process.dijet *
    process.pfemispheres *
    process.pfdijet
    )
