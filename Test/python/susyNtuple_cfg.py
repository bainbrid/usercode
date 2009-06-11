import FWCore.ParameterSet.Config as cms

# 1: GMSB,22X,Summer08Redigi 
# 2: PhotonJets,21X,Summer08Redigi 
Vers = str("1")

if Vers == str("1") :
    Source  = str("GMSB") # Options: "GMSB", "PhotonJets"
    Release = str("22X") # Options: "22X", "21X"
    Corrs   = str("Summer08Redigi") # Options: "Summer08Redigi", "Summer08", "Winter09"
elif Vers == str("2") :
    Source  = str("PhotonJets") # Options: "GMSB", "PhotonJets"
    Release = str("21X") # Options: "22X", "21X"
    Corrs   = str("Summer08Redigi") # Options: "Summer08Redigi", "Summer08", "Winter09"
else :
    print "UNKNOWN VERSION!"

Name = str(Source + "_" + Release + "_" + Corrs)

process = cms.Process("NTUPLE")

process.Timing =cms.Service("Timing")
process.Tracer = cms.Service(
    "Tracer",
    sourceSeed = cms.untracked.string("$$")
    )

process.load("DQM.SiStripCommon.MessageLogger_cfi")
process.MessageLogger.debugModules = cms.untracked.vstring('')
process.MessageLogger.destinations = cms.untracked.vstring(
    "cerr", "susyNtuple_" + Name, "info", "warning", "error"
    )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

inputFiles = cms.untracked.vstring("file:./" + "susyPatLayer1_" + Name + ".root")
process.source = cms.Source(
    "PoolSource",
    fileNames = inputFiles
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
)

process.load("SusyAnalysis.PatCrossCleaner.patCrossCleaner_cfi")
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

process.patcrosscleaner.ElectronJetCrossCleaning.SusyAnalyzerCleaning = True
process.patcrosscleaner.ElectronJetCrossCleaning.deltaR_min        = 0.5
process.patcrosscleaner.ElectronJetCrossCleaning.SharedEtoJetE     = 0.7
process.patcrosscleaner.ElectronJetCrossCleaning.SharedEForNIsoEle = -1.
process.patcrosscleaner.ElectronJetCrossCleaning.IsolationKey  = 'TrackerIso'
process.patcrosscleaner.ElectronJetCrossCleaning.IsoValueCut   = 1.
process.patcrosscleaner.ElectronJetCrossCleaning.ElectronID   = 'eidRobustLoose'

process.patcrosscleaner.patPhotons = cms.InputTag("patPhotonIDProducer")
process.patcrosscleaner.PhotonJetCrossCleaning.deltaR_min   = 0.4
process.patcrosscleaner.PhotonJetCrossCleaning.IsoValueCut  = 1000.
process.patcrosscleaner.PhotonJetCrossCleaning.IsolationKey = 'CaloIso'
process.patcrosscleaner.PhotonJetCrossCleaning.PhotonID = 'TightPhoton'

process.patcrosscleaner.MuonJetCrossCleaning.deltaR_min   = 0.2
process.patcrosscleaner.MuonJetCrossCleaning.caloIso_max  = 10.0
process.patcrosscleaner.MuonJetCrossCleaning.trackIso_max = 10.0
process.patcrosscleaner.MuonJetCrossCleaning.MuonID = 'AllGlobalMuons'

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

if Corrs == str("Summer08Redigi") :
    process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08Redigi_cff")
elif Corrs == str("Summer08") :
    process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08_cff")
elif Corrs == str("Winter09") :
    process.load("JetMETCorrections.Configuration.L2L3Corrections_Winter09_cff")
else :
    print "UNKNOWN CORRECTIONS!"

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)


process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("susyNtuple_" + Name + ".root")
    )

process.load("bainbrid.Test.patPhotonIDProducer_cff")

process.common = cms.Sequence(
    process.patPhotonIDProducer *
    process.patcrosscleaner *
    process.patcrosscleanerJPT 
    )

if Release == str("22X") :
    process.GlobalTag.globaltag = cms.string('IDEAL_V11::All')
    process.load("PF.Susy.hemisphere_cfi")
    process.load("bainbrid.Test.dijet_cfi")
    process.dijet.photTag = cms.InputTag("patPhotonIDProducer")
    process.load("PF.Susy.pfhemisphere_cfi")
    process.load("bainbrid.Test.pfdijet_cfi")
    process.pfdijet.photTag = cms.InputTag("patPhotonIDProducer")
    process.p = cms.Path(
        process.common *
        process.selectedLayer2Hemispheres *
        process.dijet *
        process.pfemispheres *
        process.pfdijet
        )
elif Release == str("21X") :
    process.p = cms.Path(
        process.common
        )
else :
    print "UNKNOWN RELEASE!"
