import FWCore.ParameterSet.Config as cms

process = cms.Process("patCrossCleaning")

# Message logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATLayer0Summary')
process.MessageLogger.categories.append('PatCrossCleaner')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default = cms.untracked.PSet( limit = cms.untracked.int32(0) ),
    PATLayer0Summary = cms.untracked.PSet( limit = cms.untracked.int32(-1) ),
    PatCrossCleaner = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)

## Uncomment this to enable debug output for CrossCleaner
#process.MessageLogger.debugModules = ['*'] 
#process.MessageLogger.cout = cms.untracked.PSet(
#    threshold = cms.untracked.string('DEBUG'),
#    default = cms.untracked.PSet( limit = cms.untracked.int32(0) ),
#    PatCrossCleaner = cms.untracked.PSet ( limit = cms.untracked.int32(-1) )
#)

# Turn on summaries at the end of job
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# Number of events to analyzer
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) ) 

# Input file
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
				'/store/relval/CMSSW_2_2_0_pre1/RelValZTT/GEN-SIM-RECO/STARTUP_V7_v1/0000/08E9EC80-C7AE-DD11-9F8A-000423D99CEE.root'
        )
)

# Output file
process.out = cms.OutputModule("PoolOutputModule",
    fileName       = cms.untracked.string('PATcrosscleaned.root'),
    verbose        = cms.untracked.bool(False),
    outputCommands = cms.untracked.vstring('keep *_*_*_USER')
)

# CaloTowerConstituentsMap needed for Electron/Photon-Jet cleaning
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
)

# Cross-cleaner setup
process.load("SusyAnalysis.PatCrossCleaner.patCrossCleaner_cfi")
# Switch on/off some components
process.patcrosscleaner.doMuonJetCC        = True
process.patcrosscleaner.doElectronJetCC    = True
process.patcrosscleaner.doPhotonJetCC      = False
process.patcrosscleaner.doElectronPhotonCC = False
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
process.patcrosscleaner.MuonJetCrossCleaning.MuonID = 'TM2DCompatibilityTight'
process.patcrosscleaner.MuonJetCrossCleaning.modifyJetEnergy = True

# Include PAT Layer 0 & 1 if not running on pattified data
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V4::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")
## Necessary fixes to run 2.2.X on 2.1.X data
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
run22XonSummer08AODSIM(process)


# Finally: set the paths
process.path = cms.Path(
  process.patLayer0 +
  process.patLayer1 +
  process.patcrosscleaner
)
process.outpath = cms.EndPath(process.out)
# save PAT Layer 1 output
process.load("PhysicsTools.PatAlgos.patLayer1_EventContent_cff")
process.out.outputCommands.extend(process.patLayer1EventContent.outputCommands)

