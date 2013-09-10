import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATLayer0Summary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    PATLayer0Summary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring('file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTBar-2_2_X_2008-11-03-STARTUP_V7-AODSIM.100.root')
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")

## Load additional RECO config
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V4::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

from PhysicsTools.PatAlgos.tools.jetTools import *

## ==== Example with CaloJets
switchJetCollection(process, 
        'sisCone5CaloJets',    # Jet collection; must be already in the event when patLayer0 sequence is executed
        layers=[0,1],          # If you're not runnint patLayer1, set 'layers=[0]' 
        runCleaner="CaloJet",  # =None if not to clean
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
        jetCorrLabel=('SC5','Calo'), # example jet correction name; set to None for no JEC
        doType1MET=True)       # recompute Type1 MET using these jets

## ==== FOR BASIC JETS
### make some basic jets for testing
# from RecoJets.JetProducers.BasicJetIcone5_cfi import iterativeCone5BasicJets
# process.iterativeCone5BasicJets = iterativeCone5BasicJets.clone(src = cms.InputTag("towerMaker"))
#
### configure PAT
# switchJetCollection(process, 
#        'iterativeCone5BasicJets', # Jet collection; must be already in the event patLayer0 sequence is executed
#        layers=[0,1],        
#        runCleaner="BasicJet", # =None for no cleaning
#        doJTA=True,
#        doBTagging=True,
#        jetCorrLabel=None,#=('S5','Calo'), # If you have JES corrections, you can apply them even to BasicJets
#        doType1MET=False)                  # Type1MET dows not work on BasicJets :-(

process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.p = cms.Path(
                ## process.iterativeCone5BasicJets +  ## Turn on this to run tests on BasicJets
                process.patLayer0  
                #+ process.content    # uncomment to get a dump 
                + process.patLayer1
            )


# Output module configuration
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('PATLayer1_Output.fromAOD_sisCone_full.root'),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('drop *')
)
process.outpath = cms.EndPath(process.out)
# save PAT Layer 1 output
process.load("PhysicsTools.PatAlgos.patLayer1_EventContent_cff")
process.out.outputCommands.extend(process.patLayer1EventContent.outputCommands)
