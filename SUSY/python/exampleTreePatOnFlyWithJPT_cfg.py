# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(100)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('Configuration.StandardSequences.Services_cff')
process.add_( cms.Service( "TFileService",
    fileName = cms.string( 'patTree.root' ),
    closeFileFast = cms.untracked.bool(True) ) 
)

#-- Input Source --------------------------------------------------------------
### process.load("SUSYBSMAnalysis.SusyCAF.dataSamples.MinBias_Summer09_STARTUP3X_V8D_900GeV_v1_GEN_SIM_RECO_cff")
source = cms.Source('PoolSource',
fileNames = cms.untracked.vstring('/castor/cern.ch/user/n/nmohr/V7production/QCDDiJet_Pt380to470_MC_31X_V9_ReReco332.root')
)

process.maxEvents.input = 100

# Due to problem in production of LM samples: same event number appears multiple times
#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')


#-- Calibration tag -----------------------------------------------------------
# Should match input file's tag
# First Collision Data
#process.GlobalTag.globaltag = 'GR09_P_V6::All'
# MinBias MC STARTUP
process.GlobalTag.globaltag = 'STARTUP3X_V8D::All'

##### JPT Corrections on-the-fly #####
process.load("bainbrid.SUSY.JptOnTheFly_cff")

addJetCollection(process,
                 cms.InputTag('ak5sJPT'),
                 'AK5sJPT',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = None,
                 doType1MET   = False,
                 doL1Cleaning = True,
                 doL1Counters = True,
                 doJetID      = False,
                 genJetCollection = cms.InputTag("ak5GenJets")
                 )

addJetCollection(process,
                 cms.InputTag('ak5vJPT'),
                 'AK5vJPT',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = None,
                 doType1MET   = False,
                 doL1Cleaning = True,
                 doL1Counters = True,
                 doJetID      = False,
                 genJetCollection = cms.InputTag("ak5GenJets")
                 )

############################# START SUSYPAT specifics ####################################
from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands, removeMCDependence
#removeMCDependence(process)
addDefaultSUSYPAT(process, True,'HLT','900GeV')
SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process )
############################## END SUSYPAT specifics ####################################

#-- Output module configuration -----------------------------------------------
process.out.fileName = 'SUSYPAT.root'       # <-- CHANGE THIS TO SUIT YOUR NEEDS

# Custom settings
process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files???)
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
process.out.outputCommands = cms.untracked.vstring('drop *', *SUSY_pattuple_outputCommands )

# SusyCAF_nTuple
process.load('SUSYBSMAnalysis.SusyCAF.theBigNtuplePatOnFlyWithJPT_cfi')

#-- Execution path ------------------------------------------------------------
# Full path
process.p = cms.Path( process.JPT * process.seqSUSYDefaultSequence * process.theBigNtuplePat)
# if hasattr(process,"out"): # remove outpath 
#     del process.out
#     del process.outpath

