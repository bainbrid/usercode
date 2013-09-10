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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )


# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")

## Load additional RECO config
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V4::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

from PhysicsTools.PatAlgos.tools.jetTools import *

# Taking away BasicJets because RecoJets/JetProducers/python/BasicJetIcone5_cfi.py in 2.2.X is wrong
## FIXME ### make also some basic jets for testing
## FIXME from RecoJets.JetProducers.BasicJetIcone5_cfi import iterativeCone5BasicJets
## FIXME process.iterativeCone5BasicJets = iterativeCone5BasicJets.clone(src = cms.InputTag("towerMaker"))


addJetCollection(process,'sisCone5CaloJets','SC5',
                        runCleaner="CaloJet",doJTA=True,doBTagging=True,jetCorrLabel=('SC5','Calo'),doType1MET=True,doL1Counters=False)
addJetCollection(process,'sisCone7CaloJets','SC7',
                        runCleaner="CaloJet",doJTA=True,doBTagging=False,jetCorrLabel=None,doType1MET=True,doL1Counters=False)
addJetCollection(process,'kt4CaloJets','KT4',
                        runCleaner=None,doJTA=True,doBTagging=True,jetCorrLabel=('KT4','Calo'),doType1MET=True,doL1Counters=False)
addJetCollection(process,'kt6CaloJets','KT6',
                        runCleaner=None,doJTA=True,doBTagging=False,jetCorrLabel=None,doType1MET=True,doL1Counters=False)
## FIXME addJetCollection(process,'iterativeCone5BasicJets', 'BJ5',
## FIXME                        runCleaner="BasicJet",doJTA=True,doBTagging=True,jetCorrLabel=('MC5','Calo'),doType1MET=True,doL1Counters=False)
addJetCollection(process,'iterativeCone5PFJets', 'PFc',
                        runCleaner="PFJet",doJTA=True,doBTagging=True,jetCorrLabel=None,doType1MET=True,doL1Counters=False)
addJetCollection(process,'iterativeCone5PFJets', 'PFr',
                        runCleaner=None,doJTA=True,doBTagging=True,jetCorrLabel=None,doType1MET=True,doL1Counters=False)

process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.p = cms.Path(
                ## FIXME process.iterativeCone5BasicJets +  ## Turn on this to run tests on BasicJets
                process.patLayer0  
                + process.patLayer1
                + process.content    # uncomment to get a dump 
            )


# Output module configuration
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('PATLayer1_Output.fromAOD_jetSuite_full.root'),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('drop *')
)
process.outpath = cms.EndPath(process.out)
# save PAT Layer 1 output
process.load("PhysicsTools.PatAlgos.patLayer1_EventContent_cff")
process.out.outputCommands.extend(process.patLayer1EventContent.outputCommands)
process.out.outputCommands.extend(["keep *_selectedLayer1Jets*_*_*"])

#### Dump the python config
#
#f = open("patLayer1_fromAOD_jetSuite_full.dump.py", "w")
#f.write(process.dumpPython())
#f.close()
#
 
#### GraphViz dumps of sequences and modules, useful for debugging.
#### WARNING: it's not for the weak-hearted; the output plot is HUGE
#### needs a version of 'dot' that works with png graphics. 
#### in case, you can borrw mine with
####   export LD_LIBRARY_PATH=/afs/cern.ch/user/g/gpetrucc/scratch0/graphviz/lib:${LD_LIBRARY_PATH}
####   export PATH=/afs/cern.ch/user/g/gpetrucc/scratch0/graphviz/bin:${PATH}
#
# from PhysicsTools.PatAlgos.tools.circuitry import *
# plotSequences(   process.p, 'patSequences.png')
# plotModuleInputs(process.p, 'patModules.png'  ,printOuter=False,printLinkNames=True)  # no nodes for non-PAT modules
