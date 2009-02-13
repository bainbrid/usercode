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

# FAMOS source
process.load("FastSimulation.Configuration.ttbar_cfi")
# some FAMOS setup (FIXME maybe this should go)
process.load("PhysicsTools.PatAlgos.famos.boostrapWithFamos_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")

## Load additional RECO config
# Magnetic field now needs to be in the high-level py
process.load("Configuration.StandardSequences.MagneticField_cff")

## Switch off trigger matching (just to make things easier)
from PhysicsTools.PatAlgos.tools.trigTools import switchLayer1Off
switchLayer1Off(process)

#process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.p = cms.Path(
                process.famosWithEverything        +    # run full FAMOS (ttbar events)
                #process.content                   +    # to get a dump of the output of FAMOS
                process.patLayer0_withoutTrigMatch +    # PAT Layer 0, no trigger matching
                process.patLayer1                       # PAT Layer 1
            )


# Output module configuration
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('PATLayer1_Output.fromScratch_fast.root'),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('drop *')
)
process.outpath = cms.EndPath(process.out)
# save PAT Layer 1 output
process.load("PhysicsTools.PatAlgos.patLayer1_EventContent_cff")
process.out.outputCommands.extend(process.patLayer1EventContent.outputCommands)

