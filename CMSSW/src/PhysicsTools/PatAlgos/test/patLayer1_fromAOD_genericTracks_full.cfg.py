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


process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V4::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")

from PhysicsTools.PatAlgos.tools.trackTools import *
makeTrackCandidates(process, 
        label='TrackCands',                   # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
        tracks=cms.InputTag('generalTracks'), # input track collection
        particleType="pi+",                   # particle type (for assigning a mass)
        preselection='pt > 10',               # preselection cut on candidates. Only methods of 'reco::Candidate' are available
        cleaning=True,                        # Run the PATGenericParticleCleaner ('allLayer0TrackCands')
        selection='pt > 10',                  # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
        isolation=['tracker','caloTowers'],   # Isolations to use (currently only 'tracker' and 'caloTowers' are valid options)
                                              # You can also set it to 'None' to avoid computing isolatio alltogether
        mcAs='allLayer0Muons',                # Replicate MC match as the one used for Muons
        triggerAs=['allLayer0Muons']          # Replicate trigger match as all the ones used for Muons
        );                                    #  you can specify more than one collection for this


process.p = cms.Path(
        process.patLayer0 *
        process.patLayer1 
)

# Output module configuration
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('PATLayer1_Output.fromAOD_genericTracks_full.root'),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('drop *')
)
process.outpath = cms.EndPath(process.out)
# save PAT Layer 1 output
process.load("PhysicsTools.PatAlgos.patLayer1_EventContent_cff")
process.out.outputCommands.extend(process.patLayer1EventContent.outputCommands)
process.out.outputCommands.append('keep *_selectedLayer1TrackCands_*_*')

#### Dump the python config
#
# f = open("patLayer1_fromAOD_genericTracks_full.dump.py", "w")
# f.write(process.dumpPython())
# f.close()
#
 
