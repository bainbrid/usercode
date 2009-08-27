import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.Timing = cms.Service("Timing")

process.Tracer = cms.Service(
    "Tracer",
    sourceSeed = cms.untracked.string("$$")
    )

process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_3_2_5/RelValQCD_Pt_80_120/GEN-SIM-RECO/MC_31X_V5-v2/0011/96C592E4-CF90-DE11-9D21-000423D992A4.root',
    '/store/relval/CMSSW_3_2_5/RelValQCD_Pt_80_120/GEN-SIM-RECO/MC_31X_V5-v2/0011/427A4C5E-D290-DE11-94C2-000423D6C8E6.root',
    '/store/relval/CMSSW_3_2_5/RelValQCD_Pt_80_120/GEN-SIM-RECO/MC_31X_V5-v2/0011/0A75D05C-D390-DE11-A20A-000423D98804.root',
    ),
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('MC_31X_V5::All')

process.load("JetMETCorrections.Configuration.ZSPJetCorrections219_cff")
process.load("JetMETCorrections.Configuration.JetPlusTrackCorrections_cff")
process.load("bainbrid.Test.JPTCorrections_cff")
process.JPTCorrectionIC5.Verbose = cms.bool(False)

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
    process.ZSPJetCorrections 
    * process.JetPlusTrackCorrections
    * process.JPTCorrections 
    )

process.o = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    outputCommands = cms.untracked.vstring(
    'drop *',
    'keep recoCaloJets_*_*_*'
    )
    )

process.e = cms.EndPath( process.o )
 
process.ProfilerService = cms.Service(
    "ProfilerService",
    firstEvent = cms.untracked.int32(2),
    lastEvent = cms.untracked.int32(12),
    paths = cms.untracked.vstring("p")
    )

process.MessageLogger = cms.Service(
    "MessageLogger",
    
    debug = cms.untracked.PSet(
    threshold = cms.untracked.string('DEBUG'),
    limit = cms.untracked.uint32(100000),
    noLineBreaks = cms.untracked.bool(False),
    ),

    info = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'),
    limit = cms.untracked.uint32(100000),
    noLineBreaks = cms.untracked.bool(False),
    ),

    warning = cms.untracked.PSet(
    threshold = cms.untracked.string('WARNING'),
    limit = cms.untracked.uint32(100000),
    noLineBreaks = cms.untracked.bool(False),
    ),

    error = cms.untracked.PSet(
    threshold = cms.untracked.string('ERROR'),
    limit = cms.untracked.uint32(100000),
    noLineBreaks = cms.untracked.bool(False),
    ),
    
    cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('ERROR'),
    limit = cms.untracked.uint32(100000),
    noLineBreaks = cms.untracked.bool(False),
    ),

    destinations = cms.untracked.vstring(
    'debug', 
    'info', 
    'warning', 
    'error',
    'cerr'
    ),
    
    #@@ comment to suppress debug statements!
    debugModules = cms.untracked.vstring('JPTCorrectorIC5'),
    
    # allows to suppress output from specific modules 
    suppressDebug = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    
    )


