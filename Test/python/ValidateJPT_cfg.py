import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.Timing = cms.Service("Timing")
process.Tracer = cms.Service("Tracer",sourceSeed = cms.untracked.string("$$"))

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP31X_V4::All')

# Input files: RelVal QCD 80-120 GeV, STARTUP conditions, 9000 events, from CMSSW_3_2_5 (replace with 33X when available!)
process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_3_2_5/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP31X_V4-v1/0011/2AA47C1E-828E-DE11-B3C5-001D09F34488.root',
    '/store/relval/CMSSW_3_2_5/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP31X_V4-v1/0010/DA6FF61D-3D8E-DE11-938A-003048D37538.root',
    '/store/relval/CMSSW_3_2_5/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP31X_V4-v1/0010/82CF8666-398E-DE11-8F3B-000423D94A20.root',
    '/store/relval/CMSSW_3_2_5/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP31X_V4-v1/0010/6255E85C-3F8E-DE11-B46A-000423D6B48C.root',
    '/store/relval/CMSSW_3_2_5/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP31X_V4-v1/0010/3453CD32-418E-DE11-87D2-003048D2C020.root',
    '/store/relval/CMSSW_3_2_5/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP31X_V4-v1/0010/2EC02533-3B8E-DE11-BE85-003048D37514.root',
    ),
    )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


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


