import FWCore.ParameterSet.Config as cms

process = cms.Process("ENERGYSCALE")

process.load("DQM.SiStripCommon.MessageLogger_cfi")

inputFiles = cms.untracked.vstring()
process.source = cms.Source(
    "PoolSource",
    fileNames = inputFiles
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
    )

process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.load("JetMETCorrections.JetPlusTrack.EnergyScaleHistogrammer_cfi")

from JetMETCorrections.JetPlusTrack.EnergyScaleAnalyzer_cfi import * 

process.test = energyScaleAnalyzer.clone(
    GenObjectType  = cms.string('GenJet'),
    RecoObjectType = cms.string('CaloJet'),
    GenObjectTag   = cms.InputTag('iterativeCone5GenJets'),
    RecoObjectTag  = cms.InputTag('iterativeCone5CaloJects')
    )

process.p1 = cms.Path(
    process.content
    + process.test 
    )

inputFiles.extend( [
    '/store/relval/CMSSW_2_2_1/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP_V7_v2/0001/D862836E-8EC4-DD11-92F2-001617E30D00.root',
    '/store/relval/CMSSW_2_2_1/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP_V7_v2/0001/DEB36463-8EC4-DD11-9A07-001617DBD332.root',
    '/store/relval/CMSSW_2_2_1/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP_V7_v2/0002/9AC413EC-98C4-DD11-8C80-001D09F24047.root',
    '/store/relval/CMSSW_2_2_1/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP_V7_v2/0002/CA8058FA-9FC4-DD11-8A5B-001D09F24448.root',
    '/store/relval/CMSSW_2_2_1/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP_V7_v2/0003/2809699B-FDC4-DD11-9284-000423D6AF24.root'
    ] );
