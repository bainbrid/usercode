import FWCore.ParameterSet.Config as cms

process = cms.Process("JEC")

process.load("DQM.SiStripCommon.MessageLogger_cfi")
#process.MessageLogger.debugModules = '' # comment to allow LogDebug
Suppress = cms.untracked.vstring(
    "Geometry", "GeometryConfiguration", "DDLParser",
    "TrackerGeom", "PixelGeom",
    "TIBGeom", "TIDGeom", "TOBGeom", "TECGeom",
    "EcalGeom", "SFGeom", "HCalGeom"
    )
process.MessageLogger.suppressError = Suppress
process.MessageLogger.suppressDebug = Suppress
process.MessageLogger.suppressInfo = Suppress

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FakeConditions_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('STARTUP_V4::All')

process.load("JetMETCorrections.JetPlusTrack.EnergyScaleHistogrammer_cfi")

inputFiles21X = cms.untracked.vstring()
inputFiles22X = cms.untracked.vstring()
process.source = cms.Source(
    "PoolSource",
    fileNames = inputFiles22X
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
    )

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.genParticlesForJets.ignoreParticleIDs = cms.vuint32(
    12, 14, 16, 39, 1000022, 1000039, 2000012, 2000014, 2000016, 4000012, 5000039, 9900012, 9900014, 9900016
    )
process.genParticlesForJets.excludeFromResonancePids = cms.vuint32(
    12, 14, 16
    )
from RecoJets.JetProducers.iterativeCone5GenJets_cff import iterativeCone5GenJets
process.iterativeCone5GenJetsNoNuBSM = iterativeCone5GenJets.clone()

process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.correctedJetSeq = cms.Sequence(
    process.genParticlesForJets * 
    process.iterativeCone5GenJetsNoNuBSM 
    #+ process.content 
    )

process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('CorrectedJets.root'),
    outputCommands = cms.untracked.vstring(
    'drop *',
    'keep recoGenJets_iterativeCone5GenJets_*_*',
    'keep recoGenJets_iterativeCone5GenJetsNoNuBSM_*_*',
    'keep recoCaloJets_iterativeCone5CaloJets_*_*',
    )
    )

inputFiles22X.extend( [
    '/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0000/048342F6-8ACB-DD11-A608-001D0967D616.root',
    '/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0000/18F38FAD-6BCB-DD11-B56F-0019B9E4FC0D.root',
    '/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0000/44F46DAC-6BCB-DD11-BACD-001D0967C0A9.root',
    '/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0000/789EFAAA-6BCB-DD11-99C2-0019B9E13D4D.root',
    '/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0000/DAEA3773-4ECB-DD11-9A20-001D0967D062.root',
    '/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0001/7E9256AA-ACCB-DD11-90F2-0019B9E71500.root',
    '/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0001/A810340A-EACB-DD11-B8ED-001D0967A19D.root',
    '/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0001/F6D9AA91-F1CC-DD11-A364-001D0967CE96.root',
    '/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0003/32DE4C7A-F8CD-DD11-86C9-0019B9E48FFC.root',
    '/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0007/14D7A555-BBCE-DD11-92D6-0019B9E50063.root',
    '/store/mc/Summer08/QCDDiJetPt80to120/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0007/706E94B7-7CCE-DD11-82EA-0019B9E70855.root',
    ])

inputFiles21X.extend( [
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/000AD2A4-6E86-DD11-AA99-000423D9863C.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/02D641CC-6D86-DD11-B1AA-001617C3B64C.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/0876950D-6F86-DD11-B4B8-000423D6A6F4.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/0AC366B7-6D86-DD11-8695-001617E30D0A.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/1415D0A5-6E86-DD11-8F0C-001617C3B778.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/145FBEB9-6D86-DD11-BD98-001617DF785A.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/1815C2A0-6E86-DD11-A160-001617E30F58.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/1869C2C0-6D86-DD11-ADB6-001617C3B76E.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/186A3FC0-6D86-DD11-B020-000423D6CA72.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/1A688BC9-6D86-DD11-8C60-001617DBD5AC.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/20EBA2C9-6D86-DD11-88EC-001617E30D40.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/24822ABB-6D86-DD11-B76B-001617E30F4C.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/24E0B0C7-6D86-DD11-8C09-0019DB29C5FC.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/2652503A-6F86-DD11-9552-001617E30F4C.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/26E79ACA-6D86-DD11-8A10-001617E30D00.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/2A0341A9-6E86-DD11-9AE4-000423D6CA02.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/3038D377-6E86-DD11-A974-001617C3B5E4.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/369AEBC5-6D86-DD11-B911-001617E30D4A.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/3841A6BC-6D86-DD11-9D84-001617E30D52.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/3C049BBF-6D86-DD11-8D0D-001617DBD556.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/3EBAB9A6-6E86-DD11-AB4D-001617C3B70E.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/3EBCA837-6F86-DD11-8892-001617E30D4A.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/3EDB9AD6-6E86-DD11-9EE3-000423D992A4.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/42CED2B7-6D86-DD11-B551-001617C3B6E2.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/441FDE81-6D86-DD11-863A-000423D985E4.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/48659740-6F86-DD11-926F-001617E30D12.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/488BFB15-6F86-DD11-B4A2-000423D6C8EE.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/4AA2ACC0-6D86-DD11-8CC9-001617E30E2C.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/5619FDAA-6E86-DD11-85B7-000423D94700.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/5A1FDE7E-6E86-DD11-BAFF-000423D6AF24.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/5A8019BB-6D86-DD11-80F7-001617C3B78C.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/64769ED1-6D86-DD11-9684-000423D6CAF2.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/64A80339-6F86-DD11-8266-001617C3B73A.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/66F1CAA9-6E86-DD11-8BE7-000423D992DC.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/6E70BEBC-6D86-DD11-97A0-001617E30CA4.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/7479B5CA-6D86-DD11-8901-001617E30D12.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/7A9AB3C3-6D86-DD11-BCA2-001617E30CC8.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/7E668FC2-6D86-DD11-8260-001617C3B6E8.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/801375AA-6E86-DD11-97AF-000423D9880C.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/801C1FA2-6E86-DD11-BAE4-001617C3B69C.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/808501C6-6D86-DD11-AB39-001617E30E28.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/80BEEB00-6F86-DD11-8777-0016177CA7A0.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/82C51735-6F86-DD11-8646-001617DBD230.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/8638B0BD-6D86-DD11-99C5-001617DBD230.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/8668573A-6F86-DD11-BB1F-0019DB29C5FC.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/8687E8C2-6D86-DD11-AA19-0019DB2F3F9B.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/86A60FB4-6D86-DD11-B3CD-001617DBD316.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/9A32E433-6F86-DD11-B6BF-001617DF785A.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/9C65C37A-6E86-DD11-A876-000423D6B5C4.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/9E9705BA-6D86-DD11-81E4-001617DC1F70.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/A26163C0-6D86-DD11-984A-001617E30CE8.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/A8CAC58A-6D86-DD11-9805-000423D992A4.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/A8D988C2-6D86-DD11-8342-001617E30F48.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/AA49AFC0-6D86-DD11-9D71-001617DBD5B2.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/AED3FD3C-6F86-DD11-A9F7-001617C3B6CC.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/B2AFB941-6F86-DD11-82E7-001617C3B77C.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/BA7B56B9-6D86-DD11-957E-001617E30D38.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/BAA819DC-6E86-DD11-8139-000423D98DB4.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/BC325410-6F86-DD11-9E26-000423D6CAF2.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/C02991CA-6D86-DD11-A4A6-001617E30D06.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/CAF574C0-6D86-DD11-8D8E-001617DBCF90.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/D06ECFC8-6D86-DD11-A19B-001617C3B66C.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/D2DB72C0-6D86-DD11-BFBE-001617C3B5F4.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/D4700338-6F86-DD11-9F9D-001617C3B6DC.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/D67B9ABD-6D86-DD11-A6A2-001617DBD288.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/E47BD0AA-6E86-DD11-B89E-000423DD2F34.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/EA298EBC-6D86-DD11-9FAD-001617C3B73A.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/EEE12E83-6D86-DD11-835A-000423D9880C.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/F0875383-6E86-DD11-BF96-001617E30F50.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/F0A4EEA2-6E86-DD11-B683-001617C3B710.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/F2A287BF-6D86-DD11-8DD6-001617C3B77C.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/F47E188D-6D86-DD11-B87C-000423D992DC.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/F66E700A-6F86-DD11-B510-000423D6CA42.root',
    '/store/relval/CMSSW_2_1_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0001/FC15908F-6D86-DD11-AF1B-000423D9853C.root'
    ] );
