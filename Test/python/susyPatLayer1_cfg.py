import FWCore.ParameterSet.Config as cms

# 1: GMSB,22X,Summer08Redigi 
# 2: PhotonJets,21X,Summer08Redigi 
# 3: PhotonJets,22X,Summer08Redigi(Fall08) 
Vers = str("1")

if Vers == str("1") :
    Source  = str("GMSB") # Options: "GMSB", "PhotonJets"
    Release = str("22X") # Options: "22X", "21X"
    Corrs   = str("Summer08Redigi") # Options: "Summer08Redigi", "Summer08", "Winter09"
elif Vers == str("2") :
    Source  = str("PhotonJets") # Options: "GMSB", "PhotonJets"
    Release = str("21X") # Options: "22X", "21X"
    Corrs   = str("Summer08Redigi") # Options: "Summer08Redigi", "Summer08", "Winter09"
elif Vers == str("3") :
    Source  = str("PhotonJets") # Options: "GMSB", "PhotonJets"
    Release = str("22X") # Options: "22X", "21X"
    Corrs   = str("Summer08Redigi") # Options: "Summer08Redigi", "Summer08", "Winter09"
else :
    print "UNKNOWN VERSION!"

Name = str(Source + "_" + Release + "_" + Corrs)

process = cms.Process("PAT")

process.Timing =cms.Service("Timing")
process.Tracer = cms.Service(
    "Tracer",
    sourceSeed = cms.untracked.string("$$")
    )

process.load("DQM.SiStripCommon.MessageLogger_cfi")
process.MessageLogger.debugModules = cms.untracked.vstring('')
process.MessageLogger.destinations = cms.untracked.vstring(
    "cerr", "susyPatLayer1_" + Name, "info", "warning", "error"
    )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

inputFiles = cms.untracked.vstring()
if Source == str("GMSB") :
    inputFiles.extend([
        'file:/home/bainbrid/data/susy/photonid/data/Summer08.Exotica_GMSB_GM1g.GEN-SIM-RECO.IDEAL_V11_redigi_v1.root'
        ])
elif Source == str("PhotonJets") :
    if Release == str("21X") :
        inputFiles.extend([
            'file:/home/bainbrid/data/susy/photonid/data/Fall08.PhotonJets200toInf-madgraph.GEN-SIM-RECO.IDEAL_V9_reco-v1.root'
            ])
    elif Release == str("22X") :
        inputFiles.extend([
            'file:/data/tom/susysamples/fall08/PhotonJetMadGraph/GEN-SIM-RECO/PhotonJet200toInfMadGraph_221_GEN-SIM-RECO_block0001.root'
            ])
    else :
        print "UNKNOWN SAMPLE AND RELEASE!"
else :
    print "UNKNOWN SAMPLE!"

if Corrs == str("Summer08Redigi") :
    process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08Redigi_cff")
elif Corrs == str("Summer08") :
    process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08_cff")
elif Corrs == str("Winter09") :
    process.load("JetMETCorrections.Configuration.L2L3Corrections_Winter09_cff")
else :
    print "UNKNOWN CORRECTIONS!"

process.source = cms.Source(
    "PoolSource",
    fileNames = inputFiles
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.load("SUSYBSMAnalysis.CorrectedECalRecHitProducer.correctedecalrechitproducer_cfi")
process.CorrectedECalRecHitProducer.recHitsEB = cms.InputTag("ecalRecHit","EcalRecHitsEB")
process.CorrectedECalRecHitProducer.recHitsEE = cms.InputTag("ecalRecHit","EcalRecHitsEE")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.towerMaker.ecalInputs = cms.VInputTag(
    cms.InputTag("CorrectedECalRecHitProducer","EcalRecHitsEBcorr","PAT"),
    cms.InputTag("CorrectedECalRecHitProducer","EcalRecHitsEEcorr","PAT")
    )
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.CaloTowerConstituentsMapBuilder = cms.ESProducer(
    "CaloTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
    )

process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")

from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(
    process,
    'iterativeCone5CaloJets',
    layers=[0,1],
    runCleaner="CaloJet",
    doJTA=True,
    doBTagging=True,
    jetCorrLabel=('IC5','Calo'),
    doType1MET=True
    )

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
    )

# Rerun PhotonID producer using RECO 
from RecoEgamma.PhotonIdentification.photonId_cfi import PhotonIDProd
process.RegeneratedPhotonID = PhotonIDProd.clone()
process.RegeneratedPhotonID.photonIDAssociationLabel = cms.string('RegeneratedPhotonAssociatedID')
process.patAODPhotonID.photonID = cms.InputTag("RegeneratedPhotonID","RegeneratedPhotonAssociatedID")

# Rerun PhotonID producer using PAT 
process.load("bainbrid.Test.patPhotonIDProducer_cff")

process.common = cms.Sequence(
    process.RegeneratedPhotonID * 
    process.CorrectedECalRecHitProducer*
    process.caloTowersRec*
    process.recoJets*
    process.metreco*
    process.patLayer0*
    process.patLayer1
    )
    
if Release == str("22X") :
    process.GlobalTag.globaltag = cms.string('IDEAL_V11::All')
    process.load("PF.Susy.patFromPF2PAT_cff")
    process.load("PF.Susy.PF2PAT_cff")
    process.p = cms.Path( 
        process.common * 
        process.PF2PAT*
        process.patFromPF2PAT *
        process.patPhotonIDProducer
        )
elif Release == str("21X") :
    process.GlobalTag.globaltag = cms.string('IDEAL_V9::All')
    from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
    run22XonSummer08AODSIM(process)
    process.p = cms.Path( 
        process.common * 
        process.patPhotonIDProducer
        )
else :
    print "UNKNOWN RELEASE!"
    
process.load("PhysicsTools.PatAlgos.patLayer0_EventContent_cff")
process.load("PF.Susy.patLayer1_EventContent_cff")
process.output = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring("keep *"),
    #process.patPFLayer1EventContent, 
    fileName = cms.untracked.string("file:/tmp/" + "susyPatLayer1_" + Name + ".root")
   )
#process.output.outputCommands += process.patLayer0EventContent.outputCommands

#process.output.outputCommands.extend(['keep patPhotons_patPhotonIDProducer_*_*'])
#process.output.outputCommands.extend(['keep recoTrackExtras_generalTracks_*_*'])
#process.output.outputCommands.extend(['keep recoSuperClusters_*_*_*'])
#process.output.outputCommands.extend(['keep recoGenParticles_genParticles_*_*'])
#process.output.outputCommands.extend(['keep recoPhotonIDs_PhotonIDProd_PhotonAssociatedID_*'])
process.e = cms.EndPath(process.output)
