import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")


process.Timing =cms.Service("Timing")

process.MessageLogger = cms.Service("MessageLogger",
                                    prod = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO')
    ),
                                    destinations = cms.untracked.vstring('prod')
                                    )





# Include PAT Layer 0 & 1 if not running on pattified data
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('IDEAL_V11::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

### Make new ECal RecHit collection ###
process.load("SUSYBSMAnalysis.CorrectedECalRecHitProducer.correctedecalrechitproducer_cfi")
process.CorrectedECalRecHitProducer.recHitsEB = cms.InputTag("ecalRecHit","EcalRecHitsEB")
process.CorrectedECalRecHitProducer.recHitsEE = cms.InputTag("ecalRecHit","EcalRecHitsEE")

## Change ECal input collection for reconstruction ###
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.towerMaker.ecalInputs = cms.VInputTag(cms.InputTag("CorrectedECalRecHitProducer","EcalRecHitsEBcorr","ANA"),
                                              cms.InputTag("CorrectedECalRecHitProducer","EcalRecHitsEEcorr","ANA"))



# CaloTowerConstituentsMap needed for Electron/Photon-Jet cleaning
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
)





process.load("PhysicsTools.PatAlgos.patLayer0_cff")

process.allLayer0METs.metSource = cms.InputTag('met','','ANA')

from PhysicsTools.PatAlgos.tools.jetTools import *

#Switching Jet collection, needed for BTag algorithm to re-run
switchJetCollection(process,
                    'iterativeCone5CaloJets',    # Jet collection; must be already in the event when patLayer0 sequence is executed
                    layers=[0,1],          # If you're not runnint patLayer1, set 'layers=[0]'
                    runCleaner="CaloJet",  # =None if not to clean
                    doJTA=True,            # Run Jet-Track association & JetCharge
                    doBTagging=True,       # Run b-tagging
                    jetCorrLabel=('IC5','Calo'), # example jet correction name; set to None for no JEC
                    doType1MET=True)       # recompute Type1 MET using these jets



process.load("PhysicsTools.PatAlgos.patLayer1_cff")

process.load("PF.Susy.patFromPF2PAT_cff")
process.load("PF.Susy.PF2PAT_cff")
process.load("PF.Susy.patLayer1_EventContent_cff")





process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     '/store/relval/CMSSW_2_2_3/RelValQCD_Pt_3000_3500/GEN-SIM-RECO/STARTUP_V7_v4/0003/243610A2-8ECB-DD11-8BAE-000423D987E0.root'
       
    )

)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)




process.p = cms.Path(process.CorrectedECalRecHitProducer*process.caloTowersRec*process.recoJets*process.metreco*process.patLayer0*process.patLayer1
                     *process.PF2PAT*process.patFromPF2PAT
                     )



process.patout = cms.OutputModule("PoolOutputModule",
    process.patPFLayer1EventContent,
    fileName = cms.untracked.string('pat.root')
)
process.outpath = cms.EndPath(  process.patout)
