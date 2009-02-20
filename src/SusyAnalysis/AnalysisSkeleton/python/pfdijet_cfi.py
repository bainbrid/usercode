import FWCore.ParameterSet.Config as cms

pfdijet = cms.EDFilter("SusyDiJetAnalysis",
    genTag = cms.InputTag("genParticles"),
    vtxTag = cms.InputTag("offlinePrimaryVertices"),                             
    tauTag = cms.InputTag("PFselectedLayer1Taus"),
    elecTag = cms.InputTag("selectedLayer1Electrons"),
    photTag = cms.InputTag("selectedLayer1Photons"),
    muonTag = cms.InputTag("PFselectedLayer1Muons"),
    jetTag = cms.InputTag("PFselectedLayer1Jets") ,
    metTag = cms.InputTag("PFselectedLayer1METs"),                             

    jptTag = cms.InputTag("patcrosscleanerJPT:ccJets"),

    ccelecTag = cms.InputTag("patcrosscleaner:ccElectrons"),                  
    ccjetTag = cms.InputTag("patcrosscleaner:ccJets"),             
    ccmuonTag = cms.InputTag("patcrosscleaner:ccMuons"),                  
    ccmetTag = cms.InputTag("patcrosscleaner:ccMETs"),
    ccphotonTag = cms.InputTag("patcrosscleaner:ccPhotons"),
  ###   selections = cms.PSet(
  #      selectionSequence = cms.vstring(
  #  'FinalDirectLeptonVeto', 
  #  'FinalMaxNumJetsSelector', 
  #  'DPhi', 
  #  'Alpha', 
  #  'Hemisphere'
  #          ),
  #      selectors = cms.PSet(
  #          FinalMaxNumJetsSelector = cms.PSet(
  #              maxEt = cms.double(30.0),
  #            #   jetTag = cms.InputTag("selectedLayer1Jets"),
  #              jetTag = cms.InputTag("PFselectedLayer1Jets"),
  #              maxNumJets = cms.uint32(100),
  #              selector = cms.string('MaxNumJetsEventSelector')
  #          ),
  #          DPhi = cms.PSet(
  #              maxDPhi = cms.double(3.15),
  #             #  jetTag = cms.InputTag("selectedLayer1Jets"),
  #             jetTag = cms.InputTag("PFselectedLayer1Jets"),
  #              selector = cms.string('DPhiEventSelector')
  #          ),
  #          Hemisphere = cms.PSet(
  #              dPhiHemispheresMin = cms.double(0.0),
  #              dPhiHemispheresMax = cms.double(3.2),
  #              hemisphereTag = cms.InputTag("pfemispheres"),
  #              selector = cms.string('HemisphereSelector')
  #          ),
  #          FinalDirectLeptonVeto = cms.PSet(
  #              #electronTag = cms.InputTag("selectedLayer1Electrons"),
  #              electronTag = cms.InputTag("selectedLayer1Electrons"),
  #              tauTag = cms.InputTag("PFselectedLayer1Taus"),
  #              minMuonEt = cms.double(30000.0),
  #              tauIsolation = cms.double(0.5),
  #              selector = cms.string('DirectLeptonVetoSelector'),
  #             muonTag = cms.InputTag("PFselectedLayer1Muons"),
  #             # muonTag = cms.InputTag("selectedLayer1Muons"),
  #              minTauEt = cms.double(30000.0),
  #              minElectronEt = cms.double(30000.0),
  #              muonIsolation = cms.double(0.5),
  #              electronIsolation = cms.double(0.5)
  #          ),
  #          Alpha = cms.PSet(
  #              minAlpha = cms.double(0.0),
  #              jetTag = cms.InputTag("PFselectedLayer1Jets"),
  #             # jetTag = cms.InputTag("selectedLayer1Jets"),
  #              selector = cms.string('AlphaSelector')
  #          )
  #      )
  #  ),
                             
                  

    pathNames = cms.vstring('HLT_Jet180','HLT_DiJetAve130','HLT_MET50','HLT_Mu9') ,    
                   
    eventWeight = cms.double(1),
  
    triggerResults = cms.InputTag("TriggerResults","","HLT"),
    plotSelection = cms.vstring(''),
                             
         
    
  
)
