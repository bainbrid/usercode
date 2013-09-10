import FWCore.ParameterSet.Config as cms

pfdijet = cms.EDFilter(
    "SusyDiJetAnalysis",

    genTag = cms.InputTag("genParticles"),
    vtxTag = cms.InputTag("offlinePrimaryVertices"),                             
    tauTag = cms.InputTag("PFselectedLayer1Taus"),
    elecTag = cms.InputTag("selectedLayer1Electrons"),
    photTag = cms.InputTag("selectedLayer1Photons"),
    muonTag = cms.InputTag("PFselectedLayer1Muons"),
    jetTag = cms.InputTag("PFselectedLayer1Jets") ,
    metTag = cms.InputTag("PFselectedLayer1METs"),                             

    ccelecTag = cms.InputTag("patcrosscleaner:ccElectrons"),                  
    ccjetTag = cms.InputTag("patcrosscleaner:ccJets"),             
    ccmuonTag = cms.InputTag("patcrosscleaner:ccMuons"),                  
    ccmetTag = cms.InputTag("patcrosscleaner:ccMETs"),
    ccphotonTag = cms.InputTag("patcrosscleaner:ccPhotons"),

    jptTag = cms.InputTag("patcrosscleanerJPT:ccJets"),

    pathNames = cms.vstring('HLT_Jet180','HLT_DiJetAve130','HLT_MET50','HLT_Mu9') ,    
                   
    eventWeight = cms.double(1.),
    
    triggerResults = cms.InputTag("TriggerResults","","HLT"),
    plotSelection = cms.vstring(''),
    
    )
