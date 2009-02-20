import FWCore.ParameterSet.Config as cms
 
patcrosscleaner = cms.EDFilter("PatCrossCleaner",
    # Main objects' sources
    patMuons      = cms.InputTag("selectedLayer1Muons"),
    patElectrons  = cms.InputTag("selectedLayer1Electrons"),
    patPhotons    = cms.InputTag("selectedLayer1Photons"),
    patTaus       = cms.InputTag("selectedLayer1Taus"),
    patJets       = cms.InputTag("selectedLayer1Jets"),
    patMets       = cms.InputTag("selectedLayer1METs"),
    # Additional objects for cross-cleaning
    patCaloTowers = cms.InputTag("towerMaker"),
    patTracks     = cms.InputTag("ctfWithMaterialTracks"),
    patVertices   = cms.InputTag("offlinePrimaryVerticesFromCTFTracks"),
    # Jet corrections (have to be redone after cleaning)
     L1JetCorrector      = cms.string('none'),
     
     L2JetCorrector      = cms.string('L2RelativeJetCorrectorIC5Calo'),
     
     L3JetCorrector      = cms.string('L3AbsoluteJetCorrectorIC5Calo'),

     L4JetCorrector      = cms.string('none'),
     
     L5udsJetCorrector   = cms.string('none'),
     L5gluonJetCorrector = cms.string('none'),
     L5cJetCorrector     = cms.string('none'),
     L5bJetCorrector     = cms.string('none'),
     
     L6JetCorrector      = cms.string('none'),
                           
     L7udsJetCorrector   = cms.string('L7PartonJetCorrectorIC5qJet'),
     L7gluonJetCorrector = cms.string('L7PartonJetCorrectorIC5gJet'),
     L7cJetCorrector     = cms.string('L7PartonJetCorrectorIC5cJet'),
     L7bJetCorrector     = cms.string('L7PartonJetCorrectorIC5bJet'),

    # Electron-jet cross-cleaning configuration
    doElectronJetCC = cms.bool(True),
    ElectronJetCrossCleaning = cms.PSet(
        SusyAnalyzerCleaning = cms.bool(True), # This is the only available at the moment anyway...
        deltaR_min        = cms.double( 0.5),
        SharedEForNIsoEle = cms.double(-1.0),
        SharedEtoJetE     = cms.double( 0.7),
        IsoValueCut       = cms.double( 1.0),
        IsolationKey      = cms.string('TrackerIso'),
	ElectronID        = cms.string('eidLoose') # Choices are: eidLoose, eidRobustHighEnergy, eid RobustLoose,
						   #              eidRobustTight, eidTight
						   # (defined in RecoEgamma/ElectronIdentification/python/electronIDSequence_cff.py
						   # No support for neural-net or likelihood methods yet
    ),
    doMuonJetCC = cms.bool(True),
    MuonJetCrossCleaning = cms.PSet(
        deltaR_min   = cms.double( 0.2),
        trackIso_max = cms.double(10.0),
        caloIso_max  = cms.double(10.0),
	MuonID = cms.string('AllGlobalMuons'), # All possible choices are defined in in DataFormats/MuonReco/interface/Muon.h
					       # See also https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookMuonAnalysis#MuonId
	modifyJetEnergy = cms.bool(True)       # option to turn of addition of muon energy to the jet for use with JPT.
    ),
    doPhotonJetCC = cms.bool(False),
    PhotonJetCrossCleaning = cms.PSet(
        deltaR_min   = cms.double(0.5),
        IsoValueCut  = cms.double(0.3),
        IsolationKey = cms.string('CaloIso'),
	PhotonID = cms.string('LooseEM') # possible choices are: LooseEM, LoosePhoton, TightPhoton or none
                                         # For further information: https://twiki.cern.ch/twiki/bin/view/CMS/PhotonIDAnalysis
    ),
    doElectronPhotonCC = cms.bool(True)
)



