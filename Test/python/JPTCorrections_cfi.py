import FWCore.ParameterSet.Config as cms

from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi import * 

JPTCorrection = cms.PSet(

    # General Configuration
    Verbose           = cms.bool(False),
    UsePatCollections = cms.bool(False),
    
    # Select correction types
    UseInConeTracks      = cms.bool(True),
    UseOutOfConeTracks   = cms.bool(True),
    UseOutOfVertexTracks = cms.bool(True),
    
    # Jet-tracks association (null value = "on-the-fly" mode)
    JetTracksAssociationAtVertex   = cms.InputTag("JetTracksAssociationAtVertex"), 
    JetTracksAssociationAtCaloFace = cms.InputTag("JetTracksAssociationAtCaloFace"),
    
    # Jet-tracks association "on-the-fly" mode
    AllowOnTheFly = cms.bool(False),
    Tracks        = cms.InputTag("generalTracks"),
    Propagator    = cms.string('SteppingHelixPropagatorAlong'),
    ConeSize      = cms.double(0.5),
    
    # Muons
    UseMuons = cms.bool(True),
    Muons    = cms.InputTag("muons"),
    
    # Electrons
    UseElectrons = cms.bool(False),
    Electrons    = cms.InputTag("pixelMatchGsfElectrons"),
    ElectronIds  = cms.InputTag("electronIdTight"),
    
    # Filtering tracks using quality
    UseTrackQuality = cms.bool(True),
    TrackQuality    = cms.string('highPurity'),
    
    # Response and efficiency maps
    ResponseMap   = cms.string("JetMETCorrections/Configuration/data/CMSSW_167_response.txt"),
    EfficiencyMap = cms.string("JetMETCorrections/Configuration/data/CMSSW_167_TrackNonEff.txt"),
    LeakageMap    = cms.string("JetMETCorrections/Configuration/data/CMSSW_167_TrackLeakage.txt"),

    )
