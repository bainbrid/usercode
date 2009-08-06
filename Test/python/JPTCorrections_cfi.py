import FWCore.ParameterSet.Config as cms

from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi import * 

JPTCorrection = cms.PSet(

    # General Configuration
    Verbose              = cms.bool(False),
    UseInConeTracks      = cms.bool(True),
    UseOutOfConeTracks   = cms.bool(True),
    UseOutOfVertexTracks = cms.bool(True),
    UseMuons             = cms.bool(True),
    UseElectrons         = cms.bool(True),
    UseTrackQuality      = cms.bool(True),

    # Jet-tracks association (null value = "on-the-fly" mode)
    JetTracksAssociationAtVertex   = cms.InputTag("JetTracksAssociationAtVertex"), 
    JetTracksAssociationAtCaloFace = cms.InputTag("JetTracksAssociationAtCaloFace"),
    
    # Jet-tracks association "on-the-fly" mode
    Tracks     = cms.InputTag("generalTracks"),
    Propagator = cms.string('SteppingHelixPropagatorAlong'),
    ConeSize   = cms.double(0.5),

    # Muons
    Muons = cms.InputTag("muons"),
    
    # Electrons
    Electrons    = cms.InputTag("pixelMatchGsfElectrons"),
    ElectronIds  = cms.InputTag("electronIdTight"),
    
    # Filtering tracks using quality
    TrackQuality    = cms.string('highPurity'),

    # Response and efficiency maps
    ResponseMap   = cms.string("JetMETCorrections/Configuration/data/CMSSW_167_response.txt"),
    EfficiencyMap = cms.string("JetMETCorrections/Configuration/data/CMSSW_167_TrackNonEff.txt"),
    LeakageMap    = cms.string("JetMETCorrections/Configuration/data/CMSSW_167_TrackLeakage.txt"),

    )
