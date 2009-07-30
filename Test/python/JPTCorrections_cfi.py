import FWCore.ParameterSet.Config as cms

from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi import * 

JPTCorrection = cms.PSet(

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
    UseElectrons = cms.bool(False),
    Electrons    = cms.InputTag("pixelMatchGsfElectrons"),
    ElectronIds  = cms.InputTag("electronIdTight"),

    # Use out-of-cone tracks
    AddOutOfConeTracks = cms.bool(True),

    # Filtering tracks using quality
    UseTrackQuality = cms.bool(True),
    TrackQuality    = cms.string('highPurity'),

    # Response and efficiency maps
    ResponseMap   = cms.string("JetMETCorrections/Configuration/data/CMSSW_167_response.txt"),
    EfficiencyMap = cms.string("JetMETCorrections/Configuration/data/CMSSW_167_TrackNonEff.txt"),
    LeakageMap    = cms.string("JetMETCorrections/Configuration/data/CMSSW_167_TrackLeakage.txt"),

    )
