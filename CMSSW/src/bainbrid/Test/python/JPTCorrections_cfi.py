import FWCore.ParameterSet.Config as cms

from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi import * 

JPTCorrection = cms.PSet(

    # General Configuration
    Verbose = cms.bool(False),

    # Vectorial corrections
    VectorialCorrection = cms.bool(True),
    JetDirFromTracks    = cms.bool(True),
    
    # Select tracks used in correction
    UseInConeTracks      = cms.bool(True),
    UseOutOfConeTracks   = cms.bool(True),
    UseOutOfVertexTracks = cms.bool(True),
    
    # Jet-tracks association
    JetTracksAssociationAtVertex   = cms.InputTag("JetTracksAssociationAtVertex"), 
    JetTracksAssociationAtCaloFace = cms.InputTag("JetTracksAssociationAtCaloFace"),
    
    # Jet merging/splitting
    JetSplitMerge = cms.int32(0),

    # Pions
    UsePions      = cms.bool(True),
    UseEfficiency = cms.bool(True),
    
    # Muons
    UseMuons = cms.bool(True),
    Muons    = cms.InputTag("muons"),
    
    # Electrons
    UseElectrons    = cms.bool(True),
    Electrons       = cms.InputTag("gsfElectrons"),
    ElectronIds     = cms.InputTag("eleIdTight"),
    
    # Filtering tracks using quality
    UseTrackQuality = cms.bool(True),
    TrackQuality    = cms.string('highPurity'),
    
    # Response and efficiency maps
    ResponseMap   = cms.string("JetMETCorrections/Configuration/data/CMSSW_167_response.txt"),
    EfficiencyMap = cms.string("JetMETCorrections/Configuration/data/CMSSW_167_TrackNonEff.txt"),
    LeakageMap    = cms.string("JetMETCorrections/Configuration/data/CMSSW_167_TrackLeakage.txt"),

    )
