import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.Timing = cms.Service("Timing")
process.Tracer = cms.Service("Tracer",sourceSeed = cms.untracked.string("$$"))

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP31X_V4::All')

# Input files: RelVal QCD 80-120 GeV, STARTUP conditions, 9000 events, from CMSSW_3_2_5 (replace with 33X when available!)
fileNames1 = cms.untracked.vstring()
fileNames2 = cms.untracked.vstring()
process.source = cms.Source(
    "PoolSource",
    fileNames = fileNames2,
    )
fileNames1.extend(
    [ '/store/relval/CMSSW_3_2_5/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP31X_V4-v1/0011/2AA47C1E-828E-DE11-B3C5-001D09F34488.root',
      '/store/relval/CMSSW_3_2_5/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP31X_V4-v1/0010/DA6FF61D-3D8E-DE11-938A-003048D37538.root',
      '/store/relval/CMSSW_3_2_5/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP31X_V4-v1/0010/82CF8666-398E-DE11-8F3B-000423D94A20.root',
      '/store/relval/CMSSW_3_2_5/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP31X_V4-v1/0010/6255E85C-3F8E-DE11-B46A-000423D6B48C.root',
      '/store/relval/CMSSW_3_2_5/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP31X_V4-v1/0010/3453CD32-418E-DE11-87D2-003048D2C020.root',
      '/store/relval/CMSSW_3_2_5/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP31X_V4-v1/0010/2EC02533-3B8E-DE11-BE85-003048D37514.root',
      ] );
fileNames2.extend(
    [ 'file:/home/bainbrid/work/src/TEST/2AA47C1E-828E-DE11-B3C5-001D09F34488.root',
      'file:/home/bainbrid/work/src/TEST/DA6FF61D-3D8E-DE11-938A-003048D37538.root',
      'file:/home/bainbrid/work/src/TEST/82CF8666-398E-DE11-8F3B-000423D94A20.root',
      'file:/home/bainbrid/work/src/TEST/6255E85C-3F8E-DE11-B46A-000423D6B48C.root',
      'file:/home/bainbrid/work/src/TEST/3453CD32-418E-DE11-87D2-003048D2C020.root',
      'file:/home/bainbrid/work/src/TEST/2EC02533-3B8E-DE11-BE85-003048D37514.root',
      ] );

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# ---------- Old and new JPT from RECO ----------

process.load("JetMETCorrections.Configuration.ZSPJetCorrections219_cff")
process.load("JetMETCorrections.Configuration.JetPlusTrackCorrections_cff")
process.JetPlusTrackZSPCorrectorIcone5.ElectronIds = cms.InputTag("eidTight")

process.load("bainbrid.Test.JPTCorrections_cff")
process.JPTCorrectionIC5.ElectronIds = cms.InputTag("eidTight")
process.JPTCorrectionIC5.VectorialCorrection = cms.bool(False)

# ---------- New JPT from RECO, all steps ----------

# None (no correction)

from bainbrid.Test.JPTCorrections_cff import *

process.JPTCorrectionNone = JPTCorrectionIC5.clone()
process.JPTCorrectionNone.label = cms.string('JPTCorrectionNone')
#process.JPTCorrectionNone.Verbose = cms.bool(True)

process.JPTCorrectionNone.ElectronIds = cms.InputTag("eidTight")

process.JPTCorrectionNone.UseInConeTracks        = cms.bool(False)
process.JPTCorrectionNone.UseOutOfConeTracks     = cms.bool(False)
process.JPTCorrectionNone.UseOutOfVertexTracks   = cms.bool(False)
process.JPTCorrectionNone.UseEfficiency          = cms.bool(False)
process.JPTCorrectionNone.UseMuons               = cms.bool(False)
process.JPTCorrectionNone.UseElectrons           = cms.bool(False)
process.JPTCorrectionNone.VectorialCorrection    = cms.bool(False)
process.JPTCorrectionNone.VecCorrUsingTracksOnly = cms.bool(False)

process.JPTCorrectorNone = JPTCorrectorIC5.clone()
process.JPTCorrectorNone.correctors = cms.vstring('JPTCorrectionNone')
process.JPTCorrectorNone.alias = cms.untracked.string('JPTCorrectorNone')

# + InCone (pions only)

from bainbrid.Test.JPTCorrections_cff import *

process.JPTCorrectionInCone = process.JPTCorrectionNone.clone()
process.JPTCorrectionInCone.label = cms.string('JPTCorrectionInCone')
process.JPTCorrectionInCone.UseInConeTracks = cms.bool(True)

process.JPTCorrectorInCone = process.JPTCorrectorIC5.clone()
process.JPTCorrectorInCone.correctors = cms.vstring('JPTCorrectionInCone')
process.JPTCorrectorInCone.alias = cms.untracked.string('JPTCorrectorInCone')

# + OutOfCone (pions only)

from bainbrid.Test.JPTCorrections_cff import *

process.JPTCorrectionOutOfCone = process.JPTCorrectionInCone.clone()
process.JPTCorrectionOutOfCone.label = cms.string('JPTCorrectionOutOfCone')
process.JPTCorrectionOutOfCone.UseOutOfConeTracks = cms.bool(True)

process.JPTCorrectorOutOfCone = process.JPTCorrectorIC5.clone()
process.JPTCorrectorOutOfCone.correctors = cms.vstring('JPTCorrectionOutOfCone')
process.JPTCorrectorOutOfCone.alias = cms.untracked.string('JPTCorrectorOutOfCone')

# + OutOfVertex (pions only)

from bainbrid.Test.JPTCorrections_cff import *

process.JPTCorrectionOutOfVertex = process.JPTCorrectionOutOfCone.clone()
process.JPTCorrectionOutOfVertex.label = cms.string('JPTCorrectionOutOfVertex')
process.JPTCorrectionOutOfVertex.UseOutOfVertexTracks = cms.bool(True)

process.JPTCorrectorOutOfVertex = process.JPTCorrectorIC5.clone()
process.JPTCorrectorOutOfVertex.correctors = cms.vstring('JPTCorrectionOutOfVertex')
process.JPTCorrectorOutOfVertex.alias = cms.untracked.string('JPTCorrectorOutOfVertex')

# + PionEff

from bainbrid.Test.JPTCorrections_cff import *

process.JPTCorrectionPionEff = process.JPTCorrectionOutOfVertex.clone()
process.JPTCorrectionPionEff.label = cms.string('JPTCorrectionPionEff')
process.JPTCorrectionPionEff.UseEfficiency = cms.bool(True)

process.JPTCorrectorPionEff = process.JPTCorrectorIC5.clone()
process.JPTCorrectorPionEff.correctors = cms.vstring('JPTCorrectionPionEff')
process.JPTCorrectorPionEff.alias = cms.untracked.string('JPTCorrectorPionEff')

# + Muons

from bainbrid.Test.JPTCorrections_cff import *

process.JPTCorrectionMuons = process.JPTCorrectionPionEff.clone()
process.JPTCorrectionMuons.label = cms.string('JPTCorrectionMuons')
process.JPTCorrectionMuons.UseMuons = cms.bool(True)

process.JPTCorrectorMuons = process.JPTCorrectorIC5.clone()
process.JPTCorrectorMuons.correctors = cms.vstring('JPTCorrectionMuons')
process.JPTCorrectorMuons.alias = cms.untracked.string('JPTCorrectorMuons')

# + Electrons

from bainbrid.Test.JPTCorrections_cff import *

process.JPTCorrectionElectrons = process.JPTCorrectionMuons.clone()
process.JPTCorrectionElectrons.label = cms.string('JPTCorrectionElectrons')
process.JPTCorrectionElectrons.UseElectrons = cms.bool(True)

process.JPTCorrectorElectrons = process.JPTCorrectorIC5.clone()
process.JPTCorrectorElectrons.correctors = cms.vstring('JPTCorrectionElectrons')
process.JPTCorrectorElectrons.alias = cms.untracked.string('JPTCorrectorElectrons')

# + VectorialCorrection

from bainbrid.Test.JPTCorrections_cff import *

process.JPTCorrectionVectorial = process.JPTCorrectionElectrons.clone()
process.JPTCorrectionVectorial.label = cms.string('JPTCorrectionVectorial')
process.JPTCorrectionVectorial.VectorialCorrection = cms.bool(True)

process.JPTCorrectorVectorial = process.JPTCorrectorIC5.clone()
process.JPTCorrectorVectorial.correctors = cms.vstring('JPTCorrectionVectorial')
process.JPTCorrectorVectorial.alias = cms.untracked.string('JPTCorrectorVectorial')

# + VecCorrUsingTracksOnly

from bainbrid.Test.JPTCorrections_cff import *

process.JPTCorrectionVecTracks = process.JPTCorrectionVectorial.clone()
process.JPTCorrectionVecTracks.label = cms.string('JPTCorrectionVecTracks')
process.JPTCorrectionVecTracks.VecCorrUsingTracksOnly = cms.bool(True)

process.JPTCorrectorVecTracks = process.JPTCorrectorIC5.clone()
process.JPTCorrectorVecTracks.correctors = cms.vstring('JPTCorrectionVecTracks')
process.JPTCorrectorVecTracks.alias = cms.untracked.string('JPTCorrectorVecTracks')


process.JPTCorrectionNone.VecCorrUsingTracksOnly = cms.bool(False)

# Sequences

process.JPTValidation = cms.Sequence(
    process.JPTCorrectorNone *
    process.JPTCorrectorInCone *
    process.JPTCorrectorOutOfCone *
    process.JPTCorrectorOutOfVertex *
    process.JPTCorrectorPionEff *
    process.JPTCorrectorMuons *
    process.JPTCorrectorElectrons * 
    process.JPTCorrectorVectorial *
    process.JPTCorrectorVecTracks 
    )

# ---------- Raw PAT Jets ----------

process.RawPATJets = cms.EDProducer(
    "RawPATJetProducer",
    JetCollection = cms.InputTag("cleanLayer1Jets"), 
    )

# ---------- ZSP from PAT ----------

process.PATZSPJetCorrectorIcone5 = cms.ESSource(
    "ZSPJetCorrectionService",
    tagName = cms.string('ZSP_CMSSW219_Iterative_Cone_05'),
    label = cms.string('PATZSPJetCorrectorIcone5'),
    )

process.PATZSPJetCorJetIcone5 = cms.EDProducer(
    "PatJetCorrectionProducer",
    src = cms.InputTag("RawPATJets"),
    correctors = cms.vstring('PATZSPJetCorrectorIcone5'),
    alias = cms.untracked.string('PATZSPJetCorJetIcone5'),
    )

# ---------- JTA from PAT ----------

process.load("JetMETCorrections.Configuration.JetCorrectionsRecord_cfi")
process.load("RecoJets.Configuration.RecoJetAssociations_cff")
process.load("RecoJets.JetAssociationProducers.iterativeCone5JTA_cff")

process.PATJtaAtVertex = process.iterativeCone5JetTracksAssociatorAtVertex.clone() 
process.PATJtaAtVertex.jets = cms.InputTag("PATZSPJetCorJetIcone5")

process.PATJtaAtCaloFace = process.iterativeCone5JetTracksAssociatorAtCaloFace.clone()
process.PATJtaAtCaloFace.jets = cms.InputTag("PATZSPJetCorJetIcone5")

process.PATJtaExtender = process.iterativeCone5JetExtender.clone() 
process.PATJtaExtender.jets = cms.InputTag("PATZSPJetCorJetIcone5")
process.PATJtaExtender.jet2TracksAtCALO = cms.InputTag("PATJtaAtCaloFace")
process.PATJtaExtender.jet2TracksAtVX = cms.InputTag("PATJtaAtVertex")

process.PATJetTrackAssociations = cms.Sequence(
    process.PATJtaAtVertex *
    process.PATJtaAtCaloFace *
    process.PATJtaExtender
    )

# ---------- Old JPT from PAT ----------

from JetMETCorrections.Configuration.JetPlusTrackCorrections_cfi import *

process.PATJetPlusTrackZSPCorrectorIcone5 = cms.ESSource(
    "JetPlusTrackCorrectionService",
    JPTZSPCorrectorICone5,
    label = cms.string('PATJetPlusTrackZSPCorrectorIcone5'),
    )

process.PATJetPlusTrackZSPCorJetIcone5 = cms.EDProducer(
    "PatJetCorrectionProducer",
    src = cms.InputTag("PATZSPJetCorJetIcone5"),
    correctors = cms.vstring('PATJetPlusTrackZSPCorrectorIcone5'),
    alias = cms.untracked.string('PATJetPlusTrackZSPCorJetIcone5'),
)

process.PATJetPlusTrackZSPCorrectorIcone5.UsePatCollections = cms.bool(True)
process.PATJetPlusTrackZSPCorrectorIcone5.Muons = cms.InputTag("cleanLayer1Muons")
process.PATJetPlusTrackZSPCorrectorIcone5.Electrons = cms.InputTag("cleanLayer1Electrons")
process.PATJetPlusTrackZSPCorrectorIcone5.ElectronIds = cms.InputTag("eidTight")
process.PATJetPlusTrackZSPCorrectorIcone5.JetTracksAssociationAtVertex = cms.InputTag("PATJtaAtVertex")
process.PATJetPlusTrackZSPCorrectorIcone5.JetTracksAssociationAtCaloFace = cms.InputTag("PATJtaAtCaloFace")
#process.PATJetPlusTrackZSPCorrectorIcone5.Verbose = cms.bool(True)

process.PATJetPlusTracksCorrections = cms.Sequence(
    process.RawPATJets *
    process.PATZSPJetCorJetIcone5 *
    process.PATJetTrackAssociations *
    process.PATJetPlusTrackZSPCorJetIcone5
    )

# ---------- New JPT from PAT ----------

from bainbrid.Test.JPTCorrections_cfi import *

process.PATJPTCorrectorIC5 = cms.ESSource(
    "JPTPatCorrectionService",
    JPTCorrection,
    # Derived
    AllowOnTheFly     = cms.bool(True),
    Tracks            = cms.InputTag("generalTracks"),
    Propagator        = cms.string('SteppingHelixPropagatorAlong'),
    ConeSize          = cms.double(0.5),
    UsePatCollections = cms.bool(True),
    label = cms.string('PATJPTCorrectorIC5'),
    )

process.PATJPTCorrectionIC5 = cms.EDProducer(
    "PatJetCorrectionProducer",
    src = cms.InputTag("PATZSPJetCorJetIcone5"),
    correctors = cms.vstring('PATJPTCorrectorIC5'),
    alias = cms.untracked.string('PATJPTCorrectionIC5'),
)

process.PATJPTCorrectorIC5.Muons = cms.InputTag("cleanLayer1Muons")
process.PATJPTCorrectorIC5.Electrons = cms.InputTag("cleanLayer1Electrons")
process.PATJPTCorrectorIC5.ElectronIds = cms.InputTag("eidTight")
process.PATJPTCorrectorIC5.JetTracksAssociationAtVertex = cms.InputTag("PATJtaAtVertex")
process.PATJPTCorrectorIC5.JetTracksAssociationAtCaloFace = cms.InputTag("PATJtaAtCaloFace")
#process.PATJPTCorrectorIC5.Verbose = cms.bool(True)

process.PATJPTCorrections = cms.Sequence(
    process.RawPATJets *
    process.PATZSPJetCorJetIcone5 *
    process.PATJetTrackAssociations *
    process.PATJPTCorrectionIC5
    )

# ---------- New JPT with JTA on-the-fly ----------

from bainbrid.Test.JPTCorrections_cfi import *

process.JTAJPTCorrectorIC5 = cms.ESSource(
    "JPTPatCorrectionService",
    JPTCorrection,
    # New
    AllowOnTheFly     = cms.bool(True),
    Tracks            = cms.InputTag("generalTracks"),
    Propagator        = cms.string('SteppingHelixPropagatorAlong'),
    ConeSize          = cms.double(0.5),
    UsePatCollections = cms.bool(True),
    label = cms.string('JTAJPTCorrectorIC5'),
    )

process.JTAJPTCorrectionIC5 = cms.EDProducer(
    "PatJetCorrectionProducer",
    src = cms.InputTag("PATZSPJetCorJetIcone5"),
    correctors = cms.vstring('JTAJPTCorrectorIC5'),
    alias = cms.untracked.string('JTAJPTCorrectionIC5'),
)

process.JTAJPTCorrectorIC5.Muons = cms.InputTag("cleanLayer1Muons")
process.JTAJPTCorrectorIC5.Electrons = cms.InputTag("cleanLayer1Electrons")
process.JTAJPTCorrectorIC5.ElectronIds = cms.InputTag("eidTight")
process.JTAJPTCorrectorIC5.JetTracksAssociationAtVertex = cms.InputTag("")
process.JTAJPTCorrectorIC5.JetTracksAssociationAtCaloFace = cms.InputTag("")
#process.JTAJPTCorrectorIC5.Verbose = cms.bool(True)

process.JTAJPTCorrections = cms.Sequence(
    process.RawPATJets *
    process.PATZSPJetCorJetIcone5 *
    process.JTAJPTCorrectionIC5
    )

# ---------- End ----------

process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.InputTagDistributorService = cms.Service("InputTagDistributorService")
process.VariableHelperService = cms.Service("VariableHelperService")
process.UpdaterService = cms.Service("UpdaterService")

kine = cms.PSet(
    energy = cms.string('energy'),
    mass   = cms.string('mass'),
    mt     = cms.string('mt'),
    et     = cms.string('et'),
    pt     = cms.string('pt'),
    p      = cms.string('p'),
    px     = cms.string('px'),
    py     = cms.string('py'),
    pz     = cms.string('pz'),
    eta    = cms.string('eta'),
    phi    = cms.string('phi'),
    theta  = cms.string('theta'),
    )

process.anal = cms.EDFilter(
    "NTuplingDevice",
    Ntupler = cms.PSet(
    ComponentName = cms.string('StringBasedNTupler'),
    useTFileService = cms.bool(True), 
    branchesPSet = cms.PSet(
    treeName = cms.string('event'),

    # reco::GenJets
    RecoGenJet = cms.PSet(
    Class = cms.string('reco::GenJet'),
    src = cms.InputTag('sort:cleanLayer1Jets'),
    leaves = cms.PSet(kine),
    ),

    # Uncorrected reco::CaloJets (IC5)
    RecoCaloJetIC5 = cms.PSet(
    Class = cms.string('reco::CaloJet'),
    src = cms.InputTag('sort:iterativeCone5CaloJets'),
    leaves = cms.PSet(kine),
    ),

    # ZSP-corrected reco::CaloJets
    RecoCaloJetZSP = cms.PSet(
    Class = cms.string('reco::CaloJet'),
    src = cms.InputTag('sort:ZSPJetCorJetIcone5'),
    leaves = cms.PSet(kine),
    ),

    # JPT-corrected reco::CaloJets (old JPT)
    RecoCaloJetJPTOld = cms.PSet(
    Class = cms.string('reco::CaloJet'),
    src = cms.InputTag('sort:JetPlusTrackZSPCorJetIcone5'),
    leaves = cms.PSet(kine),
    ),

    # JPT-corrected reco::CaloJets (new JPT, all corrections, scalar)
    RecoCaloJetJPT = cms.PSet(
    Class = cms.string('reco::CaloJet'),
    src = cms.InputTag('sort:JPTCorrectorIC5'),
    leaves = cms.PSet(kine),
    ),

    # JPT-corrected reco::CaloJets (no corrections)
    RecoCaloJetJPTNone = cms.PSet(
    Class = cms.string('reco::CaloJet'),
    src = cms.InputTag('sort:JPTCorrectorNone'),
    leaves = cms.PSet(kine),
    ),

    # JPT-corrected reco::CaloJets (in-cone only)
    RecoCaloJetJPTInCone = cms.PSet(
    Class = cms.string('reco::CaloJet'),
    src = cms.InputTag('sort:JPTCorrectorInCone'),
    leaves = cms.PSet(kine),
    ),

    # JPT-corrected reco::CaloJets (+ out-of-cone)
    RecoCaloJetJPTOutOfCone = cms.PSet(
    Class = cms.string('reco::CaloJet'),
    src = cms.InputTag('sort:JPTCorrectorOutOfCone'),
    leaves = cms.PSet(kine),
    ),

    # JPT-corrected reco::CaloJets (+ out-of-vertex)
    RecoCaloJetJPTOutOfVertex = cms.PSet(
    Class = cms.string('reco::CaloJet'),
    src = cms.InputTag('sort:JPTCorrectorOutOfVertex'),
    leaves = cms.PSet(kine),
    ),

    # JPT-corrected reco::CaloJets (+ pion efficiency)
    RecoCaloJetJPTPionEff = cms.PSet(
    Class = cms.string('reco::CaloJet'),
    src = cms.InputTag('sort:JPTCorrectorPionEff'),
    leaves = cms.PSet(kine),
    ),

    # JPT-corrected reco::CaloJets (+ muons)
    RecoCaloJetJPTMuons = cms.PSet(
    Class = cms.string('reco::CaloJet'),
    src = cms.InputTag('sort:JPTCorrectorMuons'),
    leaves = cms.PSet(kine),
    ),

    # JPT-corrected reco::CaloJets (+ electrons)
    RecoCaloJetJPTElectrons = cms.PSet(
    Class = cms.string('reco::CaloJet'),
    src = cms.InputTag('sort:JPTCorrectorElectrons'),
    leaves = cms.PSet(kine),
    ),

    # JPT-corrected reco::CaloJets (+ vectorial)
    RecoCaloJetJPTVectorial = cms.PSet(
    Class = cms.string('reco::CaloJet'),
    src = cms.InputTag('sort:JPTCorrectorVectorial'),
    leaves = cms.PSet(kine),
    ),

    # JPT-corrected reco::CaloJets (+ vectorial using tracks only)
    RecoCaloJetJPTVecTracks = cms.PSet(
    Class = cms.string('reco::CaloJet'),
    src = cms.InputTag('sort:JPTCorrectorVecTracks'),
    leaves = cms.PSet(kine),
    ),
    
    # Uncorrected pat::Jets
    PatJetIC5 = cms.PSet(
    Class = cms.string('pat::Jet'),
    src = cms.InputTag('sort:RawPATJets'),
    leaves = cms.PSet(kine),
    ),

    # ZSP-corrected pat::Jets
    PatJetZSP = cms.PSet(
    Class = cms.string('pat::Jet'),
    src = cms.InputTag('sort:PATZSPJetCorJetIcone5'),
    leaves = cms.PSet(kine),
    ),

    # JPT-corrected pat::Jets (old JPT)
    PatJetJPTOld = cms.PSet(
    Class = cms.string('pat::Jet'),
    src = cms.InputTag('sort:PATJetPlusTrackZSPCorJetIcone5'),
    leaves = cms.PSet(kine),
    ),

    # JPT-corrected pat::Jets (new JPT)
    PatJetJPT = cms.PSet(
    Class = cms.string('pat::Jet'),
    src = cms.InputTag('sort:PATJPTCorrectionIC5'),
    leaves = cms.PSet(kine),
    ),

    # JPT-corrected pat::Jets (JTA on-the-fly)
    PatJetJTA = cms.PSet(
    Class = cms.string('pat::Jet'),
    src = cms.InputTag('sort:JTAJPTCorrectionIC5'),
    leaves = cms.PSet(kine),
    ),
    
    # MC-corrected pat::Jets
    PatJetMC = cms.PSet(
    Class = cms.string('pat::Jet'),
    src = cms.InputTag('sort:cleanLayer1Jets'),
    leaves = cms.PSet(kine),
    ),

    ),
    ),
    )
    
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string('ntuple.root'),
    )

process.sort = cms.EDProducer(
    "SortJetCollectionsByGenJetPt",
    genJets  = cms.InputTag("cleanLayer1Jets"),
    caloJets = cms.VInputTag(
    "iterativeCone5CaloJets",
    "ZSPJetCorJetIcone5",
    "JetPlusTrackZSPCorJetIcone5",
    "JPTCorrectorIC5",
    "JPTCorrectorNone",
    "JPTCorrectorInCone",
    "JPTCorrectorOutOfCone",
    "JPTCorrectorOutOfVertex",
    "JPTCorrectorPionEff",
    "JPTCorrectorMuons",
    "JPTCorrectorElectrons",
    "JPTCorrectorVectorial",
    "JPTCorrectorVecTracks",
    ),
    patJets = cms.VInputTag(
    "RawPATJets",
    "PATZSPJetCorJetIcone5",
    "PATJetPlusTrackZSPCorJetIcone5",
    "PATJPTCorrectionIC5",
    "JTAJPTCorrectionIC5",
    "cleanLayer1Jets",
    ),
    )

process.p = cms.Path(
    # Old JPT from RECO
    process.ZSPJetCorJetIcone5 *
    process.ZSPiterativeCone5JetTracksAssociatorAtVertex*
    process.ZSPiterativeCone5JetTracksAssociatorAtCaloFace*
    process.ZSPiterativeCone5JetExtender *
    process.JetPlusTrackZSPCorJetIcone5 *
    # New JPT from RECO
    process.JetTrackAssociations *
    process.JPTCorrectorIC5 *
    # New JPT from RECO, all steps
    process.JPTValidation *
    # PAT sequence
    process.patDefaultSequence *
    # Old JPT from PAT
    process.PATJetPlusTracksCorrections *
    # New JPT from PAT
    process.PATJPTCorrections *
    # New JTA with JTA on-the-fly
    process.JTAJPTCorrections 
    # Ntuplizer
    * process.sort 
    * process.anal
    )

process.o = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    outputCommands = cms.untracked.vstring(
    'drop *',
    'keep recoGenJets_iterativeCone5GenJets_*_*',
    'keep recoGenJets_iterativeCone5GenJetsNoNuBSM_*_*',
    'keep recoCaloJets_iterativeCone5CaloJets_*_*',
    'keep patJets_cleanLayer1Jets_*_*',
    'keep recoCaloJets_*_*_TEST',
    'keep patJets_*_*_TEST',
    'keep *_sort_*_TEST',
    )
    )

process.e = cms.EndPath( process.o )

process.ProfilerService = cms.Service(
    "ProfilerService",
    firstEvent = cms.untracked.int32(2), # avoid first event
    lastEvent = cms.untracked.int32(12),
    paths = cms.untracked.vstring("p")
    )

process.MessageLogger = cms.Service(
    "MessageLogger",
    
    debug = cms.untracked.PSet(
    threshold = cms.untracked.string('DEBUG'),
    limit = cms.untracked.uint32(100000),
    noLineBreaks = cms.untracked.bool(False),
    ),

    info = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'),
    limit = cms.untracked.uint32(100000),
    noLineBreaks = cms.untracked.bool(False),
    ),

    warning = cms.untracked.PSet(
    threshold = cms.untracked.string('WARNING'),
    limit = cms.untracked.uint32(100000),
    noLineBreaks = cms.untracked.bool(False),
    ),

    error = cms.untracked.PSet(
    threshold = cms.untracked.string('ERROR'),
    limit = cms.untracked.uint32(100000),
    noLineBreaks = cms.untracked.bool(False),
    ),
    
    cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('ERROR'),
    limit = cms.untracked.uint32(100000),
    noLineBreaks = cms.untracked.bool(False),
    ),

    destinations = cms.untracked.vstring(
    'debug', 
    'info', 
    'warning', 
    'error',
    'cerr'
    ),
    
    #@@ comment to suppress debug statements!
    debugModules = cms.untracked.vstring('*'),
    
    # allows to suppress output from specific modules 
    suppressDebug = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    
    )



