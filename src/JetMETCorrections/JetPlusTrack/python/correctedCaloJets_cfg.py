import FWCore.ParameterSet.Config as cms

from JetMETCorrections.JetPlusTrack.correctedJets_cff import *
process.load("JetMETCorrections.JetPlusTrack.correctedCaloJets_cff")

process.output.fileName = 'CorrectedCaloJets.root'
process.output.outputCommands += process.keepCorrectedCaloJets.outputCommands

process.p = cms.Path(
    process.correctedJetSeq *
    process.correctedCaloJetSeq *
    process.correctedCaloJetHistos 
    )

process.e = cms.EndPath(
    process.output 
    )
