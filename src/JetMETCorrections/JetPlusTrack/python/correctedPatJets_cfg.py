import FWCore.ParameterSet.Config as cms

from JetMETCorrections.JetPlusTrack.correctedJets_cff import *
process.load("JetMETCorrections.JetPlusTrack.correctedPatJets_cff")

process.output.fileName = 'CorrectedPatJets.root'
process.output.outputCommands += process.keepCorrectedPatJets.outputCommands

process.p = cms.Path(
    process.correctedJetSeq *
    process.correctedPatJetSeq
    )

process.e = cms.EndPath(
    process.output +
    process.correctedPatJetHistos 
    )
