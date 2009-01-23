import FWCore.ParameterSet.Config as cms

from JetMETCorrections.JetPlusTrack.correctedJets_cff import *
process.load("JetMETCorrections.JetPlusTrack.correctedCaloJets_cff")
process.load("JetMETCorrections.JetPlusTrack.correctedPatJets_cff")

process.output.outputCommands += process.keepCorrectedCaloJets.outputCommands
process.output.outputCommands += process.keepCorrectedPatJets.outputCommands

process.p += process.correctedCaloJetSeq 
process.p += process.correctedPatJetSeq 

process.e += process.correctedCaloJetHistos 
process.e += process.correctedPatJetHistos 
