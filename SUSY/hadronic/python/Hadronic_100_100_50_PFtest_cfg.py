#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Includes

from hadronic.golden_cff import *
from icf.core import PSet
#from SweetBatchSub import *
from plottingStuff_cfg import *
from validation.cutCounters_cff import * #TW for convenience...
from brynPlots_cff import *
# -----------------------------------------------------------------------------
# Local samples

from copy import deepcopy

from hadronic.samples_cff import *

# -----------------------------------------------------------------------------
# Cut flow


#Monte-Carlo samples
cutTreeMC = Tree("MC")
cutTreeMC.Attach(count_total)

cutTreeMC.TAttach(count_total,numComJets2OrMore)
cutTreeMC.TAttach(numComJets2OrMore,count_numComJets2om)
cutTreeMC.TAttach(numComJets2OrMore,mcTrigger2om)
cutTreeMC.TAttach(mcTrigger2om,count_Trigger2om)
cutTreeMC.TAttach(mcTrigger2om,selection2om)
cutTreeMC.TAttach(selection2om,count_selection2om)
#
cutTreeMC.TAttach(selection2om,htCutBaseM250_2om)
cutTreeMC.TAttach(htCutBaseM250_2om,count_htCutBase250_2om)
cutTreeMC.TAttach(htCutBaseM250_2om,oddPhoton2om)
cutTreeMC.TAttach(oddPhoton2om,count_oddPhoton2om)
cutTreeMC.TAttach(oddPhoton2om,numComPhotons2om)
cutTreeMC.TAttach(numComPhotons2om,count_numComPhotons2om)
cutTreeMC.TAttach(numComPhotons2om,oddElectron2om)
cutTreeMC.TAttach(oddElectron2om,count_oddElectron2om)
cutTreeMC.TAttach(oddElectron2om,numComElectrons2om)
cutTreeMC.TAttach(numComElectrons2om,count_numComElectrons2om)
cutTreeMC.TAttach(numComElectrons2om,oddMuon2om)
cutTreeMC.TAttach(oddMuon2om,count_oddMuon2om)
cutTreeMC.TAttach(oddMuon2om,numComMuons2om)
cutTreeMC.TAttach(numComMuons2om,count_numComMuons2om)
cutTreeMC.TAttach(numComMuons2om,badMuonInJet2om)
cutTreeMC.TAttach(badMuonInJet2om,count_badMuonInJet2om)
cutTreeMC.TAttach(badMuonInJet2om,oddJet2om)
cutTreeMC.TAttach(oddJet2om,count_oddJet2om)
#
cutTreeMC.TAttach(oddJet2om,htCutM250_2om) #should be a dummy cut now...
cutTreeMC.TAttach(htCutM250_2om,count_htCut250_2om)
cutTreeMC.TAttach(htCutM250_2om,plotGtHT250_2om) # BRYN PLOT
cutTreeMC.TAttach(htCutM250_2om,LeadingJetEta2om)
cutTreeMC.TAttach(LeadingJetEta2om,count_LeadingJetEta2om)
cutTreeMC.TAttach(LeadingJetEta2om,secondJetET2om)
cutTreeMC.TAttach(secondJetET2om,count_secondJetET2om)
#
#
cutTreeMC.TAttach(secondJetET2om,htCutM350_2om)
cutTreeMC.TAttach(htCutM350_2om,count_HT350plus_2om)
cutTreeMC.TAttach(htCutM350_2om,alphaT2om_350plus)
cutTreeMC.TAttach(alphaT2om_350plus,count_alphaT_350plus_2om)
cutTreeMC.TAttach(alphaT2om_350plus,deadECALmc2om_350plus)
cutTreeMC.TAttach(deadECALmc2om_350plus,count_deadecal_350plus_2om)
cutTreeMC.TAttach(deadECALmc2om_350plus,plotGtHT350_deadECAL_2om) # BRYN PLOT
cutTreeMC.TAttach(deadECALmc2om_350plus,MHToverPFMET2om_350plus)
cutTreeMC.TAttach(MHToverPFMET2om_350plus,count_MHToverPFMET_350plus_2om)
cutTreeMC.TAttach(MHToverPFMET2om_350plus,plotGtHT350_AfterAllCuts_2om) # BRYN PLOT
#
cutTreeMC.TAttach(secondJetET2om,htCutM300_2om)     #--*- ensures 300 < HT < 350
cutTreeMC.TAttach(htCutM300_2om,htCutL350_2om)      # /
cutTreeMC.TAttach(htCutL350_2om,count_HT300_350_2om)#/
cutTreeMC.TAttach(htCutL350_2om,alphaT2om_300_350)
cutTreeMC.TAttach(alphaT2om_300_350,count_alphaT_300_350_2om)
cutTreeMC.TAttach(alphaT2om_300_350,deadECALmc2om_300_350)
cutTreeMC.TAttach(deadECALmc2om_300_350,count_deadecal_300_350_2om)
cutTreeMC.TAttach(deadECALmc2om_300_350,plotHT300_350_deadECAL_2om) # BRYN PLOT
cutTreeMC.TAttach(deadECALmc2om_300_350,MHToverPFMET2om_300_350)
cutTreeMC.TAttach(MHToverPFMET2om_300_350,count_MHToverPFMET_300_350_2om)
#
cutTreeMC.TAttach(secondJetET2om,htCutL300_2om)     #-*- ensures 250 < HT < 300
cutTreeMC.TAttach(htCutL300_2om,count_HT250_300_2om)#/
cutTreeMC.TAttach(htCutL300_2om,alphaT2om_250_300)
cutTreeMC.TAttach(alphaT2om_250_300,count_alphaT_250_300_2om)
cutTreeMC.TAttach(alphaT2om_250_300,deadECALmc2om_250_300)
cutTreeMC.TAttach(deadECALmc2om_250_300,count_deadecal_250_300_2om)
cutTreeMC.TAttach(deadECALmc2om_250_300,plotHT250_300_deadECAL_2om) # BRYN PLOT
cutTreeMC.TAttach(deadECALmc2om_250_300,MHToverPFMET2om_250_300)
cutTreeMC.TAttach(MHToverPFMET2om_250_300,count_MHToverPFMET_250_300_2om)
#
#
# number of common jets = 2
cutTreeMC.TAttach(count_total,numComJets2Please)
cutTreeMC.TAttach(numComJets2Please,count_numComJets2eq)
cutTreeMC.TAttach(numComJets2Please,mcTrigger2eq)
cutTreeMC.TAttach(mcTrigger2eq,count_Trigger2eq)
cutTreeMC.TAttach(mcTrigger2eq,selection2eq)
cutTreeMC.TAttach(selection2eq,count_selection2eq)
#
cutTreeMC.TAttach(selection2eq,htCutBaseM250_2eq)
cutTreeMC.TAttach(htCutBaseM250_2eq,count_htCutBase250_2eq)
cutTreeMC.TAttach(htCutBaseM250_2eq,oddPhoton2eq)
cutTreeMC.TAttach(oddPhoton2eq,count_oddPhoton2eq)
cutTreeMC.TAttach(oddPhoton2eq,numComPhotons2eq)
cutTreeMC.TAttach(numComPhotons2eq,count_numComPhotons2eq)
cutTreeMC.TAttach(numComPhotons2eq,oddElectron2eq)
cutTreeMC.TAttach(oddElectron2eq,count_oddElectron2eq)
cutTreeMC.TAttach(oddElectron2eq,numComElectrons2eq)
cutTreeMC.TAttach(numComElectrons2eq,count_numComElectrons2eq)
cutTreeMC.TAttach(numComElectrons2eq,oddMuon2eq)
cutTreeMC.TAttach(oddMuon2eq,count_oddMuon2eq)
cutTreeMC.TAttach(oddMuon2eq,numComMuons2eq)
cutTreeMC.TAttach(numComMuons2eq,count_numComMuons2eq)
cutTreeMC.TAttach(numComMuons2eq,badMuonInJet2eq)
cutTreeMC.TAttach(badMuonInJet2eq,count_badMuonInJet2eq)
cutTreeMC.TAttach(badMuonInJet2eq,oddJet2eq)
cutTreeMC.TAttach(oddJet2eq,count_oddJet2eq)
#
cutTreeMC.TAttach(oddJet2eq,htCutM250_2eq)
cutTreeMC.TAttach(htCutM250_2eq,count_htCut250_2eq)
cutTreeMC.TAttach(htCutM250_2eq,plotGtHT250_2eq) # BRYN PLOT
cutTreeMC.TAttach(htCutM250_2eq,LeadingJetEta2eq)
cutTreeMC.TAttach(LeadingJetEta2eq,count_LeadingJetEta2eq)
cutTreeMC.TAttach(LeadingJetEta2eq,secondJetET2eq)
cutTreeMC.TAttach(secondJetET2eq,count_secondJetET2eq)
#
#
cutTreeMC.TAttach(secondJetET2eq,htCutM350_2eq)
cutTreeMC.TAttach(htCutM350_2eq,count_HT350plus_2eq)
cutTreeMC.TAttach(htCutM350_2eq,alphaT2eq_350plus)
cutTreeMC.TAttach(alphaT2eq_350plus,count_alphaT_350plus_2eq)
cutTreeMC.TAttach(alphaT2eq_350plus,deadECALmc2eq_350plus)
cutTreeMC.TAttach(deadECALmc2eq_350plus,count_deadecal_350plus_2eq)
cutTreeMC.TAttach(deadECALmc2eq_350plus,plotGtHT350_deadECAL_2eq) # BRYN PLOT
cutTreeMC.TAttach(deadECALmc2eq_350plus,MHToverPFMET2eq_350plus)
cutTreeMC.TAttach(MHToverPFMET2eq_350plus,count_MHToverPFMET_350plus_2eq)
cutTreeMC.TAttach(MHToverPFMET2eq_350plus,plotGtHT350_AfterAllCuts_2eq) # BRYN PLOT
#
cutTreeMC.TAttach(secondJetET2eq,htCutM300_2eq)     #--*- ensures 300 < HT < 350
cutTreeMC.TAttach(htCutM300_2eq,htCutL350_2eq)      # /
cutTreeMC.TAttach(htCutL350_2eq,count_HT300_350_2eq)#/
cutTreeMC.TAttach(htCutL350_2eq,alphaT2eq_300_350)
cutTreeMC.TAttach(alphaT2eq_300_350,count_alphaT_300_350_2eq)
cutTreeMC.TAttach(alphaT2eq_300_350,deadECALmc2eq_300_350)
cutTreeMC.TAttach(deadECALmc2eq_300_350,count_deadecal_300_350_2eq)
cutTreeMC.TAttach(deadECALmc2eq_300_350,plotHT300_350_deadECAL_2eq) # BRYN PLOT
cutTreeMC.TAttach(deadECALmc2eq_300_350,MHToverPFMET2eq_300_350)
cutTreeMC.TAttach(MHToverPFMET2eq_300_350,count_MHToverPFMET_300_350_2eq)
#
cutTreeMC.TAttach(secondJetET2eq,htCutL300_2eq)     #-*- ensures 250 < HT < 300
cutTreeMC.TAttach(htCutL300_2eq,count_HT250_300_2eq)#/
cutTreeMC.TAttach(htCutL300_2eq,alphaT2eq_250_300)
cutTreeMC.TAttach(alphaT2eq_250_300,count_alphaT_250_300_2eq)
cutTreeMC.TAttach(alphaT2eq_250_300,deadECALmc2eq_250_300)
cutTreeMC.TAttach(deadECALmc2eq_250_300,count_deadecal_250_300_2eq)
cutTreeMC.TAttach(deadECALmc2eq_250_300,plotHT250_300_deadECAL_2eq) # BRYN PLOT
cutTreeMC.TAttach(deadECALmc2eq_250_300,MHToverPFMET2eq_250_300)
cutTreeMC.TAttach(MHToverPFMET2eq_250_300,count_MHToverPFMET_250_300_2eq)
#
#
#number of common jets >= 3
cutTreeMC.TAttach(count_total,numComJets3OrMore)
cutTreeMC.TAttach(numComJets3OrMore,count_numComJets3om)
cutTreeMC.TAttach(numComJets3OrMore,mcTrigger3om)
cutTreeMC.TAttach(mcTrigger3om,count_Trigger3om)
cutTreeMC.TAttach(mcTrigger3om,selection3om)
cutTreeMC.TAttach(selection3om,count_selection3om)
#
cutTreeMC.TAttach(selection3om,htCutBaseM250_3om)
cutTreeMC.TAttach(htCutBaseM250_3om,count_htCutBase250_3om)
cutTreeMC.TAttach(htCutBaseM250_3om,oddPhoton3om)
cutTreeMC.TAttach(oddPhoton3om,count_oddPhoton3om)
cutTreeMC.TAttach(oddPhoton3om,numComPhotons3om)
cutTreeMC.TAttach(numComPhotons3om,count_numComPhotons3om)
cutTreeMC.TAttach(numComPhotons3om,oddElectron3om)
cutTreeMC.TAttach(oddElectron3om,count_oddElectron3om)
cutTreeMC.TAttach(oddElectron3om,numComElectrons3om)
cutTreeMC.TAttach(numComElectrons3om,count_numComElectrons3om)
cutTreeMC.TAttach(numComElectrons3om,oddMuon3om)
cutTreeMC.TAttach(oddMuon3om,count_oddMuon3om)
cutTreeMC.TAttach(oddMuon3om,numComMuons3om)
cutTreeMC.TAttach(numComMuons3om,count_numComMuons3om)
cutTreeMC.TAttach(numComMuons3om,badMuonInJet3om)
cutTreeMC.TAttach(badMuonInJet3om,count_badMuonInJet3om)
cutTreeMC.TAttach(badMuonInJet3om,oddJet3om)
cutTreeMC.TAttach(oddJet3om,count_oddJet3om)
#
cutTreeMC.TAttach(oddJet3om,htCutM250_3om)
cutTreeMC.TAttach(htCutM250_3om,count_htCut250_3om)
cutTreeMC.TAttach(htCutM250_3om,plotGtHT250_3om) # BRYN PLOT
cutTreeMC.TAttach(htCutM250_3om,LeadingJetEta3om)
cutTreeMC.TAttach(LeadingJetEta3om,count_LeadingJetEta3om)
cutTreeMC.TAttach(LeadingJetEta3om,secondJetET3om)
cutTreeMC.TAttach(secondJetET3om,count_secondJetET3om)
#
#
cutTreeMC.TAttach(secondJetET3om,htCutM350_3om)
cutTreeMC.TAttach(htCutM350_3om,count_HT350plus_3om)
cutTreeMC.TAttach(htCutM350_3om,alphaT3om_350plus)
cutTreeMC.TAttach(alphaT3om_350plus,count_alphaT_350plus_3om)
cutTreeMC.TAttach(alphaT3om_350plus,deadECALmc3om_350plus)
cutTreeMC.TAttach(deadECALmc3om_350plus,count_deadecal_350plus_3om)
cutTreeMC.TAttach(deadECALmc3om_350plus,plotGtHT350_deadECAL_3om) # BRYN PLOT
cutTreeMC.TAttach(deadECALmc3om_350plus,MHToverPFMET3om_350plus)
cutTreeMC.TAttach(MHToverPFMET3om_350plus,count_MHToverPFMET_350plus_3om)
cutTreeMC.TAttach(MHToverPFMET3om_350plus,plotGtHT350_AfterAllCuts_3om) # BRYN PLOT
#
cutTreeMC.TAttach(secondJetET3om,htCutM300_3om)     #--*- ensures 300 < HT < 350
cutTreeMC.TAttach(htCutM300_3om,htCutL350_3om)      # /
cutTreeMC.TAttach(htCutL350_3om,count_HT300_350_3om)#/
cutTreeMC.TAttach(htCutL350_3om,alphaT3om_300_350)
cutTreeMC.TAttach(alphaT3om_300_350,count_alphaT_300_350_3om)
cutTreeMC.TAttach(alphaT3om_300_350,deadECALmc3om_300_350)
cutTreeMC.TAttach(deadECALmc3om_300_350,count_deadecal_300_350_3om)
cutTreeMC.TAttach(deadECALmc3om_300_350,plotHT300_350_deadECAL_3om) # BRYN PLOT
cutTreeMC.TAttach(deadECALmc3om_300_350,MHToverPFMET3om_300_350)
cutTreeMC.TAttach(MHToverPFMET3om_300_350,count_MHToverPFMET_300_350_3om)
#
cutTreeMC.TAttach(secondJetET3om,htCutL300_3om)     #-*- ensures 250 < HT < 300
cutTreeMC.TAttach(htCutL300_3om,count_HT250_300_3om)#/
cutTreeMC.TAttach(htCutL300_3om,alphaT3om_250_300)
cutTreeMC.TAttach(alphaT3om_250_300,count_alphaT_250_300_3om)
cutTreeMC.TAttach(alphaT3om_250_300,deadECALmc3om_250_300)
cutTreeMC.TAttach(deadECALmc3om_250_300,count_deadecal_250_300_3om)
cutTreeMC.TAttach(deadECALmc3om_250_300,plotHT250_300_deadECAL_3om) # BRYN PLOT
cutTreeMC.TAttach(deadECALmc3om_250_300,MHToverPFMET3om_250_300)
cutTreeMC.TAttach(MHToverPFMET3om_250_300,count_MHToverPFMET_250_300_3om)


#-----------------------------------------------------
##data samples
cutTreeData = Tree("Data")
cutTreeData.Attach(count_total)
cutTreeData.TAttach(count_total,numComJets2OrMore)
cutTreeData.TAttach(numComJets2OrMore,count_numComJets2om)
cutTreeData.TAttach(numComJets2OrMore,dataTrigger2om)
cutTreeData.TAttach(dataTrigger2om,count_Trigger2om)
cutTreeData.TAttach(dataTrigger2om,NoiseFilt2om)
cutTreeData.TAttach(NoiseFilt2om,count_NoiseFilt2om)
cutTreeData.TAttach(NoiseFilt2om,selection2om)
cutTreeData.TAttach(selection2om,count_selection2om)
#
cutTreeData.TAttach(selection2om,htCutBaseM250_2om)
cutTreeData.TAttach(htCutBaseM250_2om,count_htCutBase250_2om)
cutTreeData.TAttach(htCutBaseM250_2om,oddPhoton2om)
cutTreeData.TAttach(oddPhoton2om,count_oddPhoton2om)
cutTreeData.TAttach(oddPhoton2om,numComPhotons2om)
cutTreeData.TAttach(numComPhotons2om,count_numComPhotons2om)
cutTreeData.TAttach(numComPhotons2om,oddElectron2om)
cutTreeData.TAttach(oddElectron2om,count_oddElectron2om)
cutTreeData.TAttach(oddElectron2om,numComElectrons2om)
cutTreeData.TAttach(numComElectrons2om,count_numComElectrons2om)
cutTreeData.TAttach(numComElectrons2om,oddMuon2om)
cutTreeData.TAttach(oddMuon2om,count_oddMuon2om)
cutTreeData.TAttach(oddMuon2om,numComMuons2om)
cutTreeData.TAttach(numComMuons2om,count_numComMuons2om)
cutTreeData.TAttach(numComMuons2om,badMuonInJet2om)
cutTreeData.TAttach(badMuonInJet2om,count_badMuonInJet2om)
cutTreeData.TAttach(badMuonInJet2om,oddJet2om)
cutTreeData.TAttach(oddJet2om,count_oddJet2om)
#
cutTreeData.TAttach(oddJet2om,htCutM250_2om)
cutTreeData.TAttach(htCutM250_2om,count_htCut250_2om)
cutTreeData.TAttach(htCutM250_2om,plotGtHT250_2om) # BRYN PLOT
cutTreeData.TAttach(htCutM250_2om,LeadingJetEta2om)
cutTreeData.TAttach(LeadingJetEta2om,count_LeadingJetEta2om)
cutTreeData.TAttach(LeadingJetEta2om,secondJetET2om)
cutTreeData.TAttach(secondJetET2om,count_secondJetET2om)
#
#
cutTreeData.TAttach(secondJetET2om,htCutM350_2om)
cutTreeData.TAttach(htCutM350_2om,count_HT350plus_2om)
cutTreeData.TAttach(htCutM350_2om,alphaT2om_350plus)
cutTreeData.TAttach(alphaT2om_350plus,count_alphaT_350plus_2om)
cutTreeData.TAttach(alphaT2om_350plus,deadECALdata2om_350plus)
cutTreeData.TAttach(deadECALdata2om_350plus,count_deadecal_350plus_2om)
cutTreeData.TAttach(deadECALdata2om_350plus,plotGtHT350_deadECAL_2om) # BRYN PLOT
cutTreeData.TAttach(deadECALdata2om_350plus,MHToverPFMET2om_350plus)
cutTreeData.TAttach(MHToverPFMET2om_350plus,count_MHToverPFMET_350plus_2om)
cutTreeData.TAttach(MHToverPFMET2om_350plus,plotGtHT350_AfterAllCuts_2om) # BRYN PLOT
cutTreeData.TAttach(MHToverPFMET2om_350plus,EventNoDump_DPhi_2om) # MARKUS EVENT DUMP
#
cutTreeData.TAttach(secondJetET2om,htCutM300_2om)     #--*- ensures 300 < HT < 350
cutTreeData.TAttach(htCutM300_2om,htCutL350_2om)      # /
cutTreeData.TAttach(htCutL350_2om,count_HT300_350_2om)#/
cutTreeData.TAttach(htCutL350_2om,alphaT2om_300_350)
cutTreeData.TAttach(alphaT2om_300_350,count_alphaT_300_350_2om)
cutTreeData.TAttach(alphaT2om_300_350,deadECALdata2om_300_350)
cutTreeData.TAttach(deadECALdata2om_300_350,count_deadecal_300_350_2om)
cutTreeData.TAttach(deadECALdata2om_300_350,plotHT300_350_deadECAL_2om) # BRYN PLOT
cutTreeData.TAttach(deadECALdata2om_300_350,MHToverPFMET2om_300_350)
cutTreeData.TAttach(MHToverPFMET2om_300_350,count_MHToverPFMET_300_350_2om)
#
cutTreeData.TAttach(secondJetET2om,htCutL300_2om)     #-*- ensures 250 < HT < 300
cutTreeData.TAttach(htCutL300_2om,count_HT250_300_2om)#/
cutTreeData.TAttach(htCutL300_2om,alphaT2om_250_300)
cutTreeData.TAttach(alphaT2om_250_300,count_alphaT_250_300_2om)
cutTreeData.TAttach(alphaT2om_250_300,deadECALdata2om_250_300)
cutTreeData.TAttach(deadECALdata2om_250_300,count_deadecal_250_300_2om)
cutTreeData.TAttach(deadECALdata2om_250_300,plotHT250_300_deadECAL_2om) # BRYN PLOT
cutTreeData.TAttach(deadECALdata2om_250_300,MHToverPFMET2om_250_300)
cutTreeData.TAttach(MHToverPFMET2om_250_300,count_MHToverPFMET_250_300_2om)
#
#
# number of common jets = 2
cutTreeData.TAttach(count_total,numComJets2Please)
cutTreeData.TAttach(numComJets2Please,count_numComJets2eq)
cutTreeData.TAttach(numComJets2Please,dataTrigger2eq)
cutTreeData.TAttach(dataTrigger2eq,count_Trigger2eq)
cutTreeData.TAttach(dataTrigger2eq,NoiseFilt2eq)
cutTreeData.TAttach(NoiseFilt2eq,count_NoiseFilt2eq)
cutTreeData.TAttach(NoiseFilt2eq,selection2eq)
cutTreeData.TAttach(selection2eq,count_selection2eq)
#
cutTreeData.TAttach(selection2eq,htCutBaseM250_2eq)
cutTreeData.TAttach(htCutBaseM250_2eq,count_htCutBase250_2eq)
cutTreeData.TAttach(htCutBaseM250_2eq,oddPhoton2eq)
cutTreeData.TAttach(oddPhoton2eq,count_oddPhoton2eq)
cutTreeData.TAttach(oddPhoton2eq,numComPhotons2eq)
cutTreeData.TAttach(numComPhotons2eq,count_numComPhotons2eq)
cutTreeData.TAttach(numComPhotons2eq,oddElectron2eq)
cutTreeData.TAttach(oddElectron2eq,count_oddElectron2eq)
cutTreeData.TAttach(oddElectron2eq,numComElectrons2eq)
cutTreeData.TAttach(numComElectrons2eq,count_numComElectrons2eq)
cutTreeData.TAttach(numComElectrons2eq,oddMuon2eq)
cutTreeData.TAttach(oddMuon2eq,count_oddMuon2eq)
cutTreeData.TAttach(oddMuon2eq,numComMuons2eq)
cutTreeData.TAttach(numComMuons2eq,count_numComMuons2eq)
cutTreeData.TAttach(numComMuons2eq,badMuonInJet2eq)
cutTreeData.TAttach(badMuonInJet2eq,count_badMuonInJet2eq)
cutTreeData.TAttach(badMuonInJet2eq,oddJet2eq)
cutTreeData.TAttach(oddJet2eq,count_oddJet2eq)
#
cutTreeData.TAttach(oddJet2eq,htCutM250_2eq)
cutTreeData.TAttach(htCutM250_2eq,plotGtHT250_2eq) # BRYN PLOT
cutTreeData.TAttach(htCutM250_2eq,count_htCut250_2eq)
cutTreeData.TAttach(htCutM250_2eq,LeadingJetEta2eq)
cutTreeData.TAttach(LeadingJetEta2eq,count_LeadingJetEta2eq)
cutTreeData.TAttach(LeadingJetEta2eq,secondJetET2eq)
cutTreeData.TAttach(secondJetET2eq,count_secondJetET2eq)
#
#
cutTreeData.TAttach(secondJetET2eq,htCutM350_2eq)
cutTreeData.TAttach(htCutM350_2eq,count_HT350plus_2eq)
cutTreeData.TAttach(htCutM350_2eq,alphaT2eq_350plus)
cutTreeData.TAttach(alphaT2eq_350plus,count_alphaT_350plus_2eq)
#
cutTreeData.TAttach(alphaT2eq_350plus,deadECALdata2eq_350plus)
cutTreeData.TAttach(deadECALdata2eq_350plus,count_deadecal_350plus_2eq)
cutTreeData.TAttach(deadECALdata2eq_350plus,plotGtHT350_deadECAL_2eq) # BRYN PLOT
cutTreeData.TAttach(deadECALdata2eq_350plus,MHToverPFMET2eq_350plus)
cutTreeData.TAttach(MHToverPFMET2eq_350plus,count_MHToverPFMET_350plus_2eq)
cutTreeData.TAttach(MHToverPFMET2eq_350plus,plotGtHT350_AfterAllCuts_2eq) # BRYN PLOT
cutTreeData.TAttach(MHToverPFMET2eq_350plus,EventNoDump_DPhi_2eq) # MARKUS EVENT DUMP
#
cutTreeData.TAttach(secondJetET2eq,htCutM300_2eq)     #--*- ensures 300 < HT < 350
cutTreeData.TAttach(htCutM300_2eq,htCutL350_2eq)      # /
cutTreeData.TAttach(htCutL350_2eq,count_HT300_350_2eq)#/
cutTreeData.TAttach(htCutL350_2eq,alphaT2eq_300_350)
cutTreeData.TAttach(alphaT2eq_300_350,count_alphaT_300_350_2eq)
cutTreeData.TAttach(alphaT2eq_300_350,deadECALdata2eq_300_350)
cutTreeData.TAttach(deadECALdata2eq_300_350,count_deadecal_300_350_2eq)
cutTreeData.TAttach(deadECALdata2eq_300_350,plotHT300_350_deadECAL_2eq) # BRYN PLOT
cutTreeData.TAttach(deadECALdata2eq_300_350,MHToverPFMET2eq_300_350)
cutTreeData.TAttach(MHToverPFMET2eq_300_350,count_MHToverPFMET_300_350_2eq)
#
cutTreeData.TAttach(secondJetET2eq,htCutL300_2eq)     #-*- ensures 250 < HT < 300
cutTreeData.TAttach(htCutL300_2eq,count_HT250_300_2eq)#/
cutTreeData.TAttach(htCutL300_2eq,alphaT2eq_250_300)
cutTreeData.TAttach(alphaT2eq_250_300,count_alphaT_250_300_2eq)
cutTreeData.TAttach(alphaT2eq_250_300,deadECALdata2eq_250_300)
cutTreeData.TAttach(deadECALdata2eq_250_300,count_deadecal_250_300_2eq)
cutTreeData.TAttach(deadECALdata2eq_250_300,plotHT250_300_deadECAL_2eq) # BRYN PLOT
cutTreeData.TAttach(deadECALdata2eq_250_300,MHToverPFMET2eq_250_300)
cutTreeData.TAttach(MHToverPFMET2eq_250_300,count_MHToverPFMET_250_300_2eq)
#
#
#number of common jets >= 3
cutTreeData.TAttach(count_total,numComJets3OrMore)
cutTreeData.TAttach(numComJets3OrMore,count_numComJets3om)
cutTreeData.TAttach(numComJets3OrMore,dataTrigger3om)
cutTreeData.TAttach(dataTrigger3om,count_Trigger3om)
cutTreeData.TAttach(dataTrigger3om,NoiseFilt3om)
cutTreeData.TAttach(NoiseFilt3om,count_NoiseFilt3om)
cutTreeData.TAttach(NoiseFilt3om,selection3om)
cutTreeData.TAttach(selection3om,count_selection3om)
#
#insert here
cutTreeData.TAttach(selection3om,htCutBaseM250_3om)
cutTreeData.TAttach(htCutBaseM250_3om,count_htCutBase250_3om)
cutTreeData.TAttach(htCutBaseM250_3om,oddPhoton3om)
cutTreeData.TAttach(oddPhoton3om,count_oddPhoton3om)
cutTreeData.TAttach(oddPhoton3om,numComPhotons3om)
cutTreeData.TAttach(numComPhotons3om,count_numComPhotons3om)
cutTreeData.TAttach(numComPhotons3om,oddElectron3om)
cutTreeData.TAttach(oddElectron3om,count_oddElectron3om)
cutTreeData.TAttach(oddElectron3om,numComElectrons3om)
cutTreeData.TAttach(numComElectrons3om,count_numComElectrons3om)
cutTreeData.TAttach(numComElectrons3om,oddMuon3om)
cutTreeData.TAttach(oddMuon3om,count_oddMuon3om)
cutTreeData.TAttach(oddMuon3om,numComMuons3om)
cutTreeData.TAttach(numComMuons3om,count_numComMuons3om)
cutTreeData.TAttach(numComMuons3om,badMuonInJet3om)
cutTreeData.TAttach(badMuonInJet3om,count_badMuonInJet3om)
cutTreeData.TAttach(badMuonInJet3om,oddJet3om)
cutTreeData.TAttach(oddJet3om,count_oddJet3om)
#
cutTreeData.TAttach(oddJet3om,htCutM250_3om)
cutTreeData.TAttach(htCutM250_3om,count_htCut250_3om)
cutTreeData.TAttach(htCutM250_3om,plotGtHT250_3om) # BRYN PLOT
cutTreeData.TAttach(htCutM250_3om,LeadingJetEta3om)
cutTreeData.TAttach(LeadingJetEta3om,count_LeadingJetEta3om)
cutTreeData.TAttach(LeadingJetEta3om,secondJetET3om)
cutTreeData.TAttach(secondJetET3om,count_secondJetET3om)
#
#
cutTreeData.TAttach(secondJetET3om,htCutM350_3om)
cutTreeData.TAttach(htCutM350_3om,count_HT350plus_3om)
cutTreeData.TAttach(htCutM350_3om,alphaT3om_350plus)
cutTreeData.TAttach(alphaT3om_350plus,count_alphaT_350plus_3om)
cutTreeData.TAttach(alphaT3om_350plus,deadECALdata3om_350plus)
#
cutTreeData.TAttach(deadECALdata3om_350plus,count_deadecal_350plus_3om)
cutTreeData.TAttach(deadECALdata3om_350plus,plotGtHT350_deadECAL_3om) # BRYN PLOT
cutTreeData.TAttach(deadECALdata3om_350plus,MHToverPFMET3om_350plus)
cutTreeData.TAttach(MHToverPFMET3om_350plus,count_MHToverPFMET_350plus_3om)
cutTreeData.TAttach(MHToverPFMET3om_350plus,plotGtHT350_AfterAllCuts_3om) # BRYN PLOT
cutTreeData.TAttach(MHToverPFMET3om_350plus,EventNoDump_DPhi_3om) # MARKUS EVENT DUMP
#
cutTreeData.TAttach(secondJetET3om,htCutM300_3om)     #--*- ensures 300 < HT < 350
cutTreeData.TAttach(htCutM300_3om,htCutL350_3om)      # /
cutTreeData.TAttach(htCutL350_3om,count_HT300_350_3om)#/
cutTreeData.TAttach(htCutL350_3om,alphaT3om_300_350)
cutTreeData.TAttach(alphaT3om_300_350,count_alphaT_300_350_3om)
cutTreeData.TAttach(alphaT3om_300_350,deadECALdata3om_300_350)
cutTreeData.TAttach(deadECALdata3om_300_350,count_deadecal_300_350_3om)
cutTreeData.TAttach(deadECALdata3om_300_350,plotHT300_350_deadECAL_3om) # BRYN PLOT
cutTreeData.TAttach(deadECALdata3om_300_350,MHToverPFMET3om_300_350)
cutTreeData.TAttach(MHToverPFMET3om_300_350,count_MHToverPFMET_300_350_3om)
#
cutTreeData.TAttach(secondJetET3om,htCutL300_3om)     #-*- ensures 250 < HT < 300
cutTreeData.TAttach(htCutL300_3om,count_HT250_300_3om)#/
cutTreeData.TAttach(htCutL300_3om,alphaT3om_250_300)
cutTreeData.TAttach(alphaT3om_250_300,count_alphaT_250_300_3om)
cutTreeData.TAttach(alphaT3om_250_300,deadECALdata3om_250_300)
cutTreeData.TAttach(deadECALdata3om_250_300,count_deadecal_250_300_3om)
cutTreeData.TAttach(deadECALdata3om_250_300,plotHT250_300_deadECAL_3om) # BRYN PLOT
cutTreeData.TAttach(deadECALdata3om_250_300,MHToverPFMET3om_250_300)
cutTreeData.TAttach(MHToverPFMET3om_250_300,count_MHToverPFMET_250_300_3om)


# Object ID filters
#-----------------------------------------------------------------------
from ra1objectid.vbtfElectronId_cff import *
from ra1objectid.vbtfMuonId_cff import *
from ra1objectid.ra3PhotonId_cff import *

vbtfElectronIdFilter = Electron_IDFilter( vbtfelectronidWP95ps.ps() )
#vbtfElectronIdFilter = Electron_IDFilter( vbtfelectronidWP90ps.ps() )
#vbtMuonIdFilter      = Muon_IDFilter( vbtfmuonidps.ps() )
ra3PhotonIdFilter    = Photon_IDFilter( ra3photonidps.ps() )

def addCutFlowMC(a) :
    a.AddPhotonFilter("PreCC",ra3PhotonIdFilter)
    a.AddElectronFilter("PreCC",vbtfElectronIdFilter)
#    a.AddMuonFilter("PreCC",vbtMuonIdFilter)
    a+=cutTreeMC

def addCutFlowData(b) :
    b.AddJetFilter("PreCC",JetCorrections)
    b.AddPhotonFilter("PreCC",ra3PhotonIdFilter)
    b.AddElectronFilter("PreCC",vbtfElectronIdFilter)
#    b.AddMuonFilter("PreCC",vbtMuonIdFilter)
    b+=cutTreeData


#    a+=debug_hadronic

# -----------------------------------------------------------------------------
# Definition of analyses


#MC
anal_ak5_caloMC=Analysis("AK5Calo")
addCutFlowMC(anal_ak5_caloMC)

anal_ak7_caloMC=Analysis("AK7Calo")
addCutFlowMC(anal_ak7_caloMC)

anal_ak5_pfMC=Analysis("AK5PF")
addCutFlowMC(anal_ak5_pfMC)

#data

anal_ak5_caloData=Analysis("AK5Calo")
addCutFlowData(anal_ak5_caloData)

anal_ak7_caloData=Analysis("AK7Calo")
addCutFlowData(anal_ak7_caloData)

anal_ak5_pfData=Analysis("AK5PF")
addCutFlowData(anal_ak5_pfData)

#config modifications
conf_ak5_pf.XCleaning.Muons.ModifyJetEnergy = False

#anal_ak5_pfData.Run("../results_100_100_50_PFtest",conf_ak5_pf,[Jet_15pb_WithTP_json221010]) #15.0865/pb
#anal_ak5_pfData.Run("../results_100_100_50_PFtest",conf_ak5_pf,[MultiJet_json291010]) #+7/pb
#anal_ak5_pfData.Run("../results_100_100_50_PFtest",conf_ak5_pf,[MultiJet_Run2010B_PromptReco_v2]) #+7/pb from 05/11/2010 json

# -----------------------------------------------------------------------------
# Run analyses

#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf_msugra,[tanB3])
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf_msugra,[tanB10])
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,[LM0])
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,[LM1])
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,[LM2])
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,[LM3])
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,[LM4])
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,[LM5])

#QCD
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,[QCD_AllPtBins_7TeV_Pythia])
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,QCD_Pythia8_ALL)
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,QCD_MadGraph_ALL)
#
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,QCD_Pythia6_384patch3_V14_00_02_ALL)
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,QCD_Pythia8_384patch3_V14_00_02_ALL)

#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,[PhotonJet_AllPtBins_7TeV_Pythia])
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,PhotonJet_MadGraph_ALL)

#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,[ttbarTauola])
#wjets_madgraph_vols.FirstEntry = 0
#wjets_madgraph_vols.LastEntry = 199999
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,[wjets_madgraph])
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,[wjets_madgraph_vols])
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,[zjets_madgraph])
#anal_ak5_caloMC.Run("../results_100_100_50_PFtest",conf_ak5_pf,[Zinvisible_jets])

#myMagical14events = PSet(#from 15pb-1
#    Name="14magicalEvents",
#    Weight=1.,
#    Format=("ICF",2),
#    File=[
#        "/vols/cms02/whyntie/public/ntuple_14eventsPassingAlphaT_15invpb.root",
#        ]
#    )
#anal_ak5_pfData.Run("../results_100_100_50_PFtest",conf_ak5_pf,[myMagical14events]) #15.0865/pb
