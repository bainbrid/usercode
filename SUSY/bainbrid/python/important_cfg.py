#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Includes

import setupSUSY
from hadronic.golden_cff import *
from libbainbrid import *
from ratio_cff import *
from SweetBatchSub import *

from montecarlo.LMx import *
from montecarlo.WJets_Madgraph import *
from montecarlo.TTBarTauola import *
from montecarlo.Zjets_madgraph import *
from montecarlo.Zinvisible_jets_pset import *
from montecarlo.QCD_Pythia6_allBins import *
from data.JetMETTau_ALL_230810 import *
from data.JetMET_Run2010A_PromptReco_v4_250810 import *

#import glob
#QCD_AllPtBins_7TeV_Pythia.File = glob.glob("/vols/cms02/gouskos/QCD_Pythia_Pt15_Jun2010/SusyCAF_Tree_*.root")[:1]

qcd=PSet(
      Name="QCD",
      Format=("ICF",2),
      File=[
      "dcap://gfe02.grid.hep.ph.ic.ac.uk:22128//pnfs/hep.ph.ic.ac.uk/data/cms/store/user/bm409//ICF/automated/2010_07_20_16_52_06///SusyCAF_Tree_10_1.root"
      ],
      Weight=1.,
      )

      #qcd,
      #,

data = [JetMET_ALL_upto230810,JetMET_Run2010A_PromptReco_v4_250810]

#MC = [LM0,LM1,wjets_madgraph,zjets_madgraph,Zinvisible_jets,ttbarTauola]
MC = [QCD_AllPtBins_7TeV_Pythia]

samples = []

if (True) :
      for kk in range(0,len(data)) :
            data[kk].LastEntry = 1000
      for kk in range(0,len(MC)) :
            MC[kk].LastEntry = 1000

NoiseFilt= OP_HadronicHBHEnoiseFilter()
selection = OP_GoodEventSelection()
JetTrigger = OP_TriggerCut("HLT_Jet50U")
spikecleaner = OP_EcalSpikeCleaner()
LeadingJetEta = OP_FirstJetEta(2.5)
unCorLeadJetCut = OP_UnCorLeadJetCut(75.)
secondJetET = OP_SecondJetEtCut(60.)
oddMuon = OP_OddMuon()
oddElectron = OP_OddElectron()
oddPhoton = OP_OddPhoton()
oddJet = OP_OddJet()
badMuonINJet = OP_BadMuonInJet()
numComLeptons = OP_NumComLeptons("<=",0)
numComPhotons = OP_NumComPhotons("<=",0)
numComJets = OP_NumComJets(">=",2)

ratioMC = []
createRatioPlots(ratioMC,True)
            
ratioData = []
createRatioPlots(ratioData,False)

cutTreeData = Tree("Data")
cutTreeData.Attach(spikecleaner)
cutTreeData.TAttach(spikecleaner,NoiseFilt)
cutTreeData.TAttach(NoiseFilt,selection)
cutTreeData.TAttach(selection,JetTrigger)
cutTreeData.TAttach(JetTrigger,oddMuon)
cutTreeData.TAttach(oddMuon,oddElectron)
cutTreeData.TAttach(oddElectron,oddPhoton)
cutTreeData.TAttach(oddPhoton,oddJet)
cutTreeData.TAttach(oddJet,badMuonINJet)
cutTreeData.TAttach(badMuonINJet,numComLeptons)
cutTreeData.TAttach(numComLeptons,numComPhotons)
cutTreeData.TAttach(numComPhotons,LeadingJetEta)
cutTreeData.TAttach(LeadingJetEta,unCorLeadJetCut)
cutTreeData.TAttach(unCorLeadJetCut,secondJetET)
cutTreeData.TAttach(secondJetET,numComJets)
attachRatioPlots(ratioData,cutTreeData,numComJets)
            
cutTreeMC = Tree("MC")
cutTreeMC.Attach(selection)
cutTreeMC.TAttach(selection,oddMuon)
cutTreeMC.TAttach(oddMuon,oddElectron)
cutTreeMC.TAttach(oddElectron,oddPhoton)
cutTreeMC.TAttach(oddPhoton,oddJet)
cutTreeMC.TAttach(oddJet,badMuonINJet)
cutTreeMC.TAttach(badMuonINJet,numComLeptons)
cutTreeMC.TAttach(numComLeptons,numComPhotons)
cutTreeMC.TAttach(numComPhotons,LeadingJetEta)
cutTreeMC.TAttach(LeadingJetEta,unCorLeadJetCut)
cutTreeMC.TAttach(unCorLeadJetCut,secondJetET)
cutTreeMC.TAttach(secondJetET,numComJets)
attachRatioPlots(ratioMC,cutTreeMC,numComJets)

def addCutFlowData(a) :
      a.AddJetFilter("PreCC",JetCorrections)
      a+=cutTreeData

def addCutFlowMC(b) :
      b+=cutTreeMC

# -----------------------------------------------------------------------------
# Analyses

conf_ak5_caloData = deepcopy(defaultConfig)
conf_ak5_caloData.Ntuple = deepcopy(ak5_calo)
conf_ak5_caloData.XCleaning = deepcopy(default_cc)
conf_ak5_caloData.Common = deepcopy(default_common)
conf_ak5_caloData.Common.print_out()
anal_ak5_caloData=Analysis("Ratio")
addCutFlowData(anal_ak5_caloData)
anal_ak5_caloData.Run("results",conf_ak5_caloData,data)

conf_ak5_caloMC = deepcopy(defaultConfig)
conf_ak5_caloMC.Ntuple = deepcopy(ak5_calo)
conf_ak5_caloMC.XCleaning = deepcopy(default_cc)
conf_ak5_caloMC.Common = deepcopy(default_common)
conf_ak5_caloMC.Common.print_out()
anal_ak5_caloMC=Analysis("Ratio")
addCutFlowMC(anal_ak5_caloMC)
anal_ak5_caloMC.Run("results",conf_ak5_caloMC,MC)

# Non-QCD
#BatchRun(MC,"anal_ak5_caloMC", "conf_ak5_caloMC","/vols/cms02/bainbrid/qcd/SUSY2/bainbrid/results/batch/",10000) 

# QCD
#BatchRun(MC,"anal_ak5_caloMC", "conf_ak5_caloMC","/vols/cms02/bainbrid/qcd/SUSY2/bainbrid/results/batch/",1) 

# Data
#BatchRun(data, "anal_ak5_caloData", "conf_ak5_caloData","/vols/cms02/bainbrid/qcd/SUSY2/bainbrid/results/batch/",100) 


