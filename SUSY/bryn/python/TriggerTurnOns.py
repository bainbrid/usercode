#!/usr/bin/env python

#!/usr/bin/env python
"""
Created by Bryn Mathias on 2010-05-07.
"""

# -----------------------------------------------------------------------------
# Necessary includes
import errno
import os
import setupSUSY
from libFrameworkSUSY import *
from libHadronic import *
from libbryn import *
from icf.core import PSet,Analysis
from icf.config import defaultConfig
from copy import deepcopy
from icf.JetCorrections import *
from icf.utils import json_to_pset

# -----------------------------------------------------------------------------
# Samples
#import yours in your running script
def ensure_dir(path):
    try:
      os.makedirs(path)
    except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST:
        pass
      else: raise

# -----------------------------------------------------------------------------
# Reading the collections from the ntuple

default_ntuple = deepcopy(defaultConfig.Ntuple)
default_ntuple.Electrons=PSet(
Prefix="electron",
Suffix="Pat",
LooseID="EIDLoose",
TightID="EIDTight",
)
default_ntuple.Muons=PSet(
Prefix="muon",
Suffix="Pat",
LooseID="IsGlobalMuon",
TightID="IDGlobalMuonPromptTight",
)
default_ntuple.SecMuons=PSet(
    Prefix="muon",
    Suffix="PF")
default_ntuple.Taus=PSet(
Prefix="tau",
Suffix="Pat",
LooseID="TauIdbyTaNCfrOnePercent",
TightID="TauIdbyTaNCfrTenthPercent"
)
default_ntuple.Jets=PSet(
Prefix="ic5Jet",
Suffix="Pat",
Uncorrected=False,
)
default_ntuple.Photons=PSet(
Prefix="photon",
Suffix="Pat",
)

ic5_calo = deepcopy(default_ntuple)
ic5_calo.Jets.Prefix="ic5Jet"

ak5_calo = deepcopy(default_ntuple)
ak5_calo.Jets.Prefix="ak5Jet"

ak5_jpt = deepcopy(default_ntuple)
ak5_jpt.Jets.Prefix="ak5JetJPT"

ak5_pf = deepcopy(default_ntuple)
ak5_pf.Jets.Prefix="ak5JetPF"
ak5_pf.TerJets.Prefix="ak5Jet"

ak7_calo = deepcopy(default_ntuple)
ak7_calo.Jets.Prefix="ak7Jet"


# -----------------------------------------------------------------------------
# Cross-cleaning settings

default_cc = deepcopy(defaultConfig.XCleaning)
default_cc.Verbose=False
default_cc.MuonJet=True
default_cc.ElectronJet=True
default_cc.PhotonJet=True
default_cc.ResolveConflicts=True
default_cc.Jets.PtCut=10.0
default_cc.Jets.EtaCut=10.0
default_cc.Muons.ModifyJetEnergy=True
default_cc.Muons.PtCut=10.0
default_cc.Muons.EtaCut=2.5
default_cc.Muons.TrkIsoCut=-1.
default_cc.Muons.CombIsoCut=0.15
default_cc.Muons.MuonJetDeltaR=0.5
default_cc.Muons.MuonIsoTypePtCutoff=0.
default_cc.Muons.RequireLooseIdForInitialFilter=False
default_cc.Electrons.PtCut=10.0
default_cc.Electrons.EtaCut=2.5
default_cc.Electrons.TrkIsoCut=-1.0
default_cc.Electrons.CombIsoCut=0.15
default_cc.Electrons.ElectronJetDeltaR=0.5
default_cc.Electrons.ElectronIsoTypePtCutoff=0.
default_cc.Electrons.ElectronLooseIdRequired=True
default_cc.Electrons.ElectronTightIdRequired=False
default_cc.Electrons.RequireLooseIdForInitialFilter=False
default_cc.Photons.EtCut=25.0
default_cc.Photons.EtaCut=2.5
default_cc.Photons.TrkIsoCut=2.0
default_cc.Photons.CaloIsoCut=0.2
default_cc.Photons.IDReq=3
default_cc.Photons.PhotonJetDeltaR=0.5
default_cc.Photons.PhotonIsoTypePtCutoff=30.
# -----------------------------------------------------------------------------
# Definition of common objects
default_common = deepcopy(defaultConfig.Common)

default_common.ApplyXCleaning=True
default_common.Jets.PtCut=50.0
default_common.Jets.EtaCut=3.0
default_common.Jets.ApplyID=True
default_common.Jets.TightID=False
default_common.Electrons.PtCut=10.0
default_common.Electrons.EtaCut=2.5
default_common.Electrons.TrkIsoCut=-1.
default_common.Electrons.CombIsoCut=0.15
default_common.Electrons.ApplyID = True
default_common.Electrons.TightID = False
default_common.Electrons.RequireLooseForOdd = True
default_common.Muons.PtCut=10.0
default_common.Muons.EtaCut=2.5
default_common.Muons.TrkIsoCut=-1.
default_common.Muons.CombIsoCut=0.15
default_common.Muons.ApplyID = True
default_common.Muons.TightID = True
default_common.Muons.RequireLooseForOdd = True
default_common.Photons.EtCut=25.0
default_common.Photons.EtaCut=2.5
# the photon cuts are NOT read anyway
# default_common.Photons.TrkIsoRel=0.
# default_common.Photons.TrkIsoCut=99999.
# default_common.Photons.EcalIsoRel=0.
# default_common.Photons.EcalIsoCut=99999.
# default_common.Photons.HcalIsoRel=0.
# default_common.Photons.HcalIsoCut=99999.
# default_common.Photons.HadOverEmCut=0.5
# default_common.Photons.SigmaIetaIetaCut=0.5
##default_common.Photons.CaloIsoCut=99999.
default_common.Photons.IDReq = 3
default_common.Photons.RequireLooseForOdd = True

# CMSSW_3_8_4_patch3 V14-00-02 samples

# -----------------------------------------------------------------------------


# Common cut definitions
#Avaiable criteria for MC and for Data are at current slightly different Hence the making of two trees
#DataOnly!

from icf.JetCorrections import *
corPset =  CorrectionPset("ResidualJetEnergyCorrections.txt")
# corPset =  CorrectionPset("Spring10DataV2_L2L3Residual_AK5PF.txt")
JetCorrections = JESCorrections( corPset.ps() ,True)
NoiseFilt= OP_HadronicHBHEnoiseFilter()
selection = OP_GoodEventSelection()


Plots_HT_Trigger    = PL_TriggerTurnOns( PSet(DirName = "HT_Trigger",MinObjects =0 ,MaxObjects = 15,Plots = True).ps())
Plots_Cross_Trigger = PL_TriggerTurnOns( PSet(DirName = "Cross_Trigger",MinObjects =0 ,MaxObjects = 15,Plots = True).ps())

HT_Trigger_PS = PSet(
    Verbose = False,
    Triggers = [
        "HLT_HT250_v1",
        "HLT_HT250_v2",
        "HLT_HT250_v3",
        "HLT_HT250_v4",
        "HLT_HT250_v5",
        "HLT_HT250_v6",
        ]
    )


Cross_Trigger_PS = PSet(
    Verbose = False,
    Triggers = [
        "HLT_HT250_AlphaT0p53_v1",
        "HLT_HT250_AlphaT0p53_v2",
        "HLT_HT250_AlphaT0p53_v3",
        "HLT_HT250_AlphaT0p53_v4",
        "HLT_HT250_AlphaT0p53_v5",
        ]
    )





HT_Trigger_Filter = OP_MultiTrigger( HT_Trigger_PS.ps() )
Cross_Trigger_Filter = OP_MultiTrigger( Cross_Trigger_PS.ps() )



#Standard Event Selection
LeadingJetEta = OP_FirstJetEta(2.5)
unCorLeadJetCut = OP_UnCorLeadJetCut(0.)
LeadingJetCut = OP_FirstJetPtCut(10.)
oddMuon = OP_OddMuon()
oddElectron = OP_OddElectron()
oddPhoton = OP_OddPhoton()
oddJet = OP_OddJet()
secondJetET = OP_SecondJetEtCut(100.)
badMuonInJet = OP_BadMuonInJet()
numComLeptons = OP_NumComLeptons("<=",0)
numComPhotons = OP_NumComPhotons("<=",0)


DeadEcalCutData = OP_DeadECALCut(0.3,0.3,0.5,30.,10,0,"./deadRegionList_GR10_P_V10.txt")
# Cross check with the allhadronic analysis
VertexPtOverHT = OP_SumVertexPtOverHT(0.1)
# -----------------------------------------------------------------------------
# Definition of analyses
# Analyses
MHT_METCut = OP_MHToverMET(1.25)
# AK5 Calo
json = JSONFilter("Json Mask", json_to_pset("goldenJSON_13May2011.txt"))

cutTreeData = Tree("Data")
cutTreeData.Attach(json)
cutTreeData.TAttach(json,NoiseFilt)
cutTreeData.TAttach(NoiseFilt,selection)
cutTreeData.TAttach(selection,oddMuon)
cutTreeData.TAttach(oddMuon,oddElectron)
cutTreeData.TAttach(oddElectron,oddPhoton)
cutTreeData.TAttach(oddPhoton,numComLeptons)
cutTreeData.TAttach(numComLeptons,numComPhotons)
cutTreeData.TAttach(numComPhotons,LeadingJetEta)
cutTreeData.TAttach(LeadingJetEta,badMuonInJet)
cutTreeData.TAttach(badMuonInJet,oddJet)
cutTreeData.TAttach(oddJet,LeadingJetCut)
cutTreeData.TAttach(LeadingJetCut,secondJetET)
##########DiJet Studies
cutTreeData.TAttach(secondJetET,VertexPtOverHT)
cutTreeData.TAttach(VertexPtOverHT,DeadEcalCutData)
cutTreeData.TAttach(DeadEcalCutData,MHT_METCut)
cutTreeData.TAttach(MHT_METCut,HT_Trigger_Filter)
cutTreeData.TAttach(HT_Trigger_Filter,Plots_HT_Trigger)
cutTreeData.TAttach(HT_Trigger_Filter,Cross_Trigger_Filter)
cutTreeData.TAttach(Cross_Trigger_Filter,Plots_Cross_Trigger)

from ra1objectid.vbtfElectronId_cff import *
from ra1objectid.vbtfMuonId_cff import *
from ra1objectid.ra3PhotonId_cff import *
vbtfMuonId_cff = Muon_IDFilter( vbtfmuonidps.ps()  )
vbtfElectronIdFilter = Electron_IDFilter( vbtfelectronidWP95ps.ps() )
ra3PhotonIdFilter    = Photon_IDFilter( ra3photonidps.ps() )


def addCutFlowData(a) :
  # a.AddWeightFilter("Weight",PreScaleWeights)
  # a.AddJetFilter("PreCC",JetCorrections)
  a.AddPhotonFilter("PreCC",ra3PhotonIdFilter)
  a.AddElectronFilter("PreCC",vbtfElectronIdFilter)
  a.AddMuonFilter("PreCC",vbtfMuonId_cff)
  a+=cutTreeData

# AK5 Calo

conf_ak5_caloData = deepcopy(defaultConfig)
conf_ak5_caloData.Ntuple = deepcopy(ak5_calo)
conf_ak5_caloData.XCleaning = deepcopy(default_cc)
conf_ak5_caloData.Common = deepcopy(default_common)
# conf_ak5_calo.Common.print_out()
anal_ak5_caloData=Analysis("AK5Calo")
addCutFlowData(anal_ak5_caloData)

# AK5 PF

conf_ak5_pfData = deepcopy(defaultConfig)
conf_ak5_pfData.Ntuple = deepcopy(ak5_pf)
conf_ak5_pfData.XCleaning = deepcopy(default_cc)
conf_ak5_pfData.Common = deepcopy(default_common)
anal_ak5_pfData=Analysis("AK5PF")
addCutFlowData(anal_ak5_pfData)


from data.Run2011.HT_Run2011_promptReco_DCS import *
outDir = "../results/TriggerTurnOns/"
ensure_dir(outDir)
anal_ak5_caloData.Run(outDir,conf_ak5_caloData,[HT_Run2011_promptReco_DCS])
# anal_ak5_pfData.Run(ourDir,conf_ak5_pfData,[HT_Run2011_promptReco_DCS])
