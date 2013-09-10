#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Includes

import setupSUSY
from libFrameworkSUSY import *
from libHadronic import *
from icf.core import PSet,Analysis
from icf.config import defaultConfig
from copy import deepcopy

# -----------------------------------------------------------------------------
# Samples

from samples_cfi import *

# two lines below should be commented!
#from local_cff import * 
#from rob_cff import *

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
default_ntuple.Taus=PSet(
    Prefix="tau",
    Suffix="Pat",
    )
default_ntuple.Jets=PSet(
    Prefix="ic5Jet",
    Suffix="Pat",
    )
default_ntuple.Photons=PSet(
    Prefix="photon",
    Suffix="Pat",
    )

ak5_calo = deepcopy(default_ntuple)
ak5_calo.Jets.Prefix="ak5Jet"

ak5_jpt = deepcopy(default_ntuple)
ak5_jpt.Jets.Prefix="ak5JetJPT"

# -----------------------------------------------------------------------------
# Cross-cleaning settings

default_cc = deepcopy(defaultConfig.XCleaning)
default_cc.Verbose=False
default_cc.MuonJet=True
default_cc.ElectronJet=True
default_cc.PhotonJet=True
default_cc.ResolveConflicts=False
default_cc.Jets.EtCut=0.
default_cc.Jets.EtaCut=1000.0
default_cc.Muons.PtCut=0.
default_cc.Muons.EtaCut=1000.
default_cc.Muons.TrkIsoCut=6.0
default_cc.Muons.CombIsoCut=0.2
default_cc.Electrons.PtCut=0.0
default_cc.Electrons.EtaCut=1000.0
default_cc.Electrons.TrkIsoCut=6.0
default_cc.Electrons.CombIsoCut=0.2
default_cc.Photons.EtCut=0.0
default_cc.Photons.EtaCut=1000.0
default_cc.Photons.TrkIsoCut=9.
default_cc.Photons.AllIsoCut=0.2
default_cc.Photons.IDReq=3

# -----------------------------------------------------------------------------
# Definition of common objects

default_common = deepcopy(defaultConfig.Common)

default_common.ApplyXCleaning=True
default_common.Jets.EtCut=50.0
default_common.Jets.EtaCut=3.0
default_common.Electrons.PtCut=10.0
default_common.Electrons.EtaCut=2.5
default_common.TrkIsoCut=6.
default_common.CombIsoCut=0.2
default_common.Muons.PtCut=10.0
default_common.Muons.EtaCut=2.5
default_common.Muons.TrkIsoCut=6.
default_common.Muons.CombIsoCut=0.2
default_common.Photons.EtCut=25.0
default_common.Photons.EtaCut=2.5
default_common.Photons.TrkIsoCut=9.
default_common.Photons.EcalIsoCut=99999.
default_common.Photons.HcalIsoCut=99999.
default_common.Photons.AllIsoCut=0.2
default_common.Photons.IDReq=3

# -----------------------------------------------------------------------------
# Cut flow and plots

numComLeptons = OP_NumComLeptons("<=",0)
numComJets = OP_NumComJets(">=",2)
numComPhotons = OP_NumComPhotons("<=",0)
oddElectron = OP_OddElectron()
oddMuon = OP_OddMuon()
oddJet = OP_OddJet()
oddPhoton = OP_OddPhoton()
badMuonInJet = OP_BadMuonInJet()
photonKilledJet = OP_PhotonKilledJet()
secondJetET = OP_SecondJetEtCut(100)
htCut = RECO_CommonHTCut(350)
missedHT = OP_MissedHTCut(1.25)
HadAlphaT = OP_HadronicAlphaT(0.55)

HadronicPlottingOps = PSet(
    DirectoryName = "Hadronic",
    MinObjects = 2,
    MaxObjects = 10,
    Dalitz = True,
    AlphaT = False,
    PtHat  = False,
    )
myplots = OP_HadronicPlottingOps( HadronicPlottingOps.ps() )

def addCutFlow(a) :
    a+=numComLeptons
    a+=numComJets
    a+=numComPhotons
    a+=oddElectron
    a+=oddMuon
    a+=oddJet
    a+=oddPhoton
    a+=badMuonInJet
    a+=photonKilledJet
    a+=secondJetET
    a+=htCut
    a+=missedHT
    a+=kinsuitecomplot
    a+=myplots
    a+=HadAlphaT
    
# -----------------------------------------------------------------------------
# Analyses

# AK5 Calo

conf_ak5_calo = deepcopy(defaultConfig)
conf_ak5_calo.Ntuple = deepcopy(ak5_calo)
conf_ak5_calo.XCleaning = deepcopy(default_cc)
conf_ak5_calo.Common = deepcopy(default_common)

anal_ak5_calo=Analysis("AK5Calo")
addCutFlow(anal_ak5_calo)

# AK5 JPT 

conf_ak5_jpt = deepcopy(defaultConfig)
conf_ak5_jpt.Ntuple = deepcopy(ak5_jpt)
conf_ak5_jpt.XCleaning = deepcopy(default_cc)
conf_ak5_jpt.Common = deepcopy(default_common)

anal_ak5_jpt=Analysis("AK5JPT")
addCutFlow(anal_ak5_jpt)

# Run analyses

#anal_ak5_calo.Run("results",conf_ak5_calo,[lm0])
#anal_ak5_calo.Run("results",conf_ak5_calo,[z_inv])

if ( True ) :
    anal_ak5_calo.Run("results",conf_ak5_calo,[lm0])
    anal_ak5_calo.Run("results",conf_ak5_calo,[lm1])
    anal_ak5_calo.Run("results",conf_ak5_calo,[qcd_pythia_merged])
    anal_ak5_calo.Run("results",conf_ak5_calo,[z_inv])
    anal_ak5_calo.Run("results",conf_ak5_calo,[w_jets])
    anal_ak5_calo.Run("results",conf_ak5_calo,[ttbar_jets])
    
