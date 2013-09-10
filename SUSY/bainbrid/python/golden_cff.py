#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Necessary includes

import setupSUSY
from libFrameworkSUSY import *
from libbainbrid import *
from icf.core import PSet,Analysis
from icf.config import defaultConfig
from copy import deepcopy

# -----------------------------------------------------------------------------
# Samples

from samples.samples import *

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
    LooseID="TauIdbyTaNCfrOnePercent",
    TightID="TauIdbyTaNCfrTenthPercent"
    )
default_ntuple.Jets=PSet(
    Prefix="ic5Jet",
    Suffix="Pat",
    Uncorrected=False,
    UseGenJets=False
    )
default_ntuple.Photons=PSet(
    Prefix="photon",
    Suffix="Pat",
    )

ic5_gen = deepcopy(default_ntuple)
ic5_gen.Jets.Prefix="ic5Jet"
ic5_gen.Jets.UseGenJets=True

ak5_gen = deepcopy(default_ntuple)
ak5_gen.Jets.Prefix="ak5Jet"
ak5_gen.Jets.UseGenJets=True

ic5_calo = deepcopy(default_ntuple)
ic5_calo.Jets.Prefix="ic5Jet"

ak5_calo = deepcopy(default_ntuple)
ak5_calo.Jets.Prefix="ak5Jet"

ak5_jpt = deepcopy(default_ntuple)
ak5_jpt.Jets.Prefix="ak5JetJPT"

ak5_pf = deepcopy(default_ntuple)
ak5_pf.Jets.Prefix="ak5JetPF"

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
default_cc.Jets.PtCut=20.0
default_cc.Jets.EtaCut=10.0
default_cc.Muons.PtCut=10.0
default_cc.Muons.EtaCut=2.5
default_cc.Muons.TrkIsoCut=6.0
default_cc.Muons.CombIsoCut=0.2
default_cc.Muons.ModifyJetEnergy=True
default_cc.Electrons.PtCut=10.0
default_cc.Electrons.EtaCut=2.5
default_cc.Electrons.TrkIsoCut=6.0
default_cc.Electrons.CombIsoCut=0.2
default_cc.Photons.EtCut=25.0
default_cc.Photons.EtaCut=2.5
default_cc.Photons.TrkIsoCut=2.0
default_cc.Photons.CaloIsoCut=0.2
default_cc.Photons.IDReq=2

gen_cc = deepcopy(default_cc)
gen_cc.MuonJet = False
gen_cc.ElectronJet = False
gen_cc.PhotonJet = False
gen_cc.ResolveConflicts = False

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
default_common.Electrons.TrkIsoCut=6.
default_common.Electrons.CombIsoCut=0.2
default_common.Muons.PtCut=10.0
default_common.Muons.EtaCut=2.5
default_common.Muons.TrkIsoCut=6.
default_common.Muons.CombIsoCut=0.2
default_common.Photons.EtCut=25.0
default_common.Photons.EtaCut=2.5
default_common.Photons.TrkIsoCut=99999.
default_common.Photons.TrkIsoRel=0.
default_common.Photons.EcalIsoCut=99999.
default_common.Photons.EcalIsoRel=0.
default_common.Photons.HcalIsoCut=99999.
default_common.Photons.HcalIsoRel=0.
default_common.Photons.HadOverEmCut=0.5
default_common.Photons.SigmaIetaIetaCut=0.5
#default_common.Photons.CaloIsoCut=0.2
default_common.Photons.IDReq=2
                                                               
gen_common = deepcopy(default_common)
gen_common.ApplyXCleaning=False
gen_common.Jets.ApplyID=False

# -----------------------------------------------------------------------------
# Definition of pre-selection

numComJets = OP_NumComJets(">=",2)
numComLeptons = OP_NumComLeptons("<=",0)
numComPhotons = OP_NumComPhotons("<=",0)
oddJet = OP_OddJet()
oddMuon = OP_OddMuon()
oddElectron = OP_OddElectron()
oddPhoton = OP_OddPhoton()
badMuonInJet = OP_BadMuonInJet()
secondJetET = OP_SecondJetEtCut(100.) 
htCut = RECO_CommonHTCut(350.)
missedHT = OP_MissedHTCut(1.25)
alphaT = OP_CommonAlphaTCut(0.55)

def addPreSelection(a) :
    a+=numComJets
    a+=numComLeptons
    a+=numComPhotons
    a+=oddJet
    a+=oddMuon
    a+=oddElectron
    a+=oddPhoton
    a+=badMuonInJet
    a+=secondJetET
    a+=htCut
    a+=missedHT

# -----------------------------------------------------------------------------
# Common plots 

plots_comjet = OP_ComJetPlots("CommonJetPlots",6)
plots_objkin = OP_ObjKinPlots("ObjectKinePlots",100,6)
plots_common = OP_CommonPlots("CommonPlots")
plots_kinsuite = OP_KinSuiteComPlot("KinSuitePlots",6,2)

# -----------------------------------------------------------------------------
# Common configurations ( IC5Gen, IC5Calo, AK5Calo, AK5JPT, AK5PF)

conf_ic5_gen = deepcopy(defaultConfig)
conf_ic5_gen.Ntuple = deepcopy(ic5_gen)
conf_ic5_gen.XCleaning = deepcopy(gen_cc) 
conf_ic5_gen.Common = deepcopy(gen_common)

conf_ak5_gen = deepcopy(defaultConfig)
conf_ak5_gen.Ntuple = deepcopy(ak5_gen)
conf_ak5_gen.XCleaning = deepcopy(gen_cc) 
conf_ak5_gen.Common = deepcopy(gen_common)

conf_ic5_calo = deepcopy(defaultConfig)
conf_ic5_calo.Ntuple = deepcopy(ic5_calo)
conf_ic5_calo.XCleaning = deepcopy(default_cc)
conf_ic5_calo.Common = deepcopy(default_common)

conf_ak5_calo = deepcopy(defaultConfig)
conf_ak5_calo.Ntuple = deepcopy(ak5_calo)
conf_ak5_calo.XCleaning = deepcopy(default_cc)
conf_ak5_calo.Common = deepcopy(default_common)

conf_ak5_jpt = deepcopy(defaultConfig)
conf_ak5_jpt.Ntuple = deepcopy(ak5_jpt)
conf_ak5_jpt.XCleaning = deepcopy(default_cc)
conf_ak5_jpt.Common = deepcopy(default_common)
conf_ak5_jpt.XCleaning.Muons.ModifyJetEnergy=False

conf_ak5_pf = deepcopy(defaultConfig)
conf_ak5_pf.Ntuple = deepcopy(ak5_pf)
conf_ak5_pf.XCleaning = deepcopy(default_cc)
conf_ak5_pf.Common = deepcopy(default_common)
conf_ak5_pf.XCleaning.Muons.ModifyJetEnergy=False


