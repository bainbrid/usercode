#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Includes

from golden_cff import *

# -----------------------------------------------------------------------------
# Turn off cross-cleaning

default_cc.Verbose=False
default_cc.MuonJet=False
default_cc.ElectronJet=False
default_cc.PhotonJet=False
default_cc.ResolveConflicts=False

# -----------------------------------------------------------------------------
# Plots

ps = PSet(
    DirName    = "Weights",
    MinObjects = 2,
    MaxObjects = 10,
    Verbose    = False,
    Summary    = False,
    CC         = False,
    Dalitz     = False,
    AlphaT     = False,
    PtHat      = True,
    MET        = False,
    Kine       = False,
    Response   = False,
    )
plots = HadronicPlottingOps( ps.ps() )
    
# -----------------------------------------------------------------------------
# Analyses

conf = deepcopy(defaultConfig)
conf.Ntuple = deepcopy(default_ntuple)
conf.XCleaning = deepcopy(default_cc)
conf.Common = deepcopy(default_common)

a=Analysis("Weights")
a+=plots
a.Run("results",conf,[qcd_pythia_merged])
