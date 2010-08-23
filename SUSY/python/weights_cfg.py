#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Includes

from golden_cff import *

# -----------------------------------------------------------------------------
# Samples

from montecarlo.LMx import *
from qcd_cff import *

# -----------------------------------------------------------------------------
# Turn off cross-cleaning

default_cc.Verbose=False
default_cc.MuonJet=False
default_cc.ElectronJet=False
default_cc.PhotonJet=False
default_cc.ResolveConflicts=False

# -----------------------------------------------------------------------------
# Plots

ps = PSet( DirName = "Weights", PtHat = True, )
plots = RobPlottingOps( ps.ps() )

# -----------------------------------------------------------------------------
# Analyses

conf = deepcopy(defaultConfig)
conf.Ntuple = deepcopy(default_ntuple)
conf.XCleaning = deepcopy(default_cc)
conf.Common = deepcopy(default_common)

a=Analysis("Weights")
a+=plots
#a.Run("../results",conf,[LM0])
#a.Run("../results",conf,[LM1])
a.Run("../results",conf,[qcd_pythia6])
#a.Run("../results",conf,qcd_pythia8_merged)
