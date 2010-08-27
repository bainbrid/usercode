#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Includes

from golden_cff import *

# -----------------------------------------------------------------------------
# Samples

from samples_cff import *

if (False) :
    new_path_qcd_pythia="/vols/cms02/bainbrid/SUSYv2/bainbrid/results/prescale100/"
    new_prefix="Skim_QCDPythia6_Pt"
    new_suffix=""
    qcd.File = [new_path_qcd_pythia+new_prefix+"*"+new_suffix+".root"]
    
# -----------------------------------------------------------------------------
# Turn off cross-cleaning

default_cc.Verbose=False
default_cc.MuonJet=False
default_cc.ElectronJet=False
default_cc.PhotonJet=False
default_cc.ResolveConflicts=False

# -----------------------------------------------------------------------------
# Plots

ops = RobOps(PSet(Algo="GenMet",Cut=10.).ps())

plots = RobPlottingOps( PSet( DirName="Weights",
                              PtHat=True,
                              GenMet=True ).ps() )

# -----------------------------------------------------------------------------
# Analyses

conf = deepcopy(defaultConfig)
conf.Ntuple = deepcopy(default_ntuple)
conf.XCleaning = deepcopy(default_cc)
conf.Common = deepcopy(default_common)

a=Analysis("Weights")
a+=ops
a+=plots
a.Run("results",conf,[qcd])
