#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Includes

from golden_cff import *

# -----------------------------------------------------------------------------
# Selection and cut flow

conf = deepcopy(conf_ic5_gen) # means no cross-cleaning!

pre_scale = RobOps(PSet(Algo="PreScale",Cut=10.).ps())

def addCutFlow(a) :
    a+=pre_scale
    
# -----------------------------------------------------------------------------
# Ntuple slimming

skim_ps = PSet(
    SkimName = "Skim",
    Branches = [""],
    DropBranches = False,
    )

# -----------------------------------------------------------------------------
# Samples

from samples_cff import *

# -----------------------------------------------------------------------------
# Definition of data sets and analyses

for bin in range(0,len(qcd6)) :
    anal = Analysis("Empty")
    addCutFlow(anal)
    skim = SkimOp( skim_ps.ps() )
    anal += skim
    anal.Run("results",conf,[qcd6[bin]])

