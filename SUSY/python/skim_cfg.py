#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Includes

from golden_cff import *

# -----------------------------------------------------------------------------
# Selection and cut flow

conf = deepcopy(conf_ic5_gen)
conf.Common.Jets.PtCut=10.0
conf.Common.Jets.EtaCut=5.0

numComJets = OP_NumComJets(">=",2)
alphaT = OP_CommonAlphaTCut(0.50)

def addCutFlow(a) :
    a+=numComJets
    a+=alphaT

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
# Definition of analyses

for bin in range(0,len(qcd6)) :
    anal = Analysis(qcd6[bin].Name)
    addCutFlow(anal)
    skim = SkimOp( skim_ps.ps() )
    anal += skim
    anal.Run("../results",conf,[qcd6[bin]])

