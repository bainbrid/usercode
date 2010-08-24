#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Includes

from golden_cff import *

# -----------------------------------------------------------------------------
# Selection and cut flow

conf = deepcopy(conf_ak5_calo)

conf.XCleaning.Verbose=False
conf.XCleaning.MuonJet=False
conf.XCleaning.ElectronJet=False
conf.XCleaning.PhotonJet=False
conf.XCleaning.ResolveConflicts=False

conf.Common.ApplyXCleaning=False
conf.Common.Jets.EtaCut=3.0 
conf.Common.Jets.PtCut=20.0

numComJets = OP_NumComJets(">=",2)
alphaT = RobOps(PSet(Algo="HadronicAlphaT",Cut=0.5).ps())

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
    anal = Analysis("Empty")
    addCutFlow(anal)
    skim = SkimOp( skim_ps.ps() )
    anal += skim
    anal.Run("results",conf,[qcd6[bin]])

