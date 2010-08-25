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
#alphaT = RobOps(PSet(Algo="HadronicAlphaT",Cut=0.5).ps())

def addCutFlow(a) :
    a+=numComJets
    #a+=alphaT

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

if (False) :
    new_path_qcd_pythia="/vols/cms02/bainbrid/SUSYv2/bainbrid/results/prescale100/"
    new_prefix="Skim_QCDPythia6_Pt"
    new_suffix=""
    for bin in range(0,len(qcd6)) :
        print qcd6[bin].File
        qcd6[bin].File = qcd6[bin].File.replace(path_qcd_pythia,new_path_qcd_pythia)
        qcd6[bin].File = qcd6[bin].File.replace(prefix,new_prefix)
        qcd6[bin].File = qcd6[bin].File.replace(suffix,new_suffix)
        print qcd6[bin].File
        
# -----------------------------------------------------------------------------
# Definition of analyses

for bin in range(0,len(qcd6)) :
    anal = Analysis("Empty")
    addCutFlow(anal)
    skim = SkimOp( skim_ps.ps() )
    anal += skim
    anal.Run("results",conf,[qcd6[bin]])

