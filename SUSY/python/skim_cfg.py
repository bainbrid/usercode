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

qcd=PSet(
    Name="QCDPythia6",
    File=["/vols/cms02/bm409/QCD_Pythia6_*GeV.root"],
    Weights = PSet(
    CrossSection = [ 8.762e+08, 6.041e+07, 9.238e+05, 2.547e+04, 1.256e+03, 8.798e+01, 2.186e+00, 1.122e-02 ],
    Events       = [ 6246300,   5228992,   3203440,   3132800,   3274202,   2143390,   2143921,   1184123   ],
    PtBin        = [ 15.,       30.,       80.,       170.,      300.,      470.,      800.,      1400.     ],
    ),
    Format=("ICF",2),
    )
#File=["results/prescale100/Skim_QCDPythia6_Pt*.root"]

# -----------------------------------------------------------------------------
# Definition of analyses

anal = Analysis("QCDPythia6")
addCutFlow(anal)
skim = SkimOp( skim_ps.ps() )
anal += skim
anal.Run("../results",conf,[qcd])

