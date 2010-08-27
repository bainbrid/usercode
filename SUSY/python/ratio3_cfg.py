#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Includes

from golden_cff import *

# -----------------------------------------------------------------------------
# Selection and cut flow

conf_reco = deepcopy(conf_ak5_calo)
conf_reco.XCleaning.Verbose=False
conf_reco.XCleaning.MuonJet=False
conf_reco.XCleaning.ElectronJet=False
conf_reco.XCleaning.PhotonJet=False
conf_reco.XCleaning.ResolveConflicts=False
conf_reco.Common.ApplyXCleaning=False
conf_reco.Common.Jets.EtaCut=3.0 

conf_reco20 = deepcopy(conf_reco)
conf_reco20.Common.Jets.PtCut=20.0

conf_reco30 = deepcopy(conf_reco)
conf_reco30.Common.Jets.PtCut=30.0

conf_reco40 = deepcopy(conf_reco)
conf_reco40.Common.Jets.PtCut=40.0

conf_reco50 = deepcopy(conf_reco)
conf_reco50.Common.Jets.PtCut=50.0

# Misc

numComJets = OP_NumComJets(">=",2)

ptbin = RobPlottingOps(PSet(DirName="GenPtHat",
                            PtHat=True).ps())

ratio = []
for ii in range(0,11) :
    val = 500 + ii * 5
    ratio.append( RobPlottingOps( PSet(DirName = "Ratio"+str(val),
                                       MinObjects=2,
                                       MaxObjects=8,
                                       Ratio = True,
                                       UseGen = True,
                                       AlphaTcut = (val/1000.),
                                       ).ps() ) )
    
def addCutFlow(a) :
    a+=ptbin
    a+=numComJets
    for jj in range(0,len(ratio)) :
        a+=ratio[jj] 
    
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

qcd300=PSet(
    Name="QCDPythia6_Pt300",
    File="results/Skim_QCDPythia6_Pt300.root",
    CrossSection=1.256e+03,
    Format=("ICF",2),
    )

#qcd.File=["results/prescale100/Skim_QCDPythia6_Pt*.root"]
#qcd.File=["results/skim/Skim_QCDPythia6_Pt*.root"]

# -----------------------------------------------------------------------------
# Analyses

anal_reco20 = Analysis("Reco20")
addCutFlow(anal_reco20)
anal_reco20.Run("results",conf_reco20,[qcd])

anal_reco30 = Analysis("Reco30")
addCutFlow(anal_reco30)
anal_reco30.Run("results",conf_reco30,[qcd])

anal_reco40 = Analysis("Reco40")
addCutFlow(anal_reco40)
anal_reco40.Run("results",conf_reco40,[qcd])

anal_reco50 = Analysis("Reco50")
addCutFlow(anal_reco50)
anal_reco50.Run("results",conf_reco50,[qcd])
