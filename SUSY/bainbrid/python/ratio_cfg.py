#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Includes

from golden_cff import *

# -----------------------------------------------------------------------------
# Selection and cut flow

# Gen

conf_gen = deepcopy(conf_ak5_gen)
conf_gen.Common.Jets.EtaCut=3.0

conf_gen20 = deepcopy(conf_gen)
conf_gen20.Common.Jets.PtCut=20.0

conf_gen30 = deepcopy(conf_gen)
conf_gen30.Common.Jets.PtCut=30.0

conf_gen40 = deepcopy(conf_gen)
conf_gen40.Common.Jets.PtCut=40.0

conf_gen50 = deepcopy(conf_gen)
conf_gen50.Common.Jets.PtCut=50.0

# Reco

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

numComJets = OP_NumComJets("==",2)
alphaT50 = OP_CommonAlphaTCut(0.50)
alphaT51 = OP_CommonAlphaTCut(0.51)
alphaT52 = OP_CommonAlphaTCut(0.52)
alphaT53 = OP_CommonAlphaTCut(0.53)
alphaT54 = OP_CommonAlphaTCut(0.54)
alphaT55 = OP_CommonAlphaTCut(0.55)

ptbin = RobPlottingOps(PSet(DirName="GenPtHat",
                            PtHat=True).ps())

pre = RobPlottingOps(PSet(DirName="PreAlphaT",
                          HT=True,
                          GenHT=True,
                          Meff=True,
                          GenMeff=True,
                          ).ps())

post50 = RobPlottingOps(PSet(DirName="PostAlphaT50",
                             HT=True,
                             GenHT=True,
                             Meff=True,
                             GenMeff=True,
                             ).ps())

post51 = RobPlottingOps(PSet(DirName="PostAlphaT51",
                             HT=True,
                             GenHT=True,
                             Meff=True,
                             GenMeff=True,
                             ).ps())

post52 = RobPlottingOps(PSet(DirName="PostAlphaT52",
                             HT=True,
                             GenHT=True,
                             Meff=True,
                             GenMeff=True,
                             ).ps())

post53 = RobPlottingOps(PSet(DirName="PostAlphaT53",
                             HT=True,
                             GenHT=True,
                             Meff=True,
                             GenMeff=True,
                             ).ps())

post54 = RobPlottingOps(PSet(DirName="PostAlphaT54",
                             HT=True,
                             GenHT=True,
                             Meff=True,
                             GenMeff=True,
                             ).ps())

post55 = RobPlottingOps(PSet(DirName="PostAlphaT55",
                             HT=True,
                             GenHT=True,
                             Meff=True,
                             GenMeff=True,
                             ).ps())

def addCutFlow(a) :
    a+=ptbin
    a+=numComJets
    a+=pre
    a+=alphaT50
    a+=post50
    a+=alphaT51
    a+=post51
    a+=alphaT52
    a+=post52
    a+=alphaT53
    a+=post53
    a+=alphaT54
    a+=post54
    a+=alphaT55
    a+=post55
    
# -----------------------------------------------------------------------------
# Samples

qcd=PSet(
    Name="QCDPythia6",
    File=["results/prescale100/Skim_QCDPythia6_Pt*.root"],
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

qcd.File=["/vols/cms02/bm409/QCD_Pythia6_*GeV.root"]
#qcd.File=["results/All/Skim_QCDPythia6_Pt*.root"],
#qcd.File=["results/prescale100/Skim_QCDPythia6_Pt170.root"]

# -----------------------------------------------------------------------------
# Analyses

anal_gen20 = Analysis("Gen20")
addCutFlow(anal_gen20)
anal_gen20.Run("results",conf_gen20,[qcd])

anal_gen30 = Analysis("Gen30")
addCutFlow(anal_gen30)
anal_gen30.Run("results",conf_gen30,[qcd])

anal_gen40 = Analysis("Gen40")
addCutFlow(anal_gen40)
anal_gen40.Run("results",conf_gen40,[qcd])

anal_gen50 = Analysis("Gen50")
addCutFlow(anal_gen50)
anal_gen50.Run("results",conf_gen50,[qcd])

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



