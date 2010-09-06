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
conf.Common.Jets.PtCut=5.0

numComJets = OP_NumComJets(">=",2)

genMet = RobOps(PSet(Algo="GenMet",Cut=10.).ps())

ptpre = RobPlottingOps(PSet(DirName="GenPtHatPre",
                            PtHat=True).ps())

ptpost = RobPlottingOps(PSet(DirName="GenPtHatPost",
                             PtHat=True).ps())

ratio = []
for ii in range(0,11) :
    val = 500 + ii * 5
    psets = []
    for jj in range(0,3) :
        pt = 20. + jj * 10.
        ratio.append( RobPlottingOps( PSet(DirName = "Ratio"+str(val)+"Pt"+str(pt),
                                           MinObjects=2,
                                           MaxObjects=8,
                                           Ratio = True,
                                           MinJetPt = pt,
                                           UseGen = True,
                                           AlphaTcut = (val/1000.),
                                           ).ps() ) )
        
def addCutFlow(a) :
    a+=ptpre
    a+=genMet
    a+=ptpost
    a+=numComJets
    for kk in range(0,len(ratio)) :
        a+=ratio[kk] 
    
# -----------------------------------------------------------------------------
# Samples

from samples_cff import *

if (True) :
    new_path_qcd_pythia="/vols/cms02/bainbrid/SUSYv2/bainbrid/results/prescale100/"
    new_prefix="Skim_QCDPythia6_Pt"
    new_suffix=""
    qcd.File = [new_path_qcd_pythia+new_prefix+"*"+new_suffix+".root"]

# -----------------------------------------------------------------------------
# Analyses

a = Analysis("Ratio")
addCutFlow(a)
a.Run("results",conf,[qcd])
