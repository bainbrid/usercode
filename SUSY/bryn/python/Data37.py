#!/usr/bin/env python

import setupSUSY
from libFrameworkSUSY import *
from libbryn import *
from libHadronic import *
from icf.core import PSet,Analysis
from batchGolden import *
from ra1objectid.vbtfElectronId_cff import *
from ra1objectid.vbtfMuonId_cff import *
from ra1objectid.ra3PhotonId_cff import *
vbtfMuonId_cff = Muon_IDFilter( vbtfmuonidps.ps()  )
vbtfElectronIdFilter = Electron_IDFilter( vbtfelectronidWP95ps.ps() )
ra3PhotonIdFilter    = Photon_IDFilter( ra3photonidps.ps() )
PreScaleWeights = PreScaleReweighting( datatriggerps.ps() )
#  Change the settings from golden to use the lowest scaled bin.
# secondJetET = OP_SecondJetEtCut(73.3)
default_common.Jets.PtCut=36.7
cutTree,blah = MakeDataTree(73.3)
print cutTree

def addCutFlowData(a) :
  a.AddWeightFilter("Weight",PreScaleWeights)
  # a.AddJetFilter("PreCC",JetCorrections)
  a.AddPhotonFilter("PreCC",ra3PhotonIdFilter)
  a.AddElectronFilter("PreCC",vbtfElectronIdFilter)
  a.AddMuonFilter("PreCC",vbtfMuonId_cff)
  a+=cutTree

# AK5 Calo

conf_ak5_caloData = deepcopy(defaultConfig)
conf_ak5_caloData.Ntuple = deepcopy(ak5_calo)
conf_ak5_caloData.XCleaning = deepcopy(default_cc)
conf_ak5_caloData.Common = deepcopy(default_common)
# conf_ak5_calo.Common.print_out()
anal_ak5_caloData=Analysis("AK5Calo")
addCutFlowData(anal_ak5_caloData)

# AK5 PF

conf_ak5_pfData = deepcopy(defaultConfig)
conf_ak5_pfData.Ntuple = deepcopy(ak5_pf)
conf_ak5_pfData.XCleaning = deepcopy(default_cc)
conf_ak5_pfData.Common = deepcopy(default_common)
anal_ak5_pfData=Analysis("AK5PF")
addCutFlowData(anal_ak5_pfData)


tedSkim = PSet(
Name="TedSkim",
Format=("ICF",3),
Weight= 1.0,
File="../resultsWithSingleTop//Data/AK5Calo_Jets.root")
#"~elaird1/public_html/73_candidates/v3/350_bin/calo.root")


from data.Run2011.HT_Run2011_promptReco_DCS import *
from data.Run2011.HT42_incomplete import *
from data.Run2011.RA1ToBurn import *
outDir = "../results/Data37"
ensure_dir(outDir)
# HT42_incomplete.LastEntry = 100
# HT42_incomplete.File = HT42_incomplete.File[0:1]
anal_ak5_caloData.Run(outDir,conf_ak5_caloData,[HT42_incomplete])

# from data.MultiJet_Run2010B_Nov4ReReco_v1 import *

#
# anal_ak5_pfData.Run("../results/Data",conf_ak5_pfData,[HT_Run2011_promptReco_DCS])
# anal_ak5_jptData.Run("../results/",conf_ak5_jptData,data)
# anal_ak7_caloData.Run("../results/",conf_ak7_caloData,data)
