#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Includes

from golden_cff import *

# -----------------------------------------------------------------------------
# Local samples

qcd_path = "/vols/cms02/gouskos/"
qcd_pythia6 = PSet(
    Name = "QCD_Pythia6",
    Format = ("ICF",2),
    File = (
    qcd_path + "QCD_Pythia_Pt15_Jun2010/SusyCAF_Tree_*.root",
    qcd_path + "QCD_Pythia_Pt30_Jun2010/SusyCAF_Tree_*.root",
    qcd_path + "QCD_Pythia_Pt80_Jun2010/SusyCAF_Tree_*.root",
    qcd_path + "QCD_Pythia_Pt170_Jun2010/SusyCAF_Tree_*.root",
    qcd_path + "QCD_Pythia_Pt300_Jun2010/SusyCAF_Tree_*.root",
    qcd_path + "QCD_Pythia_Pt470_Jun2010/SusyCAF_Tree_*.root",
    qcd_path + "QCD_Pythia_Pt800_Jun2010/SusyCAF_Tree_*.root",
    qcd_path + "QCD_Pythia_Pt1400_Jun2010/SusyCAF_Tree_*.root",
    ),
    Weights = PSet(
    CrossSection = [ 8.762e+08, 6.041e+07, 9.238e+05, 2.547e+04, 1.256e+03, 8.798e+01, 2.186, 0.01122 ],
    Events       = [ 6095857, 5069664, 2065792, 3171950, 2976108, 2159497, 2181700, 1185024 ],
    PtBin        = [ 15., 30., 80., 170., 300., 470., 800., 1400. ],
    )
    )
    
# -----------------------------------------------------------------------------
# Cuts and plots

HadAlphaT = HadronicAlphaT(0.55)

pset1 = PSet(
    DirName    = "Hadronic",
    MinObjects = 2,
    MaxObjects = 10,
    Verbose    = False,
    Summary    = False,
    CC         = False,
    Dalitz     = True,
    AlphaT     = True,
    PtHat      = False,
    MET        = False,
    Kine       = True,
    Response   = False,
    )
plots_hadronic = HadronicPlottingOps( pset1.ps() )

pset2 = PSet(
    DirName    = "Response",
    MinObjects = 2,
    MaxObjects = 10,
    Verbose    = True,
    Summary    = True,
    CC         = True,
    Dalitz     = False,
    AlphaT     = False,
    PtHat      = False,
    MET        = False,
    Kine       = False,
    Response   = True,
    )
debug_hadronic = HadronicPlottingOps( pset2.ps() )

# -----------------------------------------------------------------------------
# Cut flow

def addCutFlow(a) :
    a+=numComJets
    a+=numComLeptons
    a+=numComPhotons
    a+=oddJet
    a+=oddMuon
    a+=oddElectron
    a+=oddPhoton
    a+=badMuonInJet
    a+=secondJetET
    a+=htCut
    a+=missedHT
    a+=plots_kinsuite
    a+=plots_hadronic
    a+=alphaT
    a+=plots_common
    a+=debug_hadronic
    
# -----------------------------------------------------------------------------
# Definition of analyses

anal_ic5_calo=Analysis("IC5Calo")
addCutFlow(anal_ic5_calo)

anal_ak5_calo=Analysis("AK5Calo")
addCutFlow(anal_ak5_calo)

anal_ak5_jpt=Analysis("AK5JPT")
addCutFlow(anal_ak5_jpt)

anal_ak5_pf=Analysis("AK5PF")
addCutFlow(anal_ak5_pf)

# -----------------------------------------------------------------------------
# Run analyses

if ( True ) :
    anal_ic5_calo.Run("results",conf_ic5_calo,[lm0])
    anal_ic5_calo.Run("results",conf_ic5_calo,[lm1])
    anal_ic5_calo.Run("results",conf_ic5_calo,[qcd_pythia_merged_example])
    anal_ic5_calo.Run("results",conf_ic5_calo,[w_jets_example])
    anal_ic5_calo.Run("results",conf_ic5_calo,[z_jets_example])
    anal_ic5_calo.Run("results",conf_ic5_calo,[ttbar_jets_example])
    anal_ic5_calo.Run("results",conf_ic5_calo,[z_inv_example])

all = False;
    
if ( all ) :
    anal_ic5_calo.Run("results",conf_ic5_calo,[lm0])
    anal_ic5_calo.Run("results",conf_ic5_calo,[lm1])
    anal_ic5_calo.Run("results",conf_ic5_calo,[qcd_pythia_merged])
    anal_ic5_calo.Run("results",conf_ic5_calo,[w_jets])
    anal_ic5_calo.Run("results",conf_ic5_calo,[z_jets])
    anal_ic5_calo.Run("results",conf_ic5_calo,[ttbar_jets])
    anal_ic5_calo.Run("results",conf_ic5_calo,[z_inv])
    
if ( all ) :
    anal_ak5_calo.Run("results",conf_ak5_calo,[lm0])
    anal_ak5_calo.Run("results",conf_ak5_calo,[lm1])
    anal_ak5_calo.Run("results",conf_ak5_calo,[qcd_pythia_merged])
    anal_ak5_calo.Run("results",conf_ak5_calo,[w_jets])
    anal_ak5_calo.Run("results",conf_ak5_calo,[z_jets])
    anal_ak5_calo.Run("results",conf_ak5_calo,[ttbar_jets])
    anal_ak5_calo.Run("results",conf_ak5_calo,[z_inv])

if ( all ) :
    anal_ak5_jpt.Run("results",conf_ak5_jpt,[lm0])
    anal_ak5_jpt.Run("results",conf_ak5_jpt,[lm1])
    anal_ak5_jpt.Run("results",conf_ak5_jpt,[qcd_pythia_merged])
    anal_ak5_jpt.Run("results",conf_ak5_jpt,[w_jets])
    anal_ak5_jpt.Run("results",conf_ak5_jpt,[z_jets])
    anal_ak5_jpt.Run("results",conf_ak5_jpt,[ttbar_jets])
    anal_ak5_jpt.Run("results",conf_ak5_jpt,[z_inv])

if ( all ) :
    anal_ak5_pf.Run("results",conf_ak5_pf,[lm0])
    anal_ak5_pf.Run("results",conf_ak5_pf,[lm1])
    anal_ak5_pf.Run("results",conf_ak5_pf,[qcd_pythia_merged])
    anal_ak5_pf.Run("results",conf_ak5_pf,[w_jets])
    anal_ak5_pf.Run("results",conf_ak5_pf,[z_jets])
    anal_ak5_pf.Run("results",conf_ak5_pf,[ttbar_jets])
    anal_ak5_pf.Run("results",conf_ak5_pf,[z_inv])
