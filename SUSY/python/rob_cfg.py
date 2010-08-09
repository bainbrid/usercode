#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Includes

from golden_cff import *

# -----------------------------------------------------------------------------
# Local samples

from copy import deepcopy

if ( True ) :
    
    new_path = "/data2/bainbrid/data/7TeV/V00-08-04-XX/"
    
    lm0.File = lm0.File.replace(path_lmx,new_path+"LMx/")
    lm1.File = lm1.File.replace(path_lmx,new_path+"LMx/")
    
    qcd_pythia_merged.File = qcd_pythia_merged.File.replace(path_qcd_pythia_merged,new_path+"QCDPythia/")
    qcd_pythia_15.File = qcd_pythia_15.File.replace(path_qcd_pythia,new_path+"QCDPythia/")
    qcd_pythia_30.File = qcd_pythia_30.File.replace(path_qcd_pythia,new_path+"QCDPythia/")
    qcd_pythia_80.File = qcd_pythia_80.File.replace(path_qcd_pythia,new_path+"QCDPythia/")
    qcd_pythia_170.File = qcd_pythia_170.File.replace(path_qcd_pythia,new_path+"QCDPythia/")
    qcd_pythia_300.File = qcd_pythia_300.File.replace(path_qcd_pythia,new_path+"QCDPythia/")
    qcd_pythia_470.File = qcd_pythia_470.File.replace(path_qcd_pythia,new_path+"QCDPythia/")
    qcd_pythia_800.File = qcd_pythia_800.File.replace(path_qcd_pythia,new_path+"QCDPythia/")
    
    w_jets.File = w_jets.File.replace(path_v_jets,new_path+"VJets/")
    z_jets.File = z_jets.File.replace(path_v_jets,new_path+"VJets/")
    ttbar_jets.File = ttbar_jets.File.replace(path_v_jets,new_path+"VJets/")
    z_inv.File = z_inv.File.replace(path_z_inv,new_path+"VJets/")
    
    qcd_pythia_merged_example=PSet(
        Name="QCD_Pythia_Merged",
        Format=("ICF",2),
        File=new_path+"QCDPythia/examples/QCDJets_Pythia.root",
        Weights = PSet(
        CrossSection = [ 8.762e+08, 6.041e+07, 9.238e+05, 2.547e+04, 1.256e+03, 8.798e+01, 2.186e+00, 1.122e-02 ],
        Events       = [ 10000,     10000,     10000,     10000,     10000,     10000,     10000,     10000     ],
        PtBin        = [ 15.,       30.,       80.,       170.,      300.,      470.,      800.,      1400.     ]
        )
        )
    
    w_jets_example = deepcopy(w_jets)
    w_jets_example.File = w_jets_example.File.replace(new_path+"VJets/",new_path+"VJets/examples/")
    z_jets_example = deepcopy(z_jets)
    z_jets_example.File = z_jets_example.File.replace(new_path+"VJets/",new_path+"VJets/examples/")
    ttbar_jets_example = deepcopy(ttbar_jets)
    ttbar_jets_example.File = ttbar_jets_example.File.replace(new_path+"VJets/",new_path+"VJets/examples/")
    z_inv_example = deepcopy(z_inv)
    z_inv_example.File = z_inv_example.File.replace(new_path+"VJets/",new_path+"VJets/examples/")

else :

    new_path = "~/Documents/WORK/susy/data/7TeV/V00-08-04-XX/"
    lm0.File = lm0.File.replace(path_lmx,new_path)
    lm1.File = lm1.File.replace(path_lmx,new_path)
    qcd_pythia_470.File = qcd_pythia_470.File.replace(path_qcd_pythia,new_path)
    z_inv.File = w_jets.File.replace(path_v_jets,new_path)
    
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
