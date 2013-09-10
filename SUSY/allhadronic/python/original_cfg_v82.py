#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Necessary includes

import setupSUSY
from libFrameworkSUSY import *
from libHadronic import *
from icf.core import PSet,Analysis
from icf.config import defaultConfig
from icf.utils import json_to_pset
from copy import deepcopy
import commands
pwd = commands.getoutput('echo $SUSY_WORKING_SW_DIR')

# -----------------------------------------------------------------------------
# Options

# 0=2010, 1=2011, 2=QCDpythia, 3=QCDmagdraph, 4=Wjets, 5=Zinv, 6=TTbar, 7=SingleTop, 8=LM1, 10=test
dataset = 1

output_dir = pwd + "/results/v82"

FilterOption = 0
OffsetOption = 0.

DataOption = False
if ( dataset == 0 or dataset == 1 or dataset == 10 ) :
    DataOption = True
    
JsonOption = False
if ( dataset == 1 ) :
    JsonOption = True
    
#JsonFileOption = "/allhadronic/python/Cert_160404-166502_7TeV_PromptReco_Collisions11_JSON.txt"
#JsonFileOption = "/allhadronic/python/Cert_160404-166861_7TeV_PromptReco_Collisions11_JSON.txt"
JsonFileOption = "/allhadronic/python/Cert_160404-167151_7TeV_PromptReco_Collisions11_JSON.txt"

VetoesOption = "NoEGamma"
AddMuonOption = False
if ( VetoesOption == "NoMuon" ) :
    AddMuonOption = True

SignalTriggersOption = [""]
MonitorTriggersOption = [""]
if ( dataset == 1 or dataset == 10 ) :
    SignalTriggersOption = [
        "HLT_HT250_MHT60_v2",
        "HLT_HT250_MHT60_v3",
        "HLT_HT250_MHT60_v4",
        "HLT_HT250_MHT60_v6",
        "HLT_HT250_MHT60_v7",
        "HLT_HT250_MHT70_v1",
        "HLT_HT250_MHT70_v3",
        "HLT_HT250_MHT70_v4",
        "HLT_HT250_MHT80_v3",
        "HLT_HT250_MHT80_v4",
        "HLT_HT250_MHT90_v1",
        "HLT_HT260_MHT60_v2",
        ""
        ]
    MonitorTriggersOption = [
        # 275
        "HLT_HT260_v2",
        "HLT_HT250_v2",
        "HLT_HT250_v3",
        "HLT_HT250_v4",
        "HLT_HT250_v5",
        "HLT_HT250_v6",
        "HLT_HT250_v7",
        "",
        # 325
        "HLT_HT300_v2",
        "HLT_HT300_v3",
        "HLT_HT300_v4",
        "HLT_HT300_v5",
        "HLT_HT300_v6",
        "HLT_HT300_v7",
        "HLT_HT300_v8",
        "",
        # 375
        "HLT_HT300_v2",
        "HLT_HT350_v2",
        "HLT_HT350_v3",
        "HLT_HT350_v4",
        "HLT_HT350_v5",
        "HLT_HT350_v6",
        "HLT_HT350_v7",
        "", 
        # 475
        "HLT_HT440_v2",
        "HLT_HT450_v2",
        "HLT_HT450_v3",
        "HLT_HT450_v4",
        "HLT_HT450_v5",
        "HLT_HT450_v6",
        "HLT_HT450_v7",
        "",
        # 575
        "HLT_HT520_v2",
        "HLT_HT550_v2",
        "HLT_HT550_v3",
        "HLT_HT550_v4",
        "HLT_HT550_v5",
        "HLT_HT550_v6",
        "HLT_HT550_v7",
        "",
        # 675
        "HLT_HT520_v2",
        "HLT_HT550_v2",
        "HLT_HT550_v3",
        "HLT_HT550_v4",
        "HLT_HT550_v5",
        "HLT_HT550_v6",
        "HLT_HT550_v7",
        "",
        # 775
        "HLT_HT520_v2",
        "HLT_HT550_v2",
        "HLT_HT550_v3",
        "HLT_HT550_v4",
        "HLT_HT550_v5",
        "HLT_HT550_v6",
        "HLT_HT550_v7",
        "",
        # 875
        "HLT_HT520_v2",
        "HLT_HT550_v2",
        "HLT_HT550_v3",
        "HLT_HT550_v4",
        "HLT_HT550_v5",
        "HLT_HT550_v6",
        "HLT_HT550_v7",
        "",
        ]
    
NoQcdOption = False
if ( dataset == 2 or dataset == 3 ) :
    NoQcdOption = True
    
# -----------------------------------------------------------------------------
# Reading the collections from the ntuple

default_ntuple = deepcopy(defaultConfig.Ntuple)
default_ntuple.Electrons=PSet(
Prefix="electron",
Suffix="Pat",
LooseID="EIDLoose",
TightID="EIDTight",
)

default_ntuple.Muons=PSet(
Prefix="muon",
Suffix="Pat",
LooseID="IsGlobalMuon",
TightID="IDGlobalMuonPromptTight",
)

default_ntuple.Taus=PSet(
Prefix="tau",
Suffix="Pat",
LooseID="TauIdbyTaNCfrOnePercent",
TightID="TauIdbyTaNCfrTenthPercent"
)

default_ntuple.Photons=PSet(
Prefix="photon",
Suffix="Pat",
)

default_ntuple.Jets=PSet(
Prefix="ic5Jet",
Suffix="Pat",
Uncorrected=False,
)

default_ntuple.SecJets=PSet(
    Prefix="",
    Suffix="",
    Uncorrected=False,
    )

default_ntuple.TerJets=PSet(
    Prefix="",
    Suffix="",
    Uncorrected=False,
    )

ak5_calo = deepcopy(default_ntuple)
ak5_calo.Jets.Prefix="ak5Jet"

# -----------------------------------------------------------------------------
# Cross-cleaning settings

default_cc = deepcopy(defaultConfig.XCleaning)
default_cc.Verbose=False
default_cc.MuonJet=True
default_cc.ElectronJet=True
default_cc.PhotonJet=True
default_cc.ResolveConflicts=True

default_cc.Jets.PtCut=10.0
default_cc.Jets.EtaCut=10.0

default_cc.Muons.ModifyJetEnergy=True
default_cc.Muons.PtCut=10.0
default_cc.Muons.EtaCut=2.5
default_cc.Muons.TrkIsoCut=-1.
default_cc.Muons.CombIsoCut=0.15
default_cc.Muons.MuonJetDeltaR=0.5
default_cc.Muons.MuonIsoTypePtCutoff=0.
default_cc.Muons.RequireLooseIdForInitialFilter=False

default_cc.Electrons.PtCut=10.0
default_cc.Electrons.EtaCut=2.5
default_cc.Electrons.TrkIsoCut=-1.0
default_cc.Electrons.CombIsoCut=0.15
default_cc.Electrons.ElectronJetDeltaR=0.5
default_cc.Electrons.ElectronIsoTypePtCutoff=0.   
default_cc.Electrons.ElectronLooseIdRequired=True
default_cc.Electrons.ElectronTightIdRequired=False
default_cc.Electrons.RequireLooseIdForInitialFilter=False

default_cc.Photons.EtCut=25.0
default_cc.Photons.EtaCut=2.5
default_cc.Photons.TrkIsoCut=2.0
default_cc.Photons.CaloIsoCut=0.2
default_cc.Photons.IDReq=3
default_cc.Photons.PhotonJetDeltaR=0.5
default_cc.Photons.PhotonIsoTypePtCutoff=30.

# -----------------------------------------------------------------------------
# Definition of common objects

default_common = deepcopy(defaultConfig.Common)

default_common.ApplyXCleaning=True

default_common.Jets.PtCut=10.0
default_common.Jets.EtaCut=3.0
default_common.Jets.ApplyID=True
default_common.Jets.TightID=False

default_common.Electrons.PtCut=10.0
default_common.Electrons.EtaCut=2.5
default_common.Electrons.TrkIsoCut=-1.
default_common.Electrons.CombIsoCut=0.15
default_common.Electrons.ApplyID = True
default_common.Electrons.TightID = False
default_common.Electrons.RequireLooseForOdd = True

default_common.Muons.PtCut=10.0
default_common.Muons.EtaCut=2.5
default_common.Muons.TrkIsoCut=-1.
default_common.Muons.CombIsoCut=0.15
default_common.Muons.ApplyID = True
default_common.Muons.TightID = True
default_common.Muons.RequireLooseForOdd = True

default_common.Photons.EtCut=25.0
default_common.Photons.EtaCut=2.5
default_common.Photons.IDReq = 3
default_common.Photons.RequireLooseForOdd = True

# -----------------------------------------------------------------------------

# Common configurations

conf_ak5_calo = deepcopy(defaultConfig)
conf_ak5_calo.Ntuple = deepcopy(ak5_calo)
conf_ak5_calo.XCleaning = deepcopy(default_cc)
conf_ak5_calo.Common.print_out()
conf_ak5_calo.Common = deepcopy(default_common)

# -----------------------------------------------------------------------------
# Common cut definitions

from ra1objectid.vbtfMuonId_cff import *
from ra1objectid.vbtfElectronId_cff import *
from ra1objectid.ra3PhotonId_cff import *

vbtf_muon = Muon_IDFilter( vbtfmuonidps.ps() )
vbtf_electron = Electron_IDFilter( vbtfelectronidWP95ps.ps() )
ra3_photon = Photon_IDFilter( ra3photonidps.ps() )

from icf.JetCorrections import *
corPset =  CorrectionPset("ResidualJetEnergyCorrections.txt")
JetCorrections = JESCorrections( corPset.ps(), True )

NoiseFilt= OP_HadronicHBHEnoiseFilter()
selection = OP_GoodEventSelection()
selection_duplicate = OP_GoodEventSelection()

oddMuon = OP_OddMuon()
oddElectron = OP_OddElectron()
oddPhoton = OP_OddPhoton()

numComElectrons = OP_NumComElectrons("<=",0)
numComMuons = OP_NumComMuons("<=",0)
numComPhotons = OP_NumComPhotons("<=",0)

DeadEcalFileOption = pwd+"/hadronic/python/deadRegionList_START38_V12.txt"
CutFlowOption="MC"
if ( DataOption ) :
    pwd+"/hadronic/python/deadRegionList_GR10_P_V10.txt"
    CutFlowOption="Data"
    
vertex_reweight = VertexReweighting(
    PSet(
    VertexWeights =[0.20, 0.63, 1.19, 1.57, 1.62, 1.42, 1.09, 0.80 ,0.57, 0.42, 0.30, 0.20]
    ).ps())

HtBinOption = 0
HtBinningOptions = []
if ( True ) :
    HtBinOption = 2
    HtBinningOptions = [275.,325.,375.,475.,575.,675.,775.,875.,975.,]
else :
    for n in range(0,31) :
        HtBinOption = 5
        HtBinningOptions.append(250+20.*n)

qcdBkgdEst = QcdBkgdEst( PSet(DirName = "QcdBkgdEst",
                              # Triggers
                              SignalTriggers = SignalTriggersOption,
                              MonitorTriggers = MonitorTriggersOption,
                              # Miscellaneous
                              Verbose=False,
                              Filter=FilterOption,
                              NoQCD=NoQcdOption,
                              MaxLeadingJetEta=2.5,
                              AddMuon=AddMuonOption,
                              # Histogramming
                              MinObjects=2,
                              MaxObjects=4,
                              HtHistos=False,
                              MhtHistos=False,
                              MultiplicityHistos=False,
                              BabyJetsHistos=False,
                              MinBiasDPhiHistos=False,
                              # Baby jets cleaning
                              NoCleaningInDenom=True,
                              BabyJetsCut=1.25,
                              # Dead ECAL cleaning
                              DeadEcalRegionDrCut=0.3,
                              DeadEcalRegionEtaCut=0.3,
                              MinBiasDeltaPhiCut=0.5,
                              MinJetPtCut=30.,
                              NBadCells=10,
                              UseDeadEcalFile=False,
                              DeadEcalFile=DeadEcalFileOption,
                              # RecHit cleaning
                              SumRecHitPtCut = 30.,
                              # Binning
                              AlphaTDefault=0.55,
                              AlphaT=[0.51,0.55,0.60,
                                      #0.51,0.52,
                                      #0.53,0.531,0.532,0.533,0.534,0.535,0.536,0.537,0.538,0.539,
                                      #0.54,0.541,0.542,0.543,0.544,0.545,0.546,0.547,0.548,0.549,
                                      #0.55,0.551,0.552,0.553,0.554,0.555,0.556,0.557,0.558,0.559,
                                      #0.56,0.561,0.562,0.563,0.564,0.565,0.566,0.567,0.568,0.569,
                                      #0.57,0.571,0.572,0.573,0.574,0.575,0.576,0.577,0.578,0.579,
                                      #0.58,0.581,0.582,0.583,0.584,0.585,0.586,0.587,0.588,0.589,
                                      #0.59,0.60
                                      ],
                              Inclusive=False,
                              HtBinning=HtBinningOptions,
                              HtNbins=len(HtBinningOptions)-1,
                              HtXlow=HtBinningOptions[0],
                              HtXhigh=HtBinningOptions[len(HtBinningOptions)-1],
                              HtOffset=0.,
                              # Thresholds and scaling
                              MinJetPt1=100.,
                              MinJetPt2=100.,
                              MinJetPt3=50.,
                              HtBin=HtBinOption,
                              Scale=True,
                              ScaleInSignalRegion=False,
                              ScaleSingleJetThreshold=True,
                              ).ps() )

# -----------------------------------------------------------------------------
# Cut flow

cut_flow = Tree(CutFlowOption)

if ( JsonOption ) : # Use JSON
    json = JSONFilter("JSON", json_to_pset(pwd+JsonFileOption))
    json_output = JSONOutput("_filtered")
    cut_flow.Attach(json)
    cut_flow.TAttach(json,json_output)
    cut_flow.TAttach(json,selection)
else :
    cut_flow.Attach(selection)
   
if ( DataOption ) : # Data only
    cut_flow.TAttach(selection,NoiseFilt)
    cut_flow.TAttach(NoiseFilt,selection_duplicate)
else:
    cut_flow.TAttach(selection,selection_duplicate)

mht_cut = OP_TriggerMHT_Emu(60.,30.)

if ( VetoesOption == "AllVetoes" ) :
    cut_flow.TAttach(selection_duplicate,numComPhotons)
    cut_flow.TAttach(numComPhotons,oddPhoton)
    cut_flow.TAttach(oddPhoton,numComElectrons)
    cut_flow.TAttach(numComElectrons,oddElectron)
    cut_flow.TAttach(oddElectron,numComMuons)
    cut_flow.TAttach(numComMuons,oddMuon)
    cut_flow.TAttach(oddMuon,qcdBkgdEst)
#    if (False) :
#        cut_flow.TAttach(oddMuon,mht_cut)
#        cut_flow.TAttach(mht_cut,qcdBkgdEst)
#    else :
#        cut_flow.TAttach(oddMuon,qcdBkgdEst)
elif ( VetoesOption == "NoMuon" ) :
    cut_flow.TAttach(selection_duplicate,numComPhotons)
    cut_flow.TAttach(numComPhotons,oddPhoton)
    cut_flow.TAttach(oddPhoton,numComElectrons)
    cut_flow.TAttach(numComElectrons,oddElectron)
    cut_flow.TAttach(oddElectron,qcdBkgdEst)
elif ( VetoesOption == "NoElectron" ) :
    cut_flow.TAttach(selection_duplicate,numComPhotons)
    cut_flow.TAttach(numComPhotons,oddPhoton)
    cut_flow.TAttach(oddPhoton,numComMuons)
    cut_flow.TAttach(numComMuons,oddMuon)
    cut_flow.TAttach(oddMuon,qcdBkgdEst)
elif ( VetoesOption == "NoPhoton" ) :
    cut_flow.TAttach(selection_duplicate,numComElectrons)
    cut_flow.TAttach(numComElectrons,oddElectron)
    cut_flow.TAttach(oddElectron,numComMuons)
    cut_flow.TAttach(numComMuons,oddMuon)
    cut_flow.TAttach(oddMuon,qcdBkgdEst)
elif ( VetoesOption == "NoEGamma" ) :
    cut_flow.TAttach(selection_duplicate,numComMuons)
    cut_flow.TAttach(numComMuons,oddMuon)
    cut_flow.TAttach(oddMuon,qcdBkgdEst)

if ( True and FilterOption > -1 ) :
    skim_ps = PSet(
        SkimName = "Skim",
        Branches = ["keep *"],
        DropBranches = False,
        )
    skim = SkimOp( skim_ps.ps() )
    cut_flow.TAttach(qcdBkgdEst,skim)
    
# -----------------------------------------------------------------------------
# Dataset

data2011=PSet(
    Name="Data",
    Format=("ICF",3),
    File= [
    #"/Users/bainbrid/Desktop/work/data/HT_Run2011_promptReco_DCS.root",
    "root://gfe02.grid.hep.ph.ic.ac.uk/store/user/henning//ICF/automated/2011_04_13_17_17_19///SusyCAF_Tree_9_1_zCe.root" ,
    #"file:./skim.root" ,
    #"file:/vols/cms02/bainbrid/qcd/trigger/SUSY2/results/v23/Ratio_HT42_incomplete.root",
    ],
    Weight=1.,
    #FirstEntry = 934,
    LastEntry = 10,
    )

path2010="/vols/cms02/elaird1/14_skims/MultiJet.Run2010B-PromptReco-v2.RECO.RAW.Burt2/"
data2010 = PSet(
    Name = 'JetMET_Run2010A_PromptReco',
    Format = ('ICF',2),
    File = [
    path2010+"MultiJet.Run2010B-PromptReco-v2.RECO.RAW.Burt2_0_skim.root",
    path2010+"MultiJet.Run2010B-PromptReco-v2.RECO.RAW.Burt2_1_skim.root",
    path2010+"MultiJet.Run2010B-PromptReco-v2.RECO.RAW.Burt2_2_skim.root",
    path2010+"MultiJet.Run2010B-PromptReco-v2.RECO.RAW.Burt2_3_skim.root",
    path2010+"MultiJet.Run2010B-PromptReco-v2.RECO.RAW.Burt2_4_skim.root",
    path2010+"MultiJet.Run2010B-PromptReco-v2.RECO.RAW.Burt2_5_skim.root",
    path2010+"MultiJet.Run2010B-PromptReco-v2.RECO.RAW.Burt2_6_skim.root",
    path2010+"MultiJet.Run2010B-PromptReco-v2.RECO.RAW.Burt2_7_skim.root",
    path2010+"MultiJet.Run2010B-PromptReco-v2.RECO.RAW.Burt2_8_skim.root",
    path2010+"MultiJet.Run2010B-PromptReco-v2.RECO.RAW.Burt2_9_skim.root",
    ],
    Weight=1.0,
    FirstEntry = 1,
    LastEntry = 10,
    )

candidates=PSet(
    Name="Candidates",
    Format=("ICF",3),
    File= [
    #"/vols/cms02/bainbrid/qcd/trigger/SUSY2/allhadronic/python/Ratio_HT_Run2011A_PromptReco_v1.root",
    "/vols/cms02/bainbrid/qcd/trigger/SUSY2/allhadronic/python/Ratio_HT_Run2011_promptReco_DCS-7bins.root",
    ],
    Weight=1.,
    )

from montecarlo.Spring11.QCD_TuneZ2_7TeV_pythia6_Spring11_PU_START311_ALL import *
from montecarlo.Spring11.QCD_TuneD6T_7TeV_madgraph_Spring11_PU_START311_ALL import *
from montecarlo.Spring11.WJetsToLNu_TuneZ2_7TeV_madgraph_tauola_Spring11_PU_S1_START311_V1G1_v1 import *
from montecarlo.Spring11.ZinvisibleJets_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1 import *
from montecarlo.Spring11.TTJets_TuneZ2_7TeV_madgraph_tauola_Spring11_PU_S1_START311_V1G1_v1 import *
from montecarlo.Spring11.TToBLNu_TuneZ2_t_channel_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1 import *
from montecarlo.Spring11.TToBLNu_TuneZ2_tW_channel_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1 import *
from montecarlo.Spring11.TToBLNu_TuneZ2_s_channel_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1 import *
from montecarlo.Spring11.LMx_SUSY_sftsht_7TeV_pythia6_Spring11_PU_S1_START311_V1G1_v1 import *

from data.Jet_35pb_WithTP_json051110 import Jet_35pb_WithTP_json051110
from data.Run2011.HT42_incomplete import HT42_incomplete
from data.Run2011.HT_Run2011_promptReco_DCS import HT_Run2011_promptReco_DCS
from data.Run2011.HT_Run2011A import HT_Run2011A
from data.Run2011.HT_Run2011A_PromptReco_v1 import HT_Run2011A_PromptReco_v1
from data.Run2011.HT_Run2011A_AllReco_17June import HT_Run2011A_AllReco_17June
from data.Run2011.HT_Run2011A_L1OffSet import HT_Run2011A_L1OffSet

# -----------------------------------------------------------------------------
# Analysis

anal=Analysis("Ratio")
if ( not DataOption ) :
    anal.AddWeightFilter("Weight", vertex_reweight)
anal.AddMuonFilter("PreCC",vbtf_muon)
anal.AddElectronFilter("PreCC",vbtf_electron)
anal.AddPhotonFilter("PreCC",ra3_photon)
#anal.AddJetFilter("PreCC",JetCorrections)
 
anal+=cut_flow

if   ( dataset == 0 ) : # Data 2010
    anal.Run(output_dir,conf_ak5_calo,[Jet_35pb_WithTP_json051110]) # 160 files
elif ( dataset == 1 ) : # Data 2011
    #anal.Run(output_dir,conf_ak5_calo,[HT_Run2011_promptReco_DCS]) # 1543 files (and growing)
    #anal.Run(output_dir,conf_ak5_calo,[HT42_incomplete]) # 2579 files (and growing)
    anal.Run(output_dir,conf_ak5_calo,[HT_Run2011A]) # 2579 files (and growing)
    #anal.Run(output_dir,conf_ak5_calo,[HT_Run2011A_AllReco_17June]) # 2579 files (and growing)
    #anal.Run(output_dir,conf_ak5_calo,[HT_Run2011A_L1OffSet]) # 2579 files (and growing)
elif ( dataset == 2 ) : # QCD Pythia
    anal.Run(output_dir,conf_ak5_calo,QCD_TuneZ2_7TeV_pythia6_Spring11_PU_START311_ALL) # 2462 files
elif ( dataset == 3 ) : # QCD Madgraph
    anal.Run(output_dir,conf_ak5_calo,QCD_TuneD6T_7TeV_madgraph_Spring11_PU_START311_ALL) # 1208 files
elif ( dataset == 4 ) :
    anal.Run(output_dir,conf_ak5_calo,[WJetsToLNu_TuneZ2_7TeV_madgraph_tauola_Spring11_PU_S1_START311_V1G1_v1]) # 758 files
elif ( dataset == 5 ) :
    anal.Run(output_dir,conf_ak5_calo,[ZinvisibleJets_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1]) # 110 files
elif ( dataset == 6 ) :
    anal.Run(output_dir,conf_ak5_calo,[TTJets_TuneZ2_7TeV_madgraph_tauola_Spring11_PU_S1_START311_V1G1_v1]) # 59 files
elif ( dataset == 7 ) :
    anal.Run(output_dir,conf_ak5_calo,[TToBLNu_TuneZ2_t_channel_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1,
                                       TToBLNu_TuneZ2_tW_channel_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1,
                                       TToBLNu_TuneZ2_s_channel_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1]) # 75 files
elif ( dataset == 8 ) :
    LMX = [LM1_SUSY_sftsht_7TeV_pythia6_Spring11_PU_S1_START311_V1G1_v1,
           LM6_SUSY_sftsht_7TeV_pythia6_Spring11_PU_S1_START311_V1G1_v1]
    anal.Run(output_dir,conf_ak5_calo,LMX) 
elif ( dataset == 10 ) : # Test on data 2011
    anal.Run(".",conf_ak5_calo,[data2011])
    #anal.Run(".",conf_ak5_calo,[candidates])
    #anal.Run(".",conf_ak5_calo,[data2010])
    
