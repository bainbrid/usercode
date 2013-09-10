#!/usr/bin/env python

#!/usr/bin/env python
"""
Created by Bryn Mathias on 2010-05-07.
"""

# -----------------------------------------------------------------------------
# Necessary includes
import errno
import os
import setupSUSY
from libFrameworkSUSY import *
from libHadronic import *
from libbryn import *
from icf.core import PSet,Analysis
from icf.config import defaultConfig
from icf.utils import json_to_pset
from copy import deepcopy
# from icf.JetCorrections import *

# -----------------------------------------------------------------------------
# Samples
#import yours in your running script
def ensure_dir(path):
    try:
      os.makedirs(path)
    except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST:
        pass
      else: raise


# -----------------------------------------------------------------------------
# lets get some samples together!!
from montecarlo.Spring11.QCD_TuneZ2_7TeV_pythia6_Spring11_PU_START311_ALL import *
from montecarlo.Spring11.TTJets_TuneZ2_7TeV_madgraph_tauola_Spring11_PU_S1_START311_V1G1_v1 import *
from montecarlo.Spring11.WJetsToLNu_TuneZ2_7TeV_madgraph_tauola_Spring11_PU_S1_START311_V1G1_v1 import *
from montecarlo.Spring11.ZinvisibleJets_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1 import *
from montecarlo.Spring11.LMx_SUSY_sftsht_7TeV_pythia6_Spring11_PU_S1_START311_V1G1_v1 import *
from montecarlo.Spring11.GJets_TuneD6T_HT_200_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1 import *
from montecarlo.Spring11.TToBLNu_TuneZ2_t_channel_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1 import *
from montecarlo.Spring11.TToBLNu_TuneZ2_tW_channel_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1 import *
from montecarlo.Spring11.TToBLNu_TuneZ2_s_channel_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1 import *

MC = QCD_TuneZ2_7TeV_pythia6_Spring11_PU_START311_ALL+[TTJets_TuneZ2_7TeV_madgraph_tauola_Spring11_PU_S1_START311_V1G1_v1,WJetsToLNu_TuneZ2_7TeV_madgraph_tauola_Spring11_PU_S1_START311_V1G1_v1,ZinvisibleJets_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1,LM0_SUSY_sftsht_7TeV_pythia6_Spring11_PU_S1_START311_V1G1_v1,LM1_SUSY_sftsht_7TeV_pythia6_Spring11_PU_S1_START311_V1G1_v1, LM2_SUSY_sftsht_7TeV_pythia6_Spring11_PU_S1_START311_V1G1_v1, LM3_SUSY_sftsht_7TeV_pythia6_Spring11_PU_S1_START311_V1G1_v1, LM4_SUSY_sftsht_7TeV_pythia6_Spring11_PU_S1_START311_V1G1_v1, LM5_SUSY_sftsht_7TeV_pythia6_Spring11_PU_S1_START311_V1G1_v1, LM6_SUSY_sftsht_7TeV_pythia6_Spring11_PU_S1_START311_V1G1_v1, LM7_SUSY_sftsht_7TeV_pythia6_Spring11_PU_S1_START311_V1G1_v1, LM8_SUSY_sftsht_7TeV_pythia6_Spring11_PU_S1_START311_V1G1_v1, LM9_SUSY_sftsht_7TeV_pythia6_Spring11_PU_S1_START311_V1G1_v1,
GJets_TuneD6T_HT_200_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1,
TToBLNu_TuneZ2_s_channel_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1,
TToBLNu_TuneZ2_tW_channel_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1,
TToBLNu_TuneZ2_t_channel_7TeV_madgraph_Spring11_PU_S1_START311_V1G1_v1
]

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
default_ntuple.SecMuons=PSet(
    Prefix="muon",
    Suffix="PF")
default_ntuple.Taus=PSet(
Prefix="tau",
Suffix="Pat",
LooseID="TauIdbyTaNCfrOnePercent",
TightID="TauIdbyTaNCfrTenthPercent"
)
default_ntuple.Jets=PSet(
Prefix="ic5Jet",
Suffix="Pat",
Uncorrected=False,
)
default_ntuple.Photons=PSet(
Prefix="photon",
Suffix="Pat",
)

ic5_calo = deepcopy(default_ntuple)
ic5_calo.Jets.Prefix="ic5Jet"

ak5_calo = deepcopy(default_ntuple)
ak5_calo.Jets.Prefix="ak5Jet"

ak5_jpt = deepcopy(default_ntuple)
ak5_jpt.Jets.Prefix="ak5JetJPT"

ak5_pf = deepcopy(default_ntuple)
ak5_pf.Jets.Prefix="ak5JetPF"
ak5_pf.TerJets.Prefix="ak5Jet"

ak7_calo = deepcopy(default_ntuple)
ak7_calo.Jets.Prefix="ak7Jet"


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
default_cc.Photons.UseID=True
default_cc.Photons.PhotonJetDeltaR=0.5
default_cc.Photons.PhotonIsoTypePtCutoff=30.
# -----------------------------------------------------------------------------
# Definition of common objects
default_common = deepcopy(defaultConfig.Common)
default_common.ApplyXCleaning=True
default_common.Jets.PtCut=50.0
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
# default_common.Photons.EtaCut=2.5
default_common.Photons.UseID=True
# the photon cuts are NOT read anyway
# default_common.Photons.TrkIsoRel=0.
# default_common.Photons.TrkIsoCut=99999.
# default_common.Photons.EcalIsoRel=0.
# default_common.Photons.EcalIsoCut=99999.
# default_common.Photons.HcalIsoRel=0.
# default_common.Photons.HcalIsoCut=99999.
# default_common.Photons.HadOverEmCut=0.5
# default_common.Photons.SigmaIetaIetaCut=0.5
##default_common.Photons.CaloIsoCut=99999.
default_common.Photons.IDReq = 3
default_common.Photons.RequireLooseForOdd = True

Events = PSet(
Run =[160939,160955,160957,160957,162811,162811,162909,162909,163046,163069,163069,163069,163071,163071,163071,163078,163078,163235,163235,163235,163235,163235,163235,163235,163255,163255,163255,163255,163255,163261,163261,163270,163270,163270,163270,163270,163270,163270,163270,163270,163270,163270,163270,163270,163270,163270,163270,163270,163286,163286,163286,163286,163286,163296,163296,163296,163296,163297,163297,163297,163300,163300,163300,163300,163300,163300,163300,163301,163301,163302,163302,163302,163302,163302,163302,163332,163332,163332,163332,163332,163332,163333,163334,163334,163334,163334,163334,163337,163337,163337,163337,163337,163337,163338,163338,163339,163339,163339,163340,163340,163340,163340,163340,163369,163370,163370,163371,163371,163371,163374,163374,163374,163374,163374,163374,163374,163376,163376,163376,163376,163378,163378,163378,163378,163378,163385,163385,163385,163387,163402,163402,163402,163402,163402,163402,163402,163402,163402,163475,163475,163475,163475,163475,163476,163476,163476,163479,163480,163482,163483,163583,163583,163583,163583,163586,163586,163588,163588,163588,163588,163588,163589,163589,163630,163657,163657,163659,163659,163659,163659,163659,163659,163659,163659,163659,163659,163659,163659,163659,163659,163659,163660,163660,163660,163661,163662,163662,163662,163664,163668,163668,163738,163738,163738,163738,163738,163738,163738,163758,163758,163758,163758,163758,163758,163758,163758,163758,163759,163759,163759,163759,163759,163760,163760,163760,163760,163760,163760,163760,163765,163765,163765,163765,163765,163795,163796,163796,163796,163796,163796,163796,163796,163817,163817,163817,163817,163817,163817,163817,163817,163817,163817,163817,163817,163817,163817,163817,163817,163817,163817,163817,165121,165121,165121,165121,165205,165208,165364,165364,165364,165364,165364,165364,165364,165364,165364,165364,165364,165364,165364,165364,165364,165364,165364,165364,165364,165364,165415,165415,165415,165415,165415,165415,165415,165415,165415,165415,165415,165415,165415,165415,165415,165415,165415,165415,165415,165415,165415,165415,165472,165472,165472,165472,165472,165472,165472,165472,165472,165472,165472,165472,165472,165472,165472,165472,165472,165486,165486,165487,165487,165487,165506,165506,165506],
Event = [24371951,18576870,96165671,491483557,91390314,125176731,48370161,65400245,119647477,107661306,265156983,25893970,56961763,38156147,48877189,6833339,4081021,63316951,17965714,183490182,220791530,228514066,28444848,33069964,95426489,101980751,239615748,304219363,400003395,55925035,60979586,133314497,200334760,205668936,215533223,255230619,259924803,282983988,293475492,305561898,342468785,348206909,394556047,399317868,418689254,50813748,489282865,505027103,93884489,97510916,108246660,161606268,162614748,104675948,189901294,333605114,52156625,22763918,56652154,4444230,57368030,5816503,107025374,138670755,191089040,219614010,228555655,5145145,86299370,89814612,92689113,15820394,18004679,22955850,48048119,182796843,229951639,290060361,310517723,334605324,453062228,21712494,116524826,169694364,182712860,263747194,269344428,61667010,73042788,111402754,130881024,146874823,4176096,19533169,41835038,4872222,76822473,40515999,4830513,60877125,16850617,200637251,29279251,7404160,89986633,48665006,161654590,166663172,52981357,67156409,177230370,242487385,349063851,375785727,378073345,381604121,80716296,23190414,35431342,35528644,56199361,58499014,126757726,166770991,1982546,136454635,213152538,41600744,25617613,65139186,113404298,186042589,294089105,348935564,359026121,402589263,35624459,479976713,83234226,88598161,150131089,193319870,44318680,110347959,3117371,57278411,4598210,5754166,4079087,9585045,126856305,165960555,21695548,41160080,10841564,36030982,137621860,162720793,168860940,185812235,293183271,55453063,63115578,76189444,12455810,5071047,99822371,104684639,110712786,185591256,194466360,218809820,24621131,256469811,27226115,300049345,375437326,380594303,399218272,419304899,496208715,20153888,2783995,40921027,693507,94830533,23177966,23950430,25050536,130748015,58368335,130926468,153349232,166379378,169183638,193537170,196931479,239901350,78306525,178588224,211524839,287243682,299874427,341969813,350570830,45063091,3354782,80917464,196972726,200893656,22548957,338699964,84895772,115506843,120288802,152683630,212944349,226711544,22637838,68673211,70541404,108574122,155024796,57070932,29194371,109508279,135236412,24000133,49641227,83624625,89854453,95897622,149299205,198315854,247061851,373312579,378084407,384623891,407058296,545248000,581734803,584317893,613490574,630579795,690340277,710876840,717817666,732789398,808113892,809391070,828622801,49191435,55002119,134938767,37485490,185695388,72542883,1184627074,1267860435,1358273514,1399958437,141284880,1477949824,1479075041,254428559,264185506,316591001,337725233,429993340,481026985,482153614,534253611,556171955,574817001,590583843,668793760,883140204,1291332300,1359273031,1410268779,113344188,1478044894,1485684156,1524657887,1573260343,251863772,450088973,479919143,504597781,516397079,542615397,556750519,667312880,708714244,751700331,39733972,774182048,816637188,955264961,133535097,173155335,14317087,18452039,305570274,436376807,468095459,483547880,717190088,841400997,879953987,885582991,897382077,902404476,941883006,957284610,98658826,35727889,79831820,136064858,76243949,77753094,141915838,146615579,50098607],
Lumi = [47,57,181,915,201,277,105,130,202,220,481,93,104,70,90,14,9,122,35,363,441,457,55,64,136,146,355,457,613,104,113,212,321,329,346,414,422,461,483,503,567,577,660,668,705,80,836,865,245,252,274,388,390,174,303,532,99,40,97,9,105,12,197,257,357,413,431,11,165,178,184,32,36,46,95,276,342,428,457,493,673,36,198,292,315,458,468,116,138,211,248,279,9,41,86,11,147,78,11,120,34,402,58,12,135,73,275,284,84,123,330,458,675,726,730,737,164,48,73,73,119,124,270,358,5,236,365,85,47,107,176,283,449,537,553,625,65,756,124,131,213,272,73,160,6,82,8,10,8,18,165,217,29,53,16,50,201,238,247,273,437,81,91,154,16,8,128,135,142,242,254,286,32,338,35,400,512,519,546,576,693,32,5,64,2,125,31,32,39,209,92,159,186,202,206,236,240,294,104,240,285,403,421,483,495,61,6,123,285,290,35,476,121,167,174,221,309,330,33,112,115,177,255,94,32,113,140,25,51,86,93,99,168,216,264,392,397,404,427,584,624,627,658,677,744,772,780,797,883,885,907,128,143,351,97,213,82,1007,1088,1169,1202,121,1267,1268,206,213,254,270,344,385,386,428,446,461,474,540,725,1094,1149,1191,119,1247,1253,1286,1328,224,384,409,430,439,462,474,568,604,644,64,664,702,832,115,149,14,18,268,386,416,430,620,717,747,752,761,765,797,809,85,41,73,106,61,62,150,153,78]
)




evFilter = OP_RunLumiEvSelector(Events.ps())


# -----------------------------------------------------------------------------
skim_ps=PSet(
    SkimName = "myskim",
    DropBranches = False,
    Branches = [
        " keep * "
        ]
)
skim = SkimOp(skim_ps.ps())


#Plot the common plots!

genericPSet = PSet(
DirName      = "275_325Gev",
MinObjects   = 0,
MaxObjects   = 15,
StandardPlots     = True,
)
pset1 = PSet(
DirName      = "275_325Gev",
MinObjects   = 1,
MaxObjects   = 2,
StandardPlots     = True,
)

Npset1 = PSet(
DirName      = "n275_325Gev",
MinObjects   = 3,
MaxObjects   = 15,
StandardPlots     = True,
)

pset2 = PSet(
DirName      = "325Gev",
MinObjects   = 2,
MaxObjects   = 2,
StandardPlots     = True,
)

Npset2 = PSet(
DirName      = "n325Gev",
MinObjects   = 3,
MaxObjects   = 15,
StandardPlots     = True,
)

pset3 = PSet(
DirName      = "375Gev",
MinObjects   = 2,
MaxObjects   = 2,
StandardPlots     = True,
)

Npset3 = PSet(
DirName      = "n375Gev",
MinObjects   = 3,
MaxObjects   = 15,
StandardPlots     = True,
)


pset5 = PSet(
DirName      = "375Gev_afterDeadEcal",
MinObjects   = 2,
MaxObjects   = 2,
StandardPlots     = True,
)

Npset5 = PSet(
DirName      = "n375Gev_afterDeadEcal",
MinObjects   = 3,
MaxObjects   = 15,
StandardPlots     = True,
)

pset4 = PSet(
DirName      = "All",
MinObjects   = 2,
MaxObjects   = 2,
StandardPlots     = True,
)

Npset4 = PSet(
DirName      = "nAll",
MinObjects   = 3,
MaxObjects   = 15,
StandardPlots     = True,
)

Npset6 = PSet(
DirName      = "nAllCuts",
MinObjects   = 3,
MaxObjects   = 15,
StandardPlots     = True,
)

pset6 = PSet(
DirName      = "AllCuts",
MinObjects   = 2,
MaxObjects   = 2,
StandardPlots     = True,
)

dalitz_plots_Inclusive = HadronicPlottingOps( PSet(
DirName    = "Dalitz_Inclusive",
MinObjects = 2,
MaxObjects = 10,
Verbose    = False,
Summary    = False,
CC         = False,
Dalitz     = True,
AlphaT     = False,
PtHat      = False,
MET        = False,
Kine       = False,
Response   = False,
).ps()
)

dalitz_plots_275_325 = HadronicPlottingOps( PSet(
DirName    = "Dalitz_275_325",
MinObjects = 2,
MaxObjects = 10,
Verbose    = False,
Summary    = False,
CC         = False,
Dalitz     = True,
AlphaT     = False,
PtHat      = False,
MET        = False,
Kine       = False,
Response   = False,
).ps()
)

dalitz_plots_325_375 = HadronicPlottingOps( PSet(
DirName    = "Dalitz_325_375",
MinObjects = 2,
MaxObjects = 10,
Verbose    = False,
Summary    = False,
CC         = False,
Dalitz     = True,
AlphaT     = False,
PtHat      = False,
MET        = False,
Kine       = False,
Response   = False,
).ps()
)

dalitz_plots_375 = HadronicPlottingOps( PSet(
DirName    = "Dalitz_375",
MinObjects = 2,
MaxObjects = 10,
Verbose    = False,
Summary    = False,
CC         = False,
Dalitz     = True,
AlphaT     = False,
PtHat      = False,
MET        = False,
Kine       = False,
Response   = False,
).ps()
)




HadStandard275_375 = WeeklyUpdatePlots(pset1.ps())
HadStandard325 = WeeklyUpdatePlots(pset2.ps())
HadStandard375 = WeeklyUpdatePlots(pset3.ps())
HadStandard375_after_DeadEcal = WeeklyUpdatePlots(pset5.ps())
HadStandardAll = WeeklyUpdatePlots(pset4.ps())
nHadStandard275_375 = WeeklyUpdatePlots(Npset1.ps())
nHadStandard325 = WeeklyUpdatePlots(Npset2.ps())
nHadStandard375 = WeeklyUpdatePlots(Npset3.ps())
nHadStandardAll = WeeklyUpdatePlots(Npset4.ps())
nHadStandard375_after_DeadEcal = WeeklyUpdatePlots(Npset5.ps())
# Common cut definitions
#Avaiable criteria for MC and for Data are at current slightly different Hence the making of two trees
#DataOnly!

# from icf.JetCorrections import *
# corPset =  CorrectionPset("ResidualJetEnergyCorrections.txt")
# corPset =  CorrectionPset("Spring10DataV2_L2L3Residual_AK5PF.txt")
# JetCorrections = JESCorrections( corPset.ps(),True )
NoiseFilt= OP_HadronicHBHEnoiseFilter()
GoodVertexMonster = OP_GoodEventSelection()

#Standard Event Selection
LeadingJetEta = OP_FirstJetEta(2.5)
unCorLeadJetCut = OP_UnCorLeadJetCut(30.)
LeadingJetPtCut = OP_FirstJetPtCut(100.)
oddMuon = OP_OddMuon()
oddElectron = OP_OddElectron()
oddPhoton = OP_OddPhoton()
oddJet = OP_OddJet()
badMuonInJet = OP_BadMuonInJet()
numComLeptons = OP_NumComLeptons("<=",0)
numComPhotons = OP_NumComPhotons("<=",0)

DiJet0 = OP_NumComJets("==",2)
DiJet1 = OP_NumComJets("==",2)
DiJet2 = OP_NumComJets("==",2)
DiJet3 = OP_NumComJets("==",2)
DiJet4 = OP_NumComJets("==",2)
NJet0 = OP_NumComJets(">=",3)
NJet1 = OP_NumComJets(">=",3)
NJet2 = OP_NumComJets(">=",3)
NJet3 = OP_NumComJets(">=",3)
NJet4 = OP_NumComJets(">=",3)
DiVertexJets = OP_NumComJets("==",2)
NVertexJets = OP_NumComJets(">=",3)



LessThan375 = RECO_CommonHTLessThanCut(375.)
ht250_Trigger = RECO_CommonHTCut(250.)
htCut275 = RECO_CommonHTCut(275.)
htCut325 = RECO_CommonHTCut(325.)
htCut375 = RECO_CommonHTCut(375.)
htCut275_2 = RECO_CommonHTCut(275.)

ht275_Fail     = RECO_CommonHTCut(275.)
ht325_Fail     = RECO_CommonHTCut(325.)
ht375_Fail     = RECO_CommonHTCut(375.)
htLess325_Fail = RECO_CommonHTLessThanCut(325.)
htLess375_Fail = RECO_CommonHTLessThanCut(375.)


htCut375GeV = RECO_CommonHTCut(375.)
htCut375All = RECO_CommonHTCut(375.)
ht250_Trigger = RECO_CommonHTCut(250.)
ht275     = RECO_CommonHTCut(275.)
ht325     = RECO_CommonHTCut(325.)
ht375     = RECO_CommonHTCut(375.)
ht475     = RECO_CommonHTCut(475.)
ht575     = RECO_CommonHTCut(575.)
ht675     = RECO_CommonHTCut(675.)
ht775     = RECO_CommonHTCut(775.)
ht875     = RECO_CommonHTCut(875.)
htLess325 = RECO_CommonHTLessThanCut(325.)
htLess375 = RECO_CommonHTLessThanCut(375.)
htLess475 = RECO_CommonHTLessThanCut(475.)
htLess575 = RECO_CommonHTLessThanCut(575.)
htLess675 = RECO_CommonHTLessThanCut(675.)
htLess775 = RECO_CommonHTLessThanCut(775.)
htLess875 = RECO_CommonHTLessThanCut(875.)
htCut275_2 = RECO_CommonHTCut(275.)
htCut375GeV = RECO_CommonHTCut(375.)
alphaT0 = OP_CommonAlphaTCut(0.55)
alphaT1 = OP_CommonAlphaTCut(0.55)
alphaT2 = OP_CommonAlphaTCut(0.55)
spikecleaner = OP_EcalSpikeCleaner()
event_display = OP_EventDisplay("EventDisplays", "common") #to draw all/common objects
alphat = OP_CommonAlphaTCut(0.55)
DeadEcalCutData = OP_DeadECALCut(0.3,0.3,0.5,30.,10,0,"./deadRegionList_GR10_P_V10.txt")
DeadEcalCutMC =   OP_DeadECALCut(0.3,0.3,0.5,30.,10,0,"./deadRegionList_START38_V12.txt")
MHTCut = OP_TriggerMHT_Emu(60.,30.)
MHT_METCut = OP_MHToverMET(1.25,50.)
NJet5 = OP_NumComJets(">=",3)
DiJet5 = OP_NumComJets("==",2)
nHadStandardAllCuts=  WeeklyUpdatePlots(Npset6.ps())
HadStandardAllCuts=  WeeklyUpdatePlots(pset6.ps())


# Cross check with the allhadronic analysis
t1 = PSet(
    DirName      = "HadronicCommon_1",
    MinObjects   = 2,
    MaxObjects   = 15,
    StandardPlots     = False,
    DeadECALPlots = True,
    CleaningControlPlots = False,
    MECPlots = False,
    DeadECALFile = "./deadRegionList_START36_V9.txt",
    DeadECAL_MinJetPtCut = 30.,
    DeadECAL_MinBiasCut = 0.5,
    DeadECAL_NBadCellsCut = 10
)

t2 = deepcopy(t1)
t2.DirName = "HadronicCommon_2"

pset = PSet(
DirName      = "275_infGev",
MinObjects   = 1,
MaxObjects   = 2,
StandardPlots     = True,
)

Npset = PSet(
DirName      = "n275_infGev",
MinObjects   = 3,
MaxObjects   = 15,
StandardPlots     = True,
)

pset2 = deepcopy(pset)
pset2.DirName = "275_375Gev"

Npset2 = deepcopy(Npset)
Npset2.DirName = "n275_375Gev"

pset3 = deepcopy(pset)
pset3.DirName = "375GeVafterDeadEcal"
Npset3 = deepcopy(Npset)
Npset3.DirName = "n375GeVafterDeadEcal"

pset4 = deepcopy(pset)
pset4.DirName = "allCuts"
Npset4 = deepcopy(Npset)
Npset4.DirName = "nAllCuts"
# Define a crap load more plotting ops, for HT exclusive bins
Plot_275_325_pset = deepcopy(genericPSet)
Plot_275_325_pset.DirName="275_325"
Plot_325_375_pset = deepcopy(genericPSet)
Plot_325_375_pset.DirName="325_375"
Plot_375_475_pset = deepcopy(genericPSet)
Plot_375_475_pset.DirName="375_475"
Plot_475_575_pset = deepcopy(genericPSet)
Plot_475_575_pset.DirName="475_575"
Plot_575_675_pset = deepcopy(genericPSet)
Plot_575_675_pset.DirName="575_675"
Plot_675_775_pset = deepcopy(genericPSet)
Plot_675_775_pset.DirName="675_775"
Plot_775_875_pset = deepcopy(genericPSet)
Plot_775_875_pset.DirName="775_875"
Plot_875_pset = deepcopy(genericPSet)
Plot_875_pset.DirName="875"

Plot_275_325_Fail_pset = deepcopy(genericPSet)
Plot_275_325_Fail_pset.DirName="275_325Fail"
Plot_325_375_Fail_pset = deepcopy(genericPSet)
Plot_325_375_Fail_pset.DirName="325_375Fail"
Plot_375_Fail_pset = deepcopy(genericPSet)
Plot_375_Fail_pset.DirName="375Fail"
Plot_275_325Fail = WeeklyUpdatePlots( Plot_275_325_Fail_pset.ps() )
Plot_325_375Fail = WeeklyUpdatePlots( Plot_325_375_Fail_pset.ps() )
Plot_375Fail     = WeeklyUpdatePlots( Plot_375_Fail_pset.ps() )




Plot_275_325 = WeeklyUpdatePlots( Plot_275_325_pset.ps() )
Plot_325_375 = WeeklyUpdatePlots( Plot_325_375_pset.ps() )
Plot_375_475 = WeeklyUpdatePlots( Plot_375_475_pset.ps() )
Plot_475_575 = WeeklyUpdatePlots( Plot_475_575_pset.ps() )
Plot_575_675 = WeeklyUpdatePlots( Plot_575_675_pset.ps() )
Plot_675_775 = WeeklyUpdatePlots( Plot_675_775_pset.ps() )
Plot_775_875 = WeeklyUpdatePlots( Plot_775_875_pset.ps() )
Plot_875     = WeeklyUpdatePlots( Plot_875_pset.ps()     )
HTplot = WeeklyUpdatePlots(pset.ps())
nHTplot = WeeklyUpdatePlots(Npset.ps())
controlRegion = WeeklyUpdatePlots(pset2.ps())
ncontrolRegion = WeeklyUpdatePlots(Npset2.ps())
afterDeadEcal = WeeklyUpdatePlots(pset3.ps())
nafterDeadEcal = WeeklyUpdatePlots(Npset3.ps())
afterAllCuts = WeeklyUpdatePlots(pset4.ps())
nafterAllCuts = WeeklyUpdatePlots(Npset4.ps())


pset2 = deepcopy(pset1)
pset2.DirName = "HadronicCommon_2"

t3 = deepcopy(t1)
t3.DirName = "HadronicCommon_3"

t4 = deepcopy(t1)
t4.DirName = "HadronicCommon_4"
#
HadStandard_1 = HadronicCommonPlots(t1.ps())
HadStandard_2 = HadronicCommonPlots(t2.ps())
HadStandard_3 = HadronicCommonPlots(t3.ps())
HadStandard_4 = HadronicCommonPlots(t4.ps())
VertexPtOverHT = OP_SumVertexPtOverHT(0.1)
# eventDump = OP_EventNoDump("mydump","mydump")
eventDump = EventDump()
datatriggerps = PSet(
    Verbose = False,
    Triggers = [
        # "HLT_HT150_v1",
        # "HLT_HT150_v2",
        # "HLT_HT150_v3",
        # "HLT_HT150_v4",
        # "HLT_HT150_v5",
        # "HLT_HT150_v6",
        # "HLT_HT150_v7",
        # "HLT_HT150_v8",
        # "HLT_HT200_v1",
        # "HLT_HT200_v2",
        # "HLT_HT200_v3",
        # "HLT_HT200_v4",
        # "HLT_HT200_v5",
        # "HLT_HT200_v6",
        # "HLT_HT200_v7",
        # "HLT_HT200_v8",
        "HLT_HT250_v1",
        "HLT_HT250_v2",
        "HLT_HT250_v3",
        "HLT_HT250_v4",
        "HLT_HT250_v5",
        "HLT_HT250_v6",
        "HLT_HT250_v7",
        "HLT_HT250_v8",
        "HLT_HT300_v1",
        "HLT_HT300_v2",
        "HLT_HT300_v3",
        "HLT_HT300_v4",
        "HLT_HT300_v5",
        "HLT_HT300_v6",
        "HLT_HT300_v7",
        "HLT_HT300_v8",
        "HLT_HT350_v1",
        "HLT_HT350_v2",
        "HLT_HT350_v3",
        "HLT_HT350_v4",
        "HLT_HT350_v5",
        "HLT_HT350_v6",
        "HLT_HT350_v7",
        "HLT_HT350_v8",
        "HLT_HT400_v1",
        "HLT_HT400_v2",
        "HLT_HT400_v3",
        "HLT_HT400_v4",
        "HLT_HT400_v5",
        "HLT_HT400_v6",
        "HLT_HT400_v7",
        "HLT_HT400_v8",
        "HLT_HT450_v1",
        "HLT_HT450_v2",
        "HLT_HT450_v3",
        "HLT_HT450_v4",
        "HLT_HT450_v5",
        "HLT_HT450_v6",
        "HLT_HT450_v7",
        "HLT_HT450_v8",
        "HLT_HT500_v1",
        "HLT_HT500_v2",
        "HLT_HT500_v3",
        "HLT_HT500_v4",
        "HLT_HT500_v5",
        "HLT_HT500_v6",
        "HLT_HT500_v7",
        "HLT_HT500_v8",
       "HLT_HT260_MHT60_v1",
       "HLT_HT260_MHT60_v2",
       "HLT_HT250_MHT60_v1",
       "HLT_HT250_MHT60_v2",
       "HLT_HT250_MHT50_v2",
       "HLT_HT250_MHT60_v3",
       "HLT_HT250_MHT70_v1",
       "HLT_HT250_MHT70_v2",
       "HLT_HT250_MHT70_v3",
       "HLT_HT250_MHT60_v3",
       "HLT_HT250_AlphaT0p53_v1",
       "HLT_HT250_AlphaT0p53_v2",
       "HLT_HT250_AlphaT0p53_v3",
        ]
    )
DataTrigger = OP_MultiTrigger( datatriggerps.ps() )

JetAdd = JetAddition(0.)
json = JSONFilter("Json Mask", json_to_pset("485pbjson.txt"))

# AlphatTriggerCut(0.52414,50)#
vertex_reweight = VertexReweighting(
PSet(
VertexWeights =[0.20, 0.63, 1.19, 1.57, 1.62, 1.42, 1.09, 0.80 ,0.57, 0.42, 0.30, 0.20]
# VertexWeights = [0.0, 0.027442995662725636, 0.12682983875287387, 0.28326829632076572, 0.40618954180036759, 0.41605144586432974, 0.33147399297403923, 0.21562021576661147, 0.1140047132529971]
).ps())

PreScaleWeights = PreScaleReweighting(datatriggerps.ps())

json_ouput = JSONOutput("filtered")
def MakeDataTree(Threshold):
  secondJetET = OP_SecondJetEtCut(Threshold)
  # from batchGolden import *
  cutTreeData = Tree("Data")
  cutTreeData.Attach(json)
  # cutTreeData.TAttach(json,evFilter)
  # cutTreeData.TAttach(evFilter,skim)
  # cutTreeData.TAttach(evFilter,eventDump)
  cutTreeData.TAttach(json,DataTrigger)
  cutTreeData.TAttach(json,json_ouput)
  #cutTreeData.Attach(DataTrigger)
  cutTreeData.TAttach(DataTrigger,NoiseFilt)
  cutTreeData.TAttach(NoiseFilt,GoodVertexMonster)
  cutTreeData.TAttach(GoodVertexMonster,LeadingJetEta)
  cutTreeData.TAttach(LeadingJetEta,secondJetET)
  cutTreeData.TAttach(secondJetET,oddJet)
  cutTreeData.TAttach(oddJet,badMuonInJet)
  cutTreeData.TAttach(badMuonInJet,oddMuon)
  cutTreeData.TAttach(oddMuon,oddElectron)
  cutTreeData.TAttach(oddElectron,oddPhoton)
  cutTreeData.TAttach(oddPhoton,numComLeptons)
  cutTreeData.TAttach(numComLeptons,numComPhotons)
  cutTreeData.TAttach(numComPhotons,VertexPtOverHT)
  cutTreeData.TAttach(VertexPtOverHT,htCut275)
  cutTreeData.TAttach(numComPhotons,ht275_Fail)
  cutTreeData.TAttach(numComPhotons,ht325_Fail)
  cutTreeData.TAttach(ht275_Fail,htLess325_Fail)
  cutTreeData.TAttach(ht325_Fail,htLess375_Fail)
  cutTreeData.TAttach(numComPhotons,ht375_Fail)
  cutTreeData.TAttach(htLess325_Fail,Plot_275_325Fail)
  cutTreeData.TAttach(htLess375_Fail,Plot_325_375Fail)
  cutTreeData.TAttach(ht375_Fail,Plot_375Fail)
  #cutTreeData.TAttach(htCut275,skim)
  #FOR HT > 275Gev Plot
  cutTreeData.TAttach(htCut275,DiJet3)
  cutTreeData.TAttach(htCut275,NJet3)
  cutTreeData.TAttach(DiJet3,HadStandardAll)
  cutTreeData.TAttach(NJet3,nHadStandardAll)
  #END HT 275GEV Plot
  #Begin MHT/MET plot inthe low region.
  cutTreeData.TAttach(htCut275,DeadEcalCutData)
  cutTreeData.TAttach(DeadEcalCutData,LessThan375)
  #cutTreeData.TAttach(DeadEcalCutData,skim)
  cutTreeData.TAttach(LessThan375,DiJet0)
  cutTreeData.TAttach(LessThan375,NJet0)
  cutTreeData.TAttach(DiJet0,HadStandard275_375)
  cutTreeData.TAttach(NJet0,nHadStandard275_375)
  #for Vertext multiplicity plot at 325geV
  # cutTreeData.TAttach(htCut275,htCut325)
  # cutTreeData.TAttach(htCut325,NVertexJets)
  # cutTreeData.TAttach(htCut325,DiVertexJets)
  # cutTreeData.TAttach(DiVertexJets,HadStandard325)
  # cutTreeData.TAttach(NVertexJets,nHadStandard325)
  # cutTreeData.TAttach(htCut375,dalitz_plots_375)
  cutTreeData.TAttach(htCut275,htCut375GeV)
  cutTreeData.TAttach(htCut375GeV,DiJet2)
  cutTreeData.TAttach(htCut375GeV,NJet2)
  # cutTreeData.TAttach(htCut375,alphaT0)
  cutTreeData.TAttach(DiJet2,HadStandard375)
  cutTreeData.TAttach(NJet2,nHadStandard375)
  cutTreeData.TAttach(DeadEcalCutData,htCut375)
  #Here be plots for baby jet MHT and MHT/MET in the signal region after dead ecal cuts
  cutTreeData.TAttach(htCut375,DiJet4)
  cutTreeData.TAttach(DiJet4,HadStandard375_after_DeadEcal)
  cutTreeData.TAttach(htCut375,NJet4)
  cutTreeData.TAttach(htCut375,alphaT1)
  # cutTreeData.TAttach(alphaT1,HadStandard_2)
  cutTreeData.TAttach(NJet4,nHadStandard375_after_DeadEcal)
  #Here be plots after all the cuts!!
  # cutTreeData.TAttach(htCut375GeV,alphaT2)
  # cutTreeData.TAttach(htCut375,MHT_METCut)
  cutTreeData.TAttach(DeadEcalCutData,MHT_METCut)
  # cutTreeData.TAttach(alphaT2,HadStandard_3)
  cutTreeData.TAttach(MHT_METCut,htCut375All)
  cutTreeData.TAttach(htCut375All,NJet5)
  cutTreeData.TAttach(htCut375All,DiJet5)
  cutTreeData.TAttach(htCut375All,alphat)
  cutTreeData.TAttach(alphat,eventDump)
  cutTreeData.TAttach(alphat,skim)
  cutTreeData.TAttach(NJet5,nHadStandardAllCuts)
  cutTreeData.TAttach(DiJet5,HadStandardAllCuts)
  # cutTreeData.TAttach(MHT_METCut, alphat)
  # cutTreeData.TAttach(alphat,eventDump)#skim)
  # avobe here does one big inclusive bin!
  # Now lets start binning in HT bins
  cutTreeData.TAttach(MHT_METCut,ht275)
  cutTreeData.TAttach(MHT_METCut,ht325)
  cutTreeData.TAttach(MHT_METCut,ht375)
  cutTreeData.TAttach(MHT_METCut,ht475)
  cutTreeData.TAttach(MHT_METCut,ht575)
  cutTreeData.TAttach(MHT_METCut,ht675)
  cutTreeData.TAttach(MHT_METCut,ht775)
  cutTreeData.TAttach(MHT_METCut,ht875)
  cutTreeData.TAttach(ht275,htLess325)
  cutTreeData.TAttach(ht325,htLess375)
  cutTreeData.TAttach(ht375,htLess475)
  cutTreeData.TAttach(ht475,htLess575)
  cutTreeData.TAttach(ht575,htLess675)
  cutTreeData.TAttach(ht675,htLess775)
  cutTreeData.TAttach(ht775,htLess875)
  cutTreeData.TAttach(htLess325,Plot_275_325)
  cutTreeData.TAttach(htLess375,Plot_325_375)
  cutTreeData.TAttach(htLess475,Plot_375_475)
  cutTreeData.TAttach(htLess575,Plot_475_575)
  cutTreeData.TAttach(htLess675,Plot_575_675)
  cutTreeData.TAttach(htLess775,Plot_675_775)
  cutTreeData.TAttach(htLess875,Plot_775_875)
  cutTreeData.TAttach(ht875,    Plot_875    )
  return (cutTreeData,secondJetET)

#Second MC!

def MakeMCTree(Threshold):
  secondJetET = OP_SecondJetEtCut(Threshold)
  cutTreeMC = Tree("MC")
  cutTreeMC.Attach(ht250_Trigger)
# ,MHTCut)
#   cutTreeMC.TAttach(MHTCut

  cutTreeMC.TAttach(ht250_Trigger,NoiseFilt)
  cutTreeMC.TAttach(NoiseFilt,GoodVertexMonster)
  cutTreeMC.TAttach(GoodVertexMonster,LeadingJetEta)
  cutTreeMC.TAttach(LeadingJetEta,secondJetET)
  cutTreeMC.TAttach(secondJetET,oddJet)
  cutTreeMC.TAttach(oddJet,badMuonInJet)
  cutTreeMC.TAttach(badMuonInJet,oddMuon)
  cutTreeMC.TAttach(oddMuon,oddElectron)
  cutTreeMC.TAttach(oddElectron,oddPhoton)
  cutTreeMC.TAttach(oddPhoton,numComLeptons)
  cutTreeMC.TAttach(numComLeptons,numComPhotons)
  cutTreeMC.TAttach(numComPhotons,VertexPtOverHT)
  cutTreeMC.TAttach(VertexPtOverHT,htCut275)
  cutTreeMC.TAttach(numComPhotons,ht275_Fail)
  cutTreeMC.TAttach(numComPhotons,ht325_Fail)
  cutTreeMC.TAttach(ht275_Fail,htLess325_Fail)
  cutTreeMC.TAttach(ht325_Fail,htLess375_Fail)
  cutTreeMC.TAttach(numComPhotons,ht375_Fail)
  cutTreeMC.TAttach(htLess325_Fail,Plot_275_325Fail)
  cutTreeMC.TAttach(htLess375_Fail,Plot_325_375Fail)
  cutTreeMC.TAttach(ht375_Fail,Plot_375Fail)
  #FOR HT > 275Gev Plot
  cutTreeMC.TAttach(htCut275,DiJet3)
  cutTreeMC.TAttach(htCut275,NJet3)
  cutTreeMC.TAttach(DiJet3,HadStandardAll)
  cutTreeMC.TAttach(NJet3,nHadStandardAll)
  #END HT 275GEV Plot
  #Begin MHT/MET plot inthe low region.
  cutTreeMC.TAttach(htCut275,DeadEcalCutMC)
  cutTreeMC.TAttach(DeadEcalCutMC,LessThan375)
  cutTreeMC.TAttach(LessThan375,DiJet0)
  cutTreeMC.TAttach(LessThan375,NJet0)
  cutTreeMC.TAttach(DiJet0,HadStandard275_375)
  cutTreeMC.TAttach(NJet0,nHadStandard275_375)

  #for Vertext multiplicity plot at 325geV
  # cutTreeMC.TAttach(htCut275,htCut325)
  # cutTreeMC.TAttach(htCut325,NVertexJets)
  # cutTreeMC.TAttach(htCut325,DiVertexJets)
  # cutTreeMC.TAttach(DiVertexJets,HadStandard325)
  # cutTreeMC.TAttach(NVertexJets,nHadStandard325)


  cutTreeMC.TAttach(DeadEcalCutMC,htCut375)
  # cutTreeMC.TAttach(htCut375,dalitz_plots_375)
  cutTreeMC.TAttach(htCut275,htCut375GeV)
  cutTreeMC.TAttach(htCut375GeV,DiJet2)
  cutTreeMC.TAttach(htCut375GeV,NJet2)
  # cutTreeMC.TAttach(htCut375,alphaT0)
  cutTreeMC.TAttach(DiJet2,HadStandard375)
  cutTreeMC.TAttach(NJet2,nHadStandard375)

  #Here be plots for baby jet MHT and MHT/MET in the signal region after dead ecal cuts
  cutTreeMC.TAttach(htCut375,DiJet4)
  cutTreeMC.TAttach(DiJet4,HadStandard375_after_DeadEcal)
  cutTreeMC.TAttach(htCut375,NJet4)
  cutTreeMC.TAttach(htCut375,alphaT1)
  # cutTreeMC.TAttach(alphaT1,HadStandard_2)
  cutTreeMC.TAttach(NJet4,nHadStandard375_after_DeadEcal)


  #Here be plots after all the cuts!!
  cutTreeMC.TAttach(DeadEcalCutMC,MHT_METCut)
  # cutTreeMC.TAttach(htCut375,MHT_METCut)
  # cutTreeMC.TAttach(alphaT2,HadStandard_3)
  cutTreeMC.TAttach(MHT_METCut,alphaT2)
  cutTreeMC.TAttach(MHT_METCut,htCut375All)
  cutTreeMC.TAttach(htCut375All,NJet5)
  cutTreeMC.TAttach(htCut375All,DiJet5)
  cutTreeMC.TAttach(NJet5,nHadStandardAllCuts)
  cutTreeMC.TAttach(DiJet5,HadStandardAllCuts)


  cutTreeMC.TAttach(MHT_METCut,ht275)
  cutTreeMC.TAttach(MHT_METCut,ht325)
  cutTreeMC.TAttach(MHT_METCut,ht375)
  cutTreeMC.TAttach(MHT_METCut,ht475)
  cutTreeMC.TAttach(MHT_METCut,ht575)
  cutTreeMC.TAttach(MHT_METCut,ht675)
  cutTreeMC.TAttach(MHT_METCut,ht775)
  cutTreeMC.TAttach(MHT_METCut,ht875)
  cutTreeMC.TAttach(ht275,htLess325)
  cutTreeMC.TAttach(ht325,htLess375)
  cutTreeMC.TAttach(ht375,htLess475)
  cutTreeMC.TAttach(ht475,htLess575)
  cutTreeMC.TAttach(ht575,htLess675)
  cutTreeMC.TAttach(ht675,htLess775)
  cutTreeMC.TAttach(ht775,htLess875)
  cutTreeMC.TAttach(htLess325,Plot_275_325)
  cutTreeMC.TAttach(htLess375,Plot_325_375)
  cutTreeMC.TAttach(htLess475,Plot_375_475)
  cutTreeMC.TAttach(htLess575,Plot_475_575)
  cutTreeMC.TAttach(htLess675,Plot_575_675)
  cutTreeMC.TAttach(htLess775,Plot_675_775)
  cutTreeMC.TAttach(htLess875,Plot_775_875)
  cutTreeMC.TAttach(ht875,    Plot_875    )
  return (cutTreeMC,secondJetET)

