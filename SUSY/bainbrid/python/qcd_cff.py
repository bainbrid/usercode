from icf.core import PSet

# Pythia 6

qcd_path = "/vols/cms02/gouskos/"

qcd_pythia6 = PSet(
    Name = "QCD",
    Format = ("ICF",2),
    File = (
    qcd_path + "QCD_Pythia_Pt470_Jun2010/SusyCAF_Tree_1_1_9gm.root",
    #qcd_path + "QCD_Pythia_Pt15_Jun2010/SusyCAF_Tree_*.root",
    #qcd_path + "QCD_Pythia_Pt30_Jun2010/SusyCAF_Tree_*.root",
    #qcd_path + "QCD_Pythia_Pt80_Jun2010/SusyCAF_Tree_*.root",
    #qcd_path + "QCD_Pythia_Pt170_Jun2010/SusyCAF_Tree_*.root",
    #qcd_path + "QCD_Pythia_Pt300_Jun2010/SusyCAF_Tree_*.root",
    #qcd_path + "QCD_Pythia_Pt470_Jun2010/SusyCAF_Tree_*.root",
    #qcd_path + "QCD_Pythia_Pt800_Jun2010/SusyCAF_Tree_*.root",
    #qcd_path + "QCD_Pythia_Pt1400_Jun2010/SusyCAF_Tree_*.root",
    ),
    Weights = PSet(
    CrossSection = [ 8.762e+08, 6.041e+07, 9.238e+05, 2.547e+04, 1.256e+03, 8.798e+01, 2.186, 0.01122 ],
    Events       = [ 6095857, 5069664, 2065792, 3171950, 2976108, 2159497, 2181700, 1185024 ],
    PtBin        = [ 15., 30., 80., 170., 300., 470., 800., 1400. ],
    )
    )

# Pythia 8

from montecarlo.QCD_Pt_0to15_7TeV_pythia8_Summer10_START36_V10_S09_v1 import *
from montecarlo.QCD_Pt_1000to1400_7TeV_pythia8_Summer10_START36_V10_S09_v2 import *
from montecarlo.QCD_Pt_120to170_7TeV_pythia8_Summer10_START36_V10_S09_v1 import *
from montecarlo.QCD_Pt_1400to1800_7TeV_pythia8_Summer10_START36_V10_S09_v2 import *
from montecarlo.QCD_Pt_15to20_7TeV_pythia8_Summer10_START36_V10_S09_v1 import *
from montecarlo.QCD_Pt_170to230_7TeV_pythia8_Summer10_START36_V10_S09_v2 import *
from montecarlo.QCD_Pt_1800to2200_7TeV_pythia8_Summer10_START36_V10_S09_v2 import *
from montecarlo.QCD_Pt_20to30_7TeV_pythia8_Summer10_START36_V10_S09_v1 import *
from montecarlo.QCD_Pt_2200to2600_7TeV_pythia8_Summer10_START36_V10_S09_v2 import *
from montecarlo.QCD_Pt_230to300_7TeV_pythia8_Summer10_START36_V10_S09_v2 import *
from montecarlo.QCD_Pt_2600to3000_7TeV_pythia8_Summer10_START36_V10_S09_v2 import *
from montecarlo.QCD_Pt_3000to3500_7TeV_pythia8_Summer10_START36_V10_S09_v2 import *
from montecarlo.QCD_Pt_300to380_7TeV_pythia8_Summer10_START36_V10_S09_v1 import *
from montecarlo.QCD_Pt_30to50_7TeV_pythia8_Summer10_START36_V10_S09_v2 import *
from montecarlo.QCD_Pt_380to470_7TeV_pythia8_Summer10_START36_V10_S09_v1 import *
from montecarlo.QCD_Pt_470to600_7TeV_pythia8_Summer10_START36_V10_S09_v1 import *
from montecarlo.QCD_Pt_50to80_7TeV_pythia8_Summer10_START36_V10_S09_v1 import *
from montecarlo.QCD_Pt_600to800_7TeV_pythia8_Summer10_START36_V10_S09_v1 import *
from montecarlo.QCD_Pt_800to1000_7TeV_pythia8_Summer10_START36_V10_S09_v2 import *
from montecarlo.QCD_Pt_80to120_7TeV_pythia8_Summer10_START36_V10_S09_v1 import *

qcd_pythia8_merged = [
    QCD_Pt_0to15_7TeV_pythia8_Summer10_START36_V10_S09_v1,
    QCD_Pt_15to20_7TeV_pythia8_Summer10_START36_V10_S09_v1,
    QCD_Pt_20to30_7TeV_pythia8_Summer10_START36_V10_S09_v1,
    QCD_Pt_30to50_7TeV_pythia8_Summer10_START36_V10_S09_v2,
    QCD_Pt_50to80_7TeV_pythia8_Summer10_START36_V10_S09_v1,
    QCD_Pt_80to120_7TeV_pythia8_Summer10_START36_V10_S09_v1,
    QCD_Pt_120to170_7TeV_pythia8_Summer10_START36_V10_S09_v1, 
    QCD_Pt_170to230_7TeV_pythia8_Summer10_START36_V10_S09_v2,
    QCD_Pt_230to300_7TeV_pythia8_Summer10_START36_V10_S09_v2,
    QCD_Pt_300to380_7TeV_pythia8_Summer10_START36_V10_S09_v1,
    QCD_Pt_380to470_7TeV_pythia8_Summer10_START36_V10_S09_v1,
    QCD_Pt_470to600_7TeV_pythia8_Summer10_START36_V10_S09_v1,
    QCD_Pt_600to800_7TeV_pythia8_Summer10_START36_V10_S09_v1,
    QCD_Pt_800to1000_7TeV_pythia8_Summer10_START36_V10_S09_v2,
    QCD_Pt_1000to1400_7TeV_pythia8_Summer10_START36_V10_S09_v2,
    QCD_Pt_1400to1800_7TeV_pythia8_Summer10_START36_V10_S09_v2,
    QCD_Pt_1800to2200_7TeV_pythia8_Summer10_START36_V10_S09_v2,
    QCD_Pt_2200to2600_7TeV_pythia8_Summer10_START36_V10_S09_v2,
    QCD_Pt_2600to3000_7TeV_pythia8_Summer10_START36_V10_S09_v2,
    QCD_Pt_3000to3500_7TeV_pythia8_Summer10_START36_V10_S09_v2,
    ]

qcd_files = (
    QCD_Pt_0to15_7TeV_pythia8_Summer10_START36_V10_S09_v1.File +
    QCD_Pt_15to20_7TeV_pythia8_Summer10_START36_V10_S09_v1.File +
    QCD_Pt_20to30_7TeV_pythia8_Summer10_START36_V10_S09_v1.File +
    QCD_Pt_30to50_7TeV_pythia8_Summer10_START36_V10_S09_v2.File +
    QCD_Pt_50to80_7TeV_pythia8_Summer10_START36_V10_S09_v1.File +
    QCD_Pt_80to120_7TeV_pythia8_Summer10_START36_V10_S09_v1.File +
    QCD_Pt_120to170_7TeV_pythia8_Summer10_START36_V10_S09_v1.File + 
    QCD_Pt_170to230_7TeV_pythia8_Summer10_START36_V10_S09_v2.File +
    QCD_Pt_230to300_7TeV_pythia8_Summer10_START36_V10_S09_v2.File +
    QCD_Pt_300to380_7TeV_pythia8_Summer10_START36_V10_S09_v1.File +
    QCD_Pt_380to470_7TeV_pythia8_Summer10_START36_V10_S09_v1.File +
    QCD_Pt_470to600_7TeV_pythia8_Summer10_START36_V10_S09_v1.File +
    QCD_Pt_600to800_7TeV_pythia8_Summer10_START36_V10_S09_v1.File +
    QCD_Pt_800to1000_7TeV_pythia8_Summer10_START36_V10_S09_v2.File +
    QCD_Pt_1000to1400_7TeV_pythia8_Summer10_START36_V10_S09_v2.File +
    QCD_Pt_1400to1800_7TeV_pythia8_Summer10_START36_V10_S09_v2.File +
    QCD_Pt_1800to2200_7TeV_pythia8_Summer10_START36_V10_S09_v2.File +
    QCD_Pt_2200to2600_7TeV_pythia8_Summer10_START36_V10_S09_v2.File +
    QCD_Pt_2600to3000_7TeV_pythia8_Summer10_START36_V10_S09_v2.File +
    QCD_Pt_3000to3500_7TeV_pythia8_Summer10_START36_V10_S09_v2.File 
    )

qcd_pythia8 = PSet(
    Name = "QCD_Pythia8",
    Format = ("ICF",2),
    File = qcd_files,
    Weights = PSet(
    CrossSection = [ 1. ],
    Events       = [ 770000, 520650, 613650, 651444, 517758, 635905, 540280, 499306, 499912, 263570,
                     217924, 252885, 206888, 163350, 165000, 164710, 164500, 109800, 109800, 107500 ],
    PtBin        = [ 0., 15., 20, 30., 50, 80., 120., 170., 230., 300.,
                     380, 470., 600., 800., 1000., 1400., 1800., 2200., 2600., 3000. ],
    )
    )

