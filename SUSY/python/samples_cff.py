from icf.core import PSet

# -----------------------------------------------------------------------------
# LMx

lm_path="dcap://gfe02.grid.hep.ph.ic.ac.uk:22128//pnfs/hep.ph.ic.ac.uk/data/cms/store/user/bainbrid/ICF/automated/"
lm_xsect = [38.93,4.888,0.6027,3.438,1.879,0.4734,0.3104,1.209,0.7300,7.134,0.04778,0.8236,4.414,6.899]
lm_points = []
for point in range(0,14) :
    name = "LM" + str(point)
    file = lm_path + name + ".Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/SusyCAF_Tree"
    lm = PSet(
        Name = name,
        File = [
        file + "_1_2.root",
        file + "_2_1.root",
        file + "_3_1.root",
        file + "_4_1.root",
        file + "_5_1.root",
        file + "_6_1.root",
        file + "_7_1.root",
        file + "_8_1.root",
        file + "_9_1.root",
        ],
        CrossSection = lm_xsect[point],
        Format = ("ICF",2),
        )
    lm_points.append(lm)
    
# -----------------------------------------------------------------------------
# QCD Pythia


path_qcd_pythia_merged="/nfs/data6/trommers/ICFNtuple/"
#path_v_jets="/castor/cern.ch/user/j/jad/ICFNtuples7TeV_take2/"
#path_qcd_pythia_merged="/castor/cern.ch/user/b/bainbrid/7TeV/V00-08-04-XX/"

qcd_pythia_merged=PSet(
    Name="QCD_Pythia_Merged",
    Format=("ICF",2),
    File=path_qcd_pythia_merged+"QCDPythia_7TeV_V00-08-04-XX_Skim.root",
    Weights = PSet(
    CrossSection = [ 8.762e+08, 6.041e+07, 9.238e+05, 2.547e+04, 1.256e+03, 8.798e+01, 2.186e+00, 1.122e-02 ],
    Events       = [ 6246300,   5228992,   3203440,   3132800,   3274202,   2143390,   2143921,   1184123   ],
    PtBin        = [ 15.,       30.,       80.,       170.,      300.,      470.,      800.,      1400.     ],
    ),
    )

path_qcd_pythia="/vols/cms02/bm409/"

qcd6 = []

qcd_pythia_15=PSet(
    Name="QCDPythia6_Pt15",
    File=path_qcd_pythia+"QCD_Pythia6_15GeV.root",
    CrossSection=8.762e+08,
    Format=("ICF",2),
    )
qcd6.append(qcd_pythia_15)

qcd_pythia_30=PSet(
    Name="QCDPythia6_Pt30",
    File=path_qcd_pythia+"QCD_Pythia6_30GeV.root",
    CrossSection=6.041e+07,
    Format=("ICF",2),
    )
qcd6.append(qcd_pythia_30)

qcd_pythia_80=PSet(
    Name="QCDPythia6_Pt80",
    File=path_qcd_pythia+"QCD_Pythia6_80GeV.root", 
    CrossSection=9.238e+05,
    Format=("ICF",2),
    )
qcd6.append(qcd_pythia_80)

qcd_pythia_170=PSet(
    Name="QCDPythia6_Pt170",
    File=path_qcd_pythia+"QCD_Pythia6_170GeV.root",
    CrossSection=2.547e+04,
    Format=("ICF",2),
    )
qcd6.append(qcd_pythia_170)

qcd_pythia_300=PSet(
    Name="QCDPythia6_Pt300",
    File=path_qcd_pythia+"QCD_Pythia6_300GeV.root",
    CrossSection=1.256e+03,
    Format=("ICF",2),
    )
qcd6.append(qcd_pythia_300)

qcd_pythia_470=PSet(
    Name="QCDPythia6_Pt470",
    File=path_qcd_pythia+"QCD_Pythia6_470GeV.root",
    CrossSection=8.798e+01,
    Format=("ICF",2),
    )
qcd6.append(qcd_pythia_470)

qcd_pythia_800=PSet(
    Name="QCDPythia6_Pt800",
    File=path_qcd_pythia+"QCD_Pythia6_800GeV.root",
    CrossSection=2.186e+00,
    Format=("ICF",2),
    )
qcd6.append(qcd_pythia_800)

qcd_pythia_1400=PSet(
    Name="QCDPythia6_Pt1400",
    File=path_qcd_pythia+"QCD_Pythia6_1400GeV.root",
    CrossSection=1.122e-02,
    Format=("ICF",2),
    )
qcd6.append(qcd_pythia_1400)

qcd_pythia=[
    qcd_pythia_15,
    qcd_pythia_30,
    qcd_pythia_80,
    qcd_pythia_170,
    qcd_pythia_300,
    qcd_pythia_470,
    qcd_pythia_800,
    qcd_pythia_1400
    ]

# -----------------------------------------------------------------------------
# V+Jets

path_v_jets="/nfs/data6/trommers/ICFNtuple/"
#path_v_jets="/castor/cern.ch/user/j/jad/ICFNtuples7TeV_take2/"

w_jets=PSet(
    Name="WJets",
    File=path_v_jets+"WJets-madgraph_ICFv2.root",
    CrossSection=24170.,
    Format=("ICF",2),
    )

z_jets=PSet(
    Name="ZJets",
    File=path_v_jets+"ZJets-madgraph_ICFv2.root",
    CrossSection=2400.,
    Format=("ICF",2),
    )

ttbar_jets=PSet(
    Name="TTbarJets",
    File=path_v_jets+"TTbarJets-madgraph_ICFv2.root",
    CrossSection=95.,
    Format=("ICF",2),
    )

path_z_inv="/nfs/data6/trommers/ICFNtuple/"

z_inv=PSet(
    Name="Zinv",
    File=path_z_inv+"Zinv_7TeV_V00-08-04-XX.root",
    CrossSection=4500.,
    Format=("ICF",2),
    )
