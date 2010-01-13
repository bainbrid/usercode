import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.SusyCAF.SusyCAF_Event_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_Track_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_Triggers_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_L1Triggers_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_L1GlobalTrigger_cfi import *
#from SUSYBSMAnalysis.SusyCAF.SusyCAF_L1CaloTrigger_cfi import *
#from SUSYBSMAnalysis.SusyCAF.SusyCAF_L1Extra_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_MET_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_Jet_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_Photon_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_Muon_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_Electron_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_BeamSpot_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_Vertex_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_HcalNoiseSummary_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_HcalNoiseRBX_cfi import *
#from SUSYBSMAnalysis.SusyCAF.SusyCAF_HcalRecHit_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_CaloTowers_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_PFTau_cfi import *

from bainbrid.SUSY.SusyCAF_Jpt_cfi import *

susyTree = cms.EDAnalyzer("SusyTree",
    outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_susycaf*_*_*'
    ))

theBigNtuplePat = cms.Sequence( (susycafevent +
                                 susycaftrack +
                                 susycafl1globaltrigger +  # to be dropped when all L1 triggers have names
                                 susycafL1triggers +       # susycafL1triggersP1 + susycafL1triggersM1 + # susycafL1triggersP2 + susycafL1triggersM2 +
   			         susycaftriggers +
                                 susycafmet +
                                 susycafmetnohf +
                                 susycafmetIC5 +
                                 susycafmetAK5 +
                                 susycafmetPF +
                                 susycafmetTC +
                                 susycafic5calojet +
                                 susycafak5calojet +
                                 susycafsc5calojet +
                                 susycafak5sjptjetreco +
                                 susycafak5vjptjetreco +
                                 susycafphoton +
                                 susycafmuon +
                                 susycafpfmuon +
                                 susycafbeamspot +
                                 susycafvertex +
                                 susycafelectron +
                                 susycafpfelectron +
                                 susycafhcalnoisesummary +
                                 susycafhcalnoiserbx +
                                 susycafcalotowers +
                                 susycaftau) *
                                 susyTree
                                )
