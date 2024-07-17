import os, sys, enum
import ROOT as r
from copy import deepcopy

class ch(enum.IntEnum):
    NoChan = 0
    ElMu   = 1
    Muon   = 2
    Elec   = 3
    ElMuFromTaus = 4
    MuonFromTaus = 5
    ElecFromTaus = 6
    ElMuMixedFromTaus = 7
    MuonMixedFromTaus = 8
    ElecMixedFromTaus = 9


class tags(enum.IntEnum):
    NoTag      = 0
    mc         = 1
    singlemuon = 2
    doublemuon = 3
    muon       = 4
    egamma     = 5
    muoneg     = 6
    jetmet     = 7

emass = 0.0005109989461

# =========================================================================================================================
# ============================================================================================ TRIGGER
# =========================================================================================================================
def _fires(ev, path):
    return getattr(ev, path)

triggerGroups = dict(
    Trigger_1e = {
        2022 : lambda ev : _fires(ev, 'HLT_Ele32_WPTight_Gsf'),
    },
    Trigger_1m = {
        2022 : lambda ev :     _fires(ev, 'HLT_IsoMu24'),
    },
    Trigger_2e = {
        2022 : lambda ev : (   _fires(ev, 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL')
                            or _fires(ev, 'HLT_DoubleEle25_CaloIdL_MW')),
    },
    Trigger_2m = {
        2022 : lambda ev :     _fires(ev, 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8'),
    },
    Trigger_em = {
        2022 : lambda ev : (   _fires(ev, 'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL')
                            or _fires(ev, 'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ')),
    },
)

from CMGTools.TTHAnalysis.tools.evtTagger import EvtTagger
Trigger_1e = lambda : EvtTagger('Trigger_1e',  [ lambda ev : triggerGroups['Trigger_1e'][ev.year](ev) ])
Trigger_1m = lambda : EvtTagger('Trigger_1m',  [ lambda ev : triggerGroups['Trigger_1m'][ev.year](ev) ])
Trigger_2e = lambda : EvtTagger('Trigger_2e',  [ lambda ev : triggerGroups['Trigger_2e'][ev.year](ev) ])
Trigger_2m = lambda : EvtTagger('Trigger_2m',  [ lambda ev : triggerGroups['Trigger_2m'][ev.year](ev) ])
Trigger_em = lambda : EvtTagger('Trigger_em',  [ lambda ev : triggerGroups['Trigger_em'][ev.year](ev) ])

remove_overlap_booleans = [ lambda ev : (
                            (  (ev.channel == ch.ElMu and (ev.Trigger_em or ev.Trigger_1m or ev.Trigger_1e))
                            or (ev.channel == ch.Muon and (ev.Trigger_1m or ev.Trigger_2m))
                            or (ev.channel == ch.Elec and (ev.Trigger_1e or ev.Trigger_2e)) )
                            if ev.datatag == tags.mc else

                            # First emu dataset
                            (  ev.channel == ch.ElMu and ev.Trigger_em)
                            if ev.datatag == tags.muoneg else

                            # Double muon
                            (   ev.channel == ch.Muon and ev.Trigger_2m)
                            if ev.datatag == tags.doublemuon else

                            # Muon
                            (  (ev.channel == ch.ElMu and (not ev.Trigger_em) and ev.Trigger_1m)
                            or (ev.channel == ch.Muon and (ev.Trigger_2m or ev.Trigger_1m)))
                            if ev.datatag == tags.muon else

                            # E gamma
                            (  (ev.channel == ch.ElMu and (not ev.Trigger_em) and (not ev.Trigger_1m) and ev.Trigger_1e)
                            or (ev.channel == ch.Elec and (ev.Trigger_2e or ev.Trigger_1e)))
                            if ev.datatag == tags.egamma else
                            
                            # Single muon
                            (  (ev.channel == ch.ElMu and (not ev.Trigger_em) and ev.Trigger_1m)
                            or (ev.channel == ch.Muon and (not ev.Trigger_2m) and ev.Trigger_1m))
                            if ev.datatag == tags.singlemuon else

                            (False)
                        )]

remove_overlap = lambda : EvtTagger('pass_trigger', remove_overlap_booleans)

triggerSeq = [Trigger_1e, Trigger_1m, Trigger_2e, Trigger_2m, Trigger_em, remove_overlap]

# =========================================================================================================================
# ============================================================================================ ANALYSIS MODULES & SETTINGS
# =========================================================================================================================

# %%%%%%%%%%%%%%%%%%%%%%%%%% tW full Run 3 inclusive & differential
#### LEPTON TREATMENTS ###
IDDict = {}
IDDict["muons"] = {
    "pt"       : 20,
    "eta"      : 2.4,
    "isorelpf" : 0.15,
}

IDDict["elecs"] = {
    "pt"   : 20,
    "eta0" : 1.4442,
    "eta1" : 1.566,
    "eta2" : 2.4,
    "etaLeak" : 1.56, ###### WARNING, HABRIA QUE VETAR ESTO EN GENERACION?
    "dxy_b" : 0.05,
    "dz_b"  : 0.10,
    "dxy_e" : 0.10,
    "dz_e"  : 0.20,
    "etasc_be" : 1.479,
}

IDDict["jets"] = {
    "pt"      : 30,
    "pt2"     : 20,
    "ptmin"   : 15,
    "ptfwdnoise" : 60,
    "eta"     : 2.4,
    "etafwd"  : 5.0,
    "etafwdnoise0" : 2.7,
    "etafwdnoise1" : 3.0,
    "jetid_2022" : 1,  # > X placeholder, it is hardcoded in the recleaner
}


dresslepID         = lambda l : ( (abs(l.eta) < IDDict["muons"]["eta"] and l.pt > IDDict["muons"]["pt"]) if (abs(l.pdgId) == 13) else
                                  (abs(l.eta) < IDDict["elecs"]["eta2"] and (abs(l.eta) < IDDict["elecs"]["eta0"] or abs(l.eta) > IDDict["elecs"]["eta1"])
                                   and l.pt > IDDict["elecs"]["pt"]) if (abs(l.pdgId) == 11) else False )
dressjetID         = lambda j : ( abs(j.eta) < 2.4 and j.pt > IDDict["jets"]["pt"] )
dressloosejetID    = lambda j : ( abs(j.eta) < 2.4 and j.pt > IDDict["jets"]["pt2"] and j.pt < IDDict["jets"]["pt"] )
dressfwdjetID      = lambda j : ( abs(j.eta) >= 2.4 and abs(j.eta) < 5.0 and
                                ( ((abs(j.eta) < 2.7 or abs(j.eta) >= 3.0) and j.pt > IDDict["jets"]["pt"]) or
                                ((abs(j.eta) >= 2.7 or abs(j.eta) < 3.0) and j.pt > IDDict["jets"]["ptfwdnoise"]) ) )
dressfwdloosejetID = lambda j : ( abs(j.eta) >= 2.4 and abs(j.eta) < 5.0 and (abs(j.eta) < 2.7 or abs(j.eta) >= 3.0) and j.pt > IDDict["jets"]["pt2"] )


#### JET TREATMENTS ###
#jecGroups = {'HF'                   : ['PileUpPtHF', 'RelativeJERHF', 'RelativePtHF'],
#             'BBEC1_{year}'         : ['RelativeJEREC1', 'RelativePtEC1', 'RelativeStatEC'],
#             'FlavorQCD'            : ['FlavorQCD'],
#             'RelativeSample_{year}': ['RelativeSample'],
#             'EC2'                  : ['PileUpPtEC2'],
#             'HF_{year}'            : ['RelativeStatHF'],
#             'RelativeBal'          : ['RelativeBal'],
#             'Absolute_{year}'      : ['AbsoluteStat', 'RelativeStatFSR', 'TimePtEta'],
#             'BBEC1'                : ['PileUpPtBB', 'PileUpPtEC1', 'RelativePtBB'],
#             'EC2_{year}'           : ['RelativeJEREC2', 'RelativePtEC2'],
#             'Absolute'             : ['AbsoluteMPFBias', 'AbsoluteScale', 'Fragmentation', 'PileUpDataMC',
#                                       'PileUpPtRef', 'RelativeFSR', 'SinglePionECAL', 'SinglePionHCAL'],
#}

#jecGroupsFull = []
#for v in jecGroups:
#    jecGroupsFull = jecGroupsFull + jecGroups[v]

from CMGTools.TTHAnalysis.tools.nanoAOD.calculateJECS import JetEnergyCorrector
from CMGTools.TTHAnalysis.tools.nanoAOD.jetMetGrouper_TopRun3 import jetMetCorrelate2022, groups

# Important, the era selects the jet veto map, the jec selects the JEC version
addJECs_2022_mc = lambda : JetEnergyCorrector(
    year = 2022, era = "CD", jec = "Summer22_22Sep2023", jer = "Summer22EEPrompt22", jecveto = "Summer22_23Sep2023", isMC = True,
    algo = "AK4PFPuppi", metbranchname = "PuppiMET", rhoBranchName = "Rho_fixedGridRhoFastjetAll",
    hjetvetomap = "jetvetomap",
    unc = "Total", saveMETUncs = ["T1", "T1Smear"], 
    splitJers = False, applyVetoMaps = True
)

addJECs_2022_data = lambda : JetEnergyCorrector(
    year = 2022, era = "CD", jec = "Summer22_22Sep2023", jer = "Summer22EEPrompt22", jecveto = "Summer22_23Sep2023", isMC = False,
    algo = "AK4PFPuppi", metbranchname = "PuppiMET", rhoBranchName = "Rho_fixedGridRhoFastjetAll",
    hjetvetomap = "jetvetomap",
    unc = "Total", saveMETUncs = ["T1", "T1Smear"], 
    splitJers = False, applyVetoMaps = True
)

addJECs_2022EE_mc = lambda : JetEnergyCorrector(
    year = 2022, era = "EFG", jec = "Summer22EE_22Sep2023", jer = "Summer22EEPrompt22", jecveto = "Summer22EE_23Sep2023", isMC = True,
    algo = "AK4PFPuppi", metbranchname = "PuppiMET", rhoBranchName = "Rho_fixedGridRhoFastjetAll",
    hjetvetomap = "jetvetomap",
    unc = "Total", saveMETUncs = ["T1", "T1Smear"], 
    splitJers = False, applyVetoMaps = True
)

addJECs_2022EE_data = lambda : JetEnergyCorrector(
    year = 2022, era = "EFG", jec = "Summer22EE_22Sep2023", jer = "Summer22EEPrompt22", jecveto = "Summer22EE_23Sep2023", isMC = False,
    algo = "AK4PFPuppi", metbranchname = "PuppiMET", rhoBranchName = "Rho_fixedGridRhoFastjetAll",
    hjetvetomap = "jetvetomap",
    unc = "Total", saveMETUncs = ["T1", "T1Smear"],
    splitJers = False, applyVetoMaps = True
)


# --- 2022: Add JECs (+ correlate in case of MC) --- #
jecs_mc_2022   = [addJECs_2022_mc, jetMetCorrelate2022]
jecs_data_2022 = [addJECs_2022_data]

# --- 2022EE: Add JECs (+ correlate in case of MC) --- #
jecs_mc_2022PostEE   = [addJECs_2022EE_mc, jetMetCorrelate2022]
jecs_data_2022PostEE = [addJECs_2022EE_data]


# Cleaning
from CMGTools.TTHAnalysis.tools.nanoAOD.pythonCleaningTopRun3 import pythonCleaningTopRun2UL

cleaning_mc_2022 = lambda : pythonCleaningTopRun2UL(label  = "Recl",
                                             jetPts = [IDDict["jets"]["pt"], IDDict["jets"]["pt2"]],
                                             jetPtNoisyFwd = IDDict["jets"]["ptfwdnoise"],
                                             jecvars   = ['jesTotal', 'jer'] + ['jes' + v for v in groups],
                                             lepenvars = ["mu","elscale","elsigma"],
                                             isMC      = True,
                                             year_     = "2022", # Selects the correct json for the btag WPs
                                             algo      = "robustParticleTransformer",
                                             btagWP_ = "M",
)
cleaning_mc_2022PostEE = lambda : pythonCleaningTopRun2UL(label  = "Recl",
                                             jetPts = [IDDict["jets"]["pt"], IDDict["jets"]["pt2"]],
                                             jetPtNoisyFwd = IDDict["jets"]["ptfwdnoise"],
                                             jecvars   = ['jesTotal', 'jer'] + ['jes' + v for v in groups],
                                             lepenvars = ["mu","elscale","elsigma"],
                                             isMC      = True,
                                             year_     = "2022PostEE", # Selects the correct json for the btag WPs
                                             algo      = "robustParticleTransformer",
                                             btagWP_ = "M",
)

cleaning_data_2022 = lambda : pythonCleaningTopRun2UL(label = "Recl",
                                               jetPts = [IDDict["jets"]["pt"], IDDict["jets"]["pt2"]],
                                               jetPtNoisyFwd = IDDict["jets"]["ptfwdnoise"],
                                               jecvars   = [], lepenvars = [], isMC = False,
                                               year_     = "2022", # Selects the correct json for the btag WPs
                                               algo      = "robustParticleTransformer",
                                               btagWP_ = "M",
)
cleaning_data_2022PostEE = lambda : pythonCleaningTopRun2UL(label = "Recl",
                                               jetPts = [IDDict["jets"]["pt"], IDDict["jets"]["pt2"]],
                                               jetPtNoisyFwd = IDDict["jets"]["ptfwdnoise"],
                                               jecvars   = [], lepenvars = [], isMC = False,
                                               year_     = "2022PostEE",  # Selects the correct json for the btag WPs
                                               algo      = "robustParticleTransformer",
                                               btagWP_ = "M",
)

#### Add Rochester corrections
from CMGTools.TTHAnalysis.tools.addExtraLepVarsForLepUncs import addExtraLepVarsForLepUncs
# from CMGTools.TTHAnalysis.tools.nanoAOD.jetMetGrouper_TopRun2UL import jetMetCorrelate_TopRun2
addLepUncsVars_mc   = lambda : addExtraLepVarsForLepUncs(elSigmaOrScale = True)
addLepUncsVars_data = lambda : addExtraLepVarsForLepUncs(isMC = False, elSigmaOrScale = True)

from CMGTools.TTHAnalysis.tools.nanoAOD.selectParticleAndPartonInfo import selectParticleAndPartonInfo
theDressAndPartInfo = lambda : selectParticleAndPartonInfo(dresslepSel_         = dresslepID,
                                                           dressjetSel_         = dressjetID,
                                                           dressloosejetSel_    = dressloosejetID,
                                                           dressfwdjetSel_      = dressfwdjetID,
                                                           dressfwdloosejetSel_ = dressfwdloosejetID)

lepsuncsAndParticle_mc   = [addLepUncsVars_mc, theDressAndPartInfo]
lepsuncsAndParticle_data = [addLepUncsVars_data]


#### EVENT VARIABLES ### FTREE 3
from CMGTools.TTHAnalysis.tools.eventVars_TopRun3 import EventVars_TopRun2UL
eventVars_mc_2022   = lambda : EventVars_TopRun2UL('', 'Recl',
                                              jecvars = ['jesTotal', 'jer'] + ['jes' + v for v in groups] + ["unclustEn"],
                                              lepvars = ["mu","elscale","elsigma"], metBranchName='PuppiMET')
eventVars_data = lambda : EventVars_TopRun2UL('', 'Recl', isMC = False,
                                              jecvars = [],
                                              lepvars = [""], metBranchName='PuppiMET')

from CMGTools.TTHAnalysis.tools.particleAndPartonVars_TopRun3 import particleAndPartonVars_TopRun2UL
theDressAndPartVars = lambda : particleAndPartonVars_TopRun2UL()

varstrigger_mc_2022         = [eventVars_mc_2022, theDressAndPartVars] + triggerSeq
varstrigger_mc_2022PostEE   = [eventVars_mc_2022, theDressAndPartVars] + triggerSeq
varstrigger_data            = [eventVars_data] + triggerSeq


### FTREE 4
from CMGTools.TTHAnalysis.tools.nanoAOD.TopPtWeight import TopPtWeight
addTopPtWeight = lambda : TopPtWeight()

# Reweight 13 to 13p6 weight
from CMGTools.TTHAnalysis.tools.nanoAOD.TopPtWeight13to13p6 import TopPtWeight13to13p6
addTopPtWeight13to13p6 = lambda : TopPtWeight13to13p6()



## PU weights
#from CMGTools.TTHAnalysis.tools.nanoAOD.applyPuWeights import puWeighter
#puweight_file = os.path.join( os.environ["CMSSW_BASE"], "src/CMGTools/TTHAnalysis/data/pileup/puWeights_UL_2022.root")
#puweighter = lambda : puWeighter(filename = puweight_file)

###### PU weights producer (I think this is another way of computing the weights starting from only the data PU profile)
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeightProducer
pufile_data2022PostEE = "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/pileup/MyDataPileupHistogram_2022PostEE.root" % os.environ[
    'CMSSW_BASE']
pufile_data2022PreEE = "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/pileup/MyDataPileupHistogram_2022PreEE.root" % os.environ[
    'CMSSW_BASE']
puweighter2022 = lambda : puWeightProducer(
    "auto", pufile_data2022PreEE, "pu_mc", "pileup", verbose=False)
puweighter2022PostEE = lambda : puWeightProducer(
    "auto", pufile_data2022PostEE, "pu_mc", "pileup", verbose=False)


from CMGTools.TTHAnalysis.tools.nanoAOD.btag_weighterRun3 import btag_weighterRun3
## b-tagging
btagpath = os.environ['CMSSW_BASE'] + "/src/CMGTools/TTHAnalysis/data/TopRun3/btagging"
btagWeights_2022       = lambda : btag_weighterRun3(json = btagpath + "/2022_Summer22/" + "btagging_v0.json",
                                            eff = btagpath + "/btagEffs_2024_01_16_btagEff.root",
                                            json_ptrel = btagpath + "/2022_Summer22/" + "btagging_v0.json",
                                            algo = 'robustParticleTransformer',
                                            jecvars   = ['jesTotal', 'jer'] + ['jes' + v for v in groups],
                                            lepenvars = ["mu","elscale","elsigma"],
                                            splitCorrelations = True, ## TEMPORAL
                                            useCombnuisances = True, ## TEMPORAL
                                            year = "2022",
                                            SFmeasReg = "comb")
btagWeights_2022PostEE = lambda : btag_weighterRun3(json = btagpath + "/2022_Summer22EE/" + "btagging_v0.json",
                                            eff = btagpath + "/btagEffs_2024_01_16_btagEff.root",
                                            json_ptrel = btagpath + "/2022_Summer22EE/" + "btagging_v0.json",
                                            algo = 'robustParticleTransformer',
                                            jecvars   = ['jesTotal', 'jer'] + ['jes' + v for v in groups],
                                            lepenvars = ["mu","elscale","elsigma"],
                                            splitCorrelations = True, ## TEMPORAL
                                            useCombnuisances = True, ## TEMPORAL
                                            year = "2022PostEE",
                                            SFmeasReg = "comb")

# Lepton & trigger SF
from CMGTools.TTHAnalysis.tools.nanoAOD.lepScaleFactors_TopRun3 import lepScaleFactors_TopRun3
leptrigSFs_2022    = lambda : lepScaleFactors_TopRun3(year_ = "2022",    lepenvars = ["mu","elscale","elsigma"])
leptrigSFs_2022PostEE    = lambda : lepScaleFactors_TopRun3(year_ = "2022PostEE",    lepenvars = ["mu","elscale","elsigma"])

sfSeq_2022            = [leptrigSFs_2022, puweighter2022, addTopPtWeight, btagWeights_2022]
sfSeq_2022PostEE      = [leptrigSFs_2022PostEE, puweighter2022PostEE, addTopPtWeight, btagWeights_2022PostEE]


### BDT
#import importlib
#
#mvas_mc = [lambda : getattr(importlib.import_module("CMGTools.TTHAnalysis.tools.nanoAOD.MVA_tWRun3"), "MVA_tWRun3_")()]
#
#tmpstr = """mvas_mc.append(lambda : getattr(importlib.import_module("CMGTools.TTHAnalysis.tools.nanoAOD.MVA_tWRun3"), "MVA_tWRun3_{v}{sv}")() )"""
#
#for v in (['jesTotal', 'jer'] + ['jes' + v for v in groups] + ["mu"] + ["unclustEn"]):
#    for sv in ["Up", "Down"]:
#        eval(tmpstr.format(v = v, sv = sv))
#        
#mvas_data = lambda : getattr(importlib.import_module("CMGTools.TTHAnalysis.tools.nanoAOD.MVA_tWRun3"), "MVA_tWRun3_")()


###### b-tagging efficiencies (Summary of JSON content: python3 -m correctionlib.cli summary btagging.json.gz)
from CMGTools.TTHAnalysis.tools.btageffVars_tWRun3 import btageffVars_tWRun3
btagEffFtree_2022 = lambda : btageffVars_tWRun3(wp_   = "M",
                                                algo_ = ['deepJet',
                                                         "particleNet",
                                                         "robustParticleTransformer"],
                                                json_  = btagpath + "/2022_Summer22/" + "btagging_v0.json",
                                                json_ptrel_ = btagpath + "/2022_Summer22/" + "btagging_v0.json",
                                                year_ = "2022",
                                                SFmeasReg = "comb")
                            
btagEffFtree_2022PostEE = lambda : btageffVars_tWRun3(wp_   = "M",
                                                algo_ = ['deepJet',
                                                         "particleNet",
                                                         "robustParticleTransformer"],
                                                json_  = btagpath + "/2022_Summer22EE/" + "btagging_v0.json",
                                                json_ptrel_ = btagpath + "/2022_Summer22EE/" + "btagging_v0.json",
                                                year_ = "2022PostEE",
                                                SFmeasReg = "comb")



###### New MVA without TMVA: twRun3_MVA_SkLearn.py
from CMGTools.TTHAnalysis.tools.nanoAOD.twRun3_MVA_SkLearn import tW_MVA
path_1j1b_newMVA = "/nfs/fanae/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/MVA-Training/onnxConverter/rf1j1b_full2022.onnx"
path_1j1b_mm_newMVA = "/nfs/fanae/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/MVA-Training/onnxConverter/rf1j1b_full2022.onnx"
path_1j1b_ee_newMVA = "/nfs/fanae/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/MVA-Training/onnxConverter/rf1j1b_full2022.onnx"
path_2j1b_newMVA = "/nfs/fanae/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/MVA-Training/onnxConverter/rf2j1b_full2022.onnx"
mvaNew_mc   = [lambda : tW_MVA('', path_1j1b_newMVA, path_2j1b_newMVA, path_1j1b_mm_newMVA, path_1j1b_ee_newMVA,
                                              jecvars = ['jesTotal', 'jer'] + ['jes' + v for v in groups] + ["unclustEn"],
                                              lepvars = ['mu', 'elscale', 'elsigma'])]
mvaNew_data = [lambda : tW_MVA('', path_1j1b_newMVA, path_2j1b_newMVA, path_1j1b_mm_newMVA, path_1j1b_ee_newMVA, isMC = False,
                                              jecvars = [],
                                              lepvars = [""])]

###### Minitrees for training
from CMGTools.TTHAnalysis.tools.nanoAOD.createTrainingMiniTree_TopRun3 import createTrainingMiniTree_TopRun3

createMVAMiniTree = lambda : createTrainingMiniTree_TopRun3()



###### Test for Joscha PU reweighting (calorew)
from CMGTools.TTHAnalysis.tools.nanoAOD.centralCaloWeightPU import centralCaloWeightPU
addcentralCaloWeightPU_2022       = lambda : centralCaloWeightPU(year = "2022")
addcentralCaloWeightPU_2022PostEE = lambda : centralCaloWeightPU(year = "2022PostEE")