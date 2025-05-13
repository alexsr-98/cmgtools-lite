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
    singleelec = 3
    doublemuon = 4
    doubleeg   = 5
    muoneg     = 6

emass = 0.0005109989461

# =========================================================================================================================
# ============================================================================================ TRIGGER
# =========================================================================================================================
#def _fires(ev, path):
#    if "/hasfiredtriggers_cc.so" not in r.gSystem.GetLibraries():
#        r.gROOT.ProcessLine(".L %s/src/CMGTools/Production/src/hasfiredtriggers.cc+" % os.environ['CMSSW_BASE'])
#    if not hasattr(ev, path): return False

#    if ev.run == 1:  # is MC
#        return getattr(ev, path)
#    else :
#        return getattr(r, 'fires_%s_%d'%(path, ev.year))( ev.run, getattr(ev, path) )


def _fires(ev, path):
    return getattr(ev, path)


triggerGroups = dict(
    Trigger_1e = {
        2016 : lambda ev : _fires(ev, 'HLT_Ele27_WPTight_Gsf'),
        2017 : lambda ev : _fires(ev, 'HLT_Ele35_WPTight_Gsf'),
        2018 : lambda ev : _fires(ev, 'HLT_Ele32_WPTight_Gsf'),
    },
    Trigger_1m = {
        2016 : lambda ev : (   _fires(ev, 'HLT_IsoMu24')
                            or _fires(ev, 'HLT_IsoTkMu24')),
        2017 : lambda ev :     _fires(ev, 'HLT_IsoMu27'),
        2018 : lambda ev :     _fires(ev, 'HLT_IsoMu24'),
    },
    Trigger_2e = {
        2016 : lambda ev : ((   _fires(ev, 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ')
                             or _fires(ev, 'HLT_DoubleEle33_CaloIdL_MW')
                             or _fires(ev, 'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL'))
                             if (ev.datatag == tags.mc or ev.run <= 280385) else
                            (   _fires(ev, 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ')
                             or _fires(ev, 'HLT_DoubleEle33_CaloIdL_MW'))),
        2017 : lambda ev : (   _fires(ev, 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL')
                            or _fires(ev, 'HLT_DoubleEle33_CaloIdL_MW')),
        2018 : lambda ev : (   _fires(ev, 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL')
                            or _fires(ev, 'HLT_DoubleEle25_CaloIdL_MW')),
    },
    Trigger_2m = {
        2016 : lambda ev : ((   _fires(ev, 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL')
                             or _fires(ev, 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL'))
                            if (ev.datatag == tags.mc or ev.run <= 280385) else
                            (  _fires(ev, 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ')
                            or _fires(ev, 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ'))),
        2017 : lambda ev : ((   _fires(ev, 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ')
                             or _fires(ev, 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8')
                             or _fires(ev, 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8'))
                            if (ev.datatag == tags.mc) else
                            (  _fires(ev, 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ'))
                            if (ev.run < 299368) else
                            (  _fires(ev, 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8')
                            or _fires(ev, 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8'))),
        2018 : lambda ev :     _fires(ev, 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8'),
    },
    Trigger_em = {
        2016 : lambda ev : ((  _fires(ev, 'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL')
                            or _fires(ev, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"))
                            if (ev.datatag == tags.mc or ev.run <= 280385) else
                            (  _fires(ev, 'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ')
                            or _fires(ev, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"))),
        2017 : lambda ev : ((   _fires(ev, 'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL')
                             or _fires(ev, 'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ'))
                             if (ev.datatag == tags.mc or ev.run > 299329) else
                             (  _fires(ev, 'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ'))),
        2018 : lambda ev : (   _fires(ev, 'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL')
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

                            (  (ev.channel == ch.ElMu and (not ev.Trigger_em) and ev.Trigger_1m)
                            or (ev.channel == ch.Muon and (not ev.Trigger_2m) and ev.Trigger_1m))
                            if ev.datatag == tags.singlemuon else

                            (   ev.channel == ch.ElMu and (not ev.Trigger_em) and (not ev.Trigger_1m) and ev.Trigger_1e)
                            or (ev.channel == ch.Elec and (not ev.Trigger_2e) and ev.Trigger_1e)
                            if (ev.datatag == tags.singleelec and not ev.year == 2018) else

                            (   ev.channel == ch.ElMu and (not ev.Trigger_em) and (not ev.Trigger_1m) and ev.Trigger_1e)
                            or (ev.channel == ch.Elec and (ev.Trigger_2e or ev.Trigger_1e))
                            if (ev.datatag == tags.singleelec and ev.year == 2018) else
                            
                            (   ev.channel == ch.Muon and ev.Trigger_2m)
                            if ev.datatag == tags.doublemuon else

                            (  ev.channel == ch.Elec and ev.Trigger_2e)
                            if ev.datatag == tags.doubleeg else

                            (  ev.channel == ch.ElMu and ev.Trigger_em)
                            if ev.datatag == tags.muoneg else

                            (False)
                        )]

remove_overlap = lambda : EvtTagger('pass_trigger', remove_overlap_booleans)

triggerSeq = [Trigger_1e, Trigger_1m, Trigger_2e, Trigger_2m, Trigger_em, remove_overlap]

# =========================================================================================================================
# ============================================================================================ ANALYSIS MODULES & SETTINGS
# =========================================================================================================================

# %%%%%%%%%%%%%%%%%%%%%%%%%% tW full Run 2 inclusive & differential
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
    "jetid_2016" : 1,  # > X
    "jetid_2017" : 1,  # > X
    "jetid_2018" : 1,  # > X
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
jecGroups = {'HF'                   : ['PileUpPtHF', 'RelativeJERHF', 'RelativePtHF'],
             'BBEC1_{year}'         : ['RelativeJEREC1', 'RelativePtEC1', 'RelativeStatEC'],
             'FlavorQCD'            : ['FlavorQCD'],
             'RelativeSample_{year}': ['RelativeSample'],
             'EC2'                  : ['PileUpPtEC2'],
             'HF_{year}'            : ['RelativeStatHF'],
             'RelativeBal'          : ['RelativeBal'],
             'Absolute_{year}'      : ['AbsoluteStat', 'RelativeStatFSR', 'TimePtEta'],
             'BBEC1'                : ['PileUpPtBB', 'PileUpPtEC1', 'RelativePtBB'],
             'EC2_{year}'           : ['RelativeJEREC2', 'RelativePtEC2'],
             'Absolute'             : ['AbsoluteMPFBias', 'AbsoluteScale', 'Fragmentation', 'PileUpDataMC',
                                       'PileUpPtRef', 'RelativeFSR', 'SinglePionECAL', 'SinglePionHCAL'],
}


from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import createJMECorrector

addJECs_2016apv = createJMECorrector(dataYear      = "UL2016_preVFP",
                                     jesUncert     = "Merged",
                                     metBranchName = "MET",
                                     splitJER      = True,
                                     applyHEMfix   = True,)
addJECs_2016    = createJMECorrector(dataYear      = "UL2016",
                                     jesUncert     = "Merged",
                                     metBranchName = "MET",
                                     splitJER      = True,
                                     applyHEMfix   = True,)
addJECs_2017    = createJMECorrector(dataYear      = "UL2017",
                                     jesUncert     = "Merged",
                                     metBranchName = "MET",
                                     splitJER      = True,
                                     applyHEMfix   = True,)
addJECs_2018    = createJMECorrector(dataYear      = "UL2018",
                                     jesUncert     = "Merged",
                                     metBranchName = "MET",
                                     splitJER      = True,
                                     applyHEMfix   = True,)

jecs_2016apv = addJECs_2016apv
jecs_2016    = addJECs_2016
jecs_2017    = addJECs_2017
jecs_2018    = addJECs_2018

# Cleaning
from CMGTools.TTHAnalysis.tools.nanoAOD.pythonCleaningTopRun2UL import pythonCleaningTopRun2UL
cleaning_mc_2016apv = lambda : pythonCleaningTopRun2UL(label  = "Recl",
                                             jetPts = [IDDict["jets"]["pt"], IDDict["jets"]["pt2"]],
                                             jetPtNoisyFwd = IDDict["jets"]["ptfwdnoise"],
                                             jecvars   = ['jesTotal', 'jer'] + ['jes%s'%v.format(year = 2016) for v in jecGroups] + ["jer%i"%i for i in range(6)] + ["jesHEMIssue"],
                                             lepenvars = ["mu", "elsigma"],
                                             isMC      = True,
                                             #debug     = True,
                                             year_     = "2016apv",
#                                             algo      = "DeepCSV",
)
cleaning_mc_2016 = lambda : pythonCleaningTopRun2UL(label  = "Recl",
                                             jetPts = [IDDict["jets"]["pt"], IDDict["jets"]["pt2"]],
                                             jetPtNoisyFwd = IDDict["jets"]["ptfwdnoise"],
                                             jecvars   = ['jesTotal', 'jer'] + ['jes%s'%v.format(year = 2016) for v in jecGroups] + ["jer%i"%i for i in range(6)] + ["jesHEMIssue"],
                                             lepenvars = ["mu", "elsigma"],
                                             isMC      = True,
                                             #debug     = True,
                                             year_     = "2016",
#                                             algo      = "DeepCSV",
)
cleaning_mc_2017 = lambda : pythonCleaningTopRun2UL(label  = "Recl",
                                             jetPts = [IDDict["jets"]["pt"], IDDict["jets"]["pt2"]],
                                             jetPtNoisyFwd = IDDict["jets"]["ptfwdnoise"],
                                             jecvars   = ['jesTotal', 'jer'] + ['jes%s'%v.format(year = 2017) for v in jecGroups] + ["jer%i"%i for i in range(6)] + ["jesHEMIssue"],
                                             lepenvars = ["mu", "elsigma"],
                                             isMC      = True,
                                             #debug     = True,
                                             year_     = "2017",
#                                             algo      = "DeepCSV",
)
cleaning_mc_2018 = lambda : pythonCleaningTopRun2UL(label  = "Recl",
                                             jetPts = [IDDict["jets"]["pt"], IDDict["jets"]["pt2"]],
                                             jetPtNoisyFwd = IDDict["jets"]["ptfwdnoise"],
                                             jecvars   = ['jesTotal', 'jer'] + ['jes%s'%v.format(year = 2018) for v in jecGroups] + ["jer%i"%i for i in range(6)] + ["jesHEMIssue"],
                                             lepenvars = ["mu", "elsigma"],
                                             isMC      = True,
                                             year_     = "2018",
                                             #debug     = True,
#                                             algo      = "DeepCSV",
)

cleaning_data_2016apv = lambda : pythonCleaningTopRun2UL(label = "Recl",
                                               jetPts = [IDDict["jets"]["pt"], IDDict["jets"]["pt2"]],
                                               jetPtNoisyFwd = IDDict["jets"]["ptfwdnoise"],
                                               jecvars   = [], lepenvars = ["elscale"], isMC = False,
                                               year_     = "2016apv",
#                                               algo      = "DeepCSV",
)
cleaning_data_2016 = lambda : pythonCleaningTopRun2UL(label = "Recl",
                                               jetPts = [IDDict["jets"]["pt"], IDDict["jets"]["pt2"]],
                                               jetPtNoisyFwd = IDDict["jets"]["ptfwdnoise"],
                                               jecvars   = [], lepenvars = ["elscale"], isMC = False,
                                               year_     = "2016",
#                                               algo      = "DeepCSV",
)
cleaning_data_2017 = lambda : pythonCleaningTopRun2UL(label = "Recl",
                                               jetPts = [IDDict["jets"]["pt"], IDDict["jets"]["pt2"]],
                                               jetPtNoisyFwd = IDDict["jets"]["ptfwdnoise"],
                                               jecvars   = [], lepenvars = ["elscale"], isMC = False,
                                               year_     = "2017",
#                                               algo      = "DeepCSV",
)
cleaning_data_2018 = lambda : pythonCleaningTopRun2UL(label = "Recl",
                                               jetPts = [IDDict["jets"]["pt"], IDDict["jets"]["pt2"]],
                                               jetPtNoisyFwd = IDDict["jets"]["ptfwdnoise"],
                                               jecvars   = [], lepenvars = ["elscale"], isMC = False,
                                               year_     = "2018",
#                                               algo      = "DeepCSV",
)

#### Add Rochester corrections
from CMGTools.TTHAnalysis.tools.addExtraLepVarsForLepUncsTopRun2UL import addExtraLepVarsForLepUncs
# from CMGTools.TTHAnalysis.tools.nanoAOD.jetMetGrouper_TopRun2UL import jetMetCorrelate_TopRun2
addLepUncsVars_mc   = lambda : addExtraLepVarsForLepUncs()
addLepUncsVars_data = lambda : addExtraLepVarsForLepUncs(isMC = False)

from CMGTools.TTHAnalysis.tools.nanoAOD.selectParticleAndPartonInfo import selectParticleAndPartonInfo
theDressAndPartInfo = lambda : selectParticleAndPartonInfo(dresslepSel_         = dresslepID,
                                                           dressjetSel_         = dressjetID,
                                                           dressloosejetSel_    = dressloosejetID,
                                                           dressfwdjetSel_      = dressfwdjetID,
                                                           dressfwdloosejetSel_ = dressfwdloosejetID)

lepsuncsAndParticle_mc   = [addLepUncsVars_mc, theDressAndPartInfo]
lepsuncsAndParticle_data = [addLepUncsVars_data]


#### EVENT VARIABLES ### FTREE 3
from CMGTools.TTHAnalysis.tools.eventVars_TopRun2UL import EventVars_TopRun2UL
eventVars_mc_2016   = lambda : EventVars_TopRun2UL('', 'Recl',
                                              jecvars = ['jesTotal', 'jer'] + ['jes%s'%v.format(year = 2016) for v in jecGroups] + ["jer%i"%i for i in range(6)] + ["unclustEn", "jesHEMIssue"],
                                              lepvars = ['mu', "elsigma"])
eventVars_mc_2017   = lambda : EventVars_TopRun2UL('', 'Recl',
                                              jecvars = ['jesTotal', 'jer'] + ['jes%s'%v.format(year = 2017) for v in jecGroups] + ["jer%i"%i for i in range(6)] + ["unclustEn", "jesHEMIssue"],
                                              lepvars = ['mu', "elsigma"])
eventVars_mc_2018   = lambda : EventVars_TopRun2UL('', 'Recl',
                                              jecvars = ['jesTotal', 'jer'] + ['jes%s'%v.format(year = 2018) for v in jecGroups] + ["jer%i"%i for i in range(6)] + ["unclustEn", "jesHEMIssue"],
                                              lepvars = ['mu', "elsigma"])
eventVars_data = lambda : EventVars_TopRun2UL('', 'Recl', isMC = False,
                                              jecvars = [],
                                              lepvars = ["elscale"])

from CMGTools.TTHAnalysis.tools.particleAndPartonVars_TopRun2UL import particleAndPartonVars_TopRun2UL
theDressAndPartVars = lambda : particleAndPartonVars_TopRun2UL()

varstrigger_mc_2016apv = [eventVars_mc_2016, theDressAndPartVars] + triggerSeq
varstrigger_mc_2016    = [eventVars_mc_2016, theDressAndPartVars] + triggerSeq
varstrigger_mc_2017    = [eventVars_mc_2017, theDressAndPartVars] + triggerSeq
varstrigger_mc_2018    = [eventVars_mc_2018, theDressAndPartVars] + triggerSeq
varstrigger_data       = [eventVars_data] + triggerSeq


### FTREE 4
from CMGTools.TTHAnalysis.tools.nanoAOD.TopPtWeight import TopPtWeight
addTopPtWeight = lambda : TopPtWeight()

## b-tagging
from CMGTools.TTHAnalysis.tools.nanoAOD.btag_weighterTopRun2UL import btag_weighterUL
btagpath = os.environ['CMSSW_BASE'] + "/src/CMGTools/TTHAnalysis/data/TopRun2UL/btagging/"
btagWeights_2016apv = lambda : btag_weighterUL(json = btagpath  + "btagging_2016ULpreVFP.json",
                                            eff = btagpath + "btagEffs_2024_09_08_btagEff.root",
                                            algo = 'deepJet',
                                            jecvars   = ['jesTotal', 'jer'] + ['jes%s'%v.format(year = 2016) for v in jecGroups] + ["jer%i"%i for i in range(6)] + ["jesHEMIssue"],
                                            lepenvars = ["mu", "elsigma"],
                                            splitCorrelations = True, ## TEMPORAL
                                            useCombnuisances = False, ## TEMPORAL
                                            year = "2016apv",
                                            SFmeasReg = "mujets")

btagWeights_2016 = lambda : btag_weighterUL(json = btagpath  + "btagging_2016ULpostVFP.json",
                                            eff = btagpath + "btagEffs_2024_09_08_btagEff.root",
                                            algo = 'deepJet',
                                            jecvars   = ['jesTotal', 'jer'] + ['jes%s'%v.format(year = 2016) for v in jecGroups] + ["jer%i"%i for i in range(6)] + ["jesHEMIssue"],
                                            lepenvars = ["mu", "elsigma"],
                                            splitCorrelations = True, ## TEMPORAL
                                            useCombnuisances = False, ## TEMPORAL
                                            year = "2016",
                                            SFmeasReg = "mujets")

btagWeights_2017 = lambda : btag_weighterUL(json = btagpath  + "btagging_2017UL.json",
                                            eff = btagpath + "btagEffs_2024_09_08_btagEff.root",
                                            algo = 'deepJet',
                                            jecvars   = ['jesTotal', 'jer'] + ['jes%s'%v.format(year = 2017) for v in jecGroups] + ["jer%i"%i for i in range(6)] + ["jesHEMIssue"],
                                            lepenvars = ["mu", "elsigma"],
                                            splitCorrelations = True, ## TEMPORAL
                                            useCombnuisances = False, ## TEMPORAL
                                            year = "2017",
                                            SFmeasReg = "mujets")


btagWeights_2018 = lambda : btag_weighterUL(json = btagpath  + "btagging_2018UL.json",
                                            eff = btagpath + "btagEffs_2024_09_08_btagEff.root",
                                            algo = 'deepJet',
                                            jecvars   = ['jesTotal', 'jer'] + ['jes%s'%v.format(year = 2018) for v in jecGroups] + ["jer%i"%i for i in range(6)] + ["jesHEMIssue"],
                                            lepenvars = ["mu", "elsigma"],
                                            splitCorrelations = True, ## TEMPORAL
                                            useCombnuisances = False, ## TEMPORAL
                                            year = "2018",
                                            SFmeasReg = "mujets")

# Lepton & trigger SF
from CMGTools.TTHAnalysis.tools.nanoAOD.lepScaleFactors_TopRun2UL import lepScaleFactors_TopRun2UL
leptrigSFs_2016apv = lambda : lepScaleFactors_TopRun2UL(year_ = "2016apv", lepenvars = ["mu", "elsigma"])
leptrigSFs_2016    = lambda : lepScaleFactors_TopRun2UL(year_ = "2016",    lepenvars = ["mu", "elsigma"])
leptrigSFs_2017    = lambda : lepScaleFactors_TopRun2UL(year_ = "2017",    lepenvars = ["mu", "elsigma"])
leptrigSFs_2018    = lambda : lepScaleFactors_TopRun2UL(year_ = "2018",    lepenvars = ["mu", "elsigma"])


### Con jet PU ID
from CMGTools.TTHAnalysis.tools.nanoAOD.jetPUid_weighterUL import jetPUid_weighterUL
jetpuidpath = os.environ['CMSSW_BASE'] + "/src/CMGTools/TTHAnalysis/data/TopRun2UL/jetPUid"
addjetPUidMod_2016apv = lambda : jetPUid_weighterUL(jetpuidpath + "/" + "jmar_2016ULpreVFP.json",
                                                    wp = "L",
                                                    jecvars   = ['jesTotal', 'jer'] + ['jes%s'%v.format(year = 2016) for v in jecGroups] + ["jer%i"%i for i in range(6)] + ["jesHEMIssue"],
                                                    lepenvars = ["mu", "elsigma"],
                                                    #year = 2016, debug = True)
                                                    year = 2016)
addjetPUidMod_2016    = lambda : jetPUid_weighterUL(jetpuidpath + "/" + "jmar_2016ULpostVFP.json",
                                                    wp = "L",
                                                    jecvars   = ['jesTotal', 'jer'] + ['jes%s'%v.format(year = 2016) for v in jecGroups] + ["jer%i"%i for i in range(6)] + ["jesHEMIssue"],
                                                    lepenvars = ["mu", "elsigma"],
                                                    #year = 2016, debug = True)
                                                    year = 2016)
addjetPUidMod_2017    = lambda : jetPUid_weighterUL(jetpuidpath + "/" + "jmar_2017UL.json",
                                                    wp = "L",
                                                    jecvars   = ['jesTotal', 'jer'] + ['jes%s'%v.format(year = 2017) for v in jecGroups] + ["jer%i"%i for i in range(6)] + ["jesHEMIssue"],
                                                    lepenvars = ["mu", "elsigma"],
                                                    #year = 2017, debug = True)
                                                    year = 2017)
addjetPUidMod_2018    = lambda : jetPUid_weighterUL(jetpuidpath + "/" + "jmar_2018UL.json",
                                                    wp = "L",
                                                    jecvars   = ['jesTotal', 'jer'] + ['jes%s'%v.format(year = 2018) for v in jecGroups] + ["jer%i"%i for i in range(6)] + ["jesHEMIssue"],
                                                    lepenvars = ["mu", "elsigma"],
                                                    #year = 2018, debug = True)
                                                    year = 2018)


sfSeq_2016apv   = [leptrigSFs_2016apv, btagWeights_2016apv, addjetPUidMod_2016apv]#addTopPtWeight
sfSeq_2016      = [leptrigSFs_2016,    btagWeights_2016,    addjetPUidMod_2016]   #addTopPtWeight
sfSeq_2017      = [leptrigSFs_2017,    btagWeights_2017,    addjetPUidMod_2017]   #addTopPtWeight
sfSeq_2018      = [leptrigSFs_2018,    btagWeights_2018,    addjetPUidMod_2018]   #addTopPtWeight



###### b-tagging efficiencies
from CMGTools.TTHAnalysis.tools.btageffVars_TopRun2UL import btageffVars_TopRun2UL
btagEffFtree_2016apv = lambda : btageffVars_TopRun2UL(wp_   = "M",
                                                algo_ = ['deepJet'],
                                                json_  = btagpath + "btagging_2016ULpreVFP.json",
                                                year_ = "2016apv",
                                                SFmeasReg = "mujets")
btagEffFtree_2016 = lambda : btageffVars_TopRun2UL(wp_   = "M",
                                                algo_ = ['deepJet'],
                                                json_  = btagpath + "btagging_2016ULpostVFP.json",
                                                year_ = "2016",
                                                SFmeasReg = "mujets")
btagEffFtree_2017 = lambda : btageffVars_TopRun2UL(wp_   = "M",
                                                algo_ = ['deepJet'],
                                                json_  = btagpath + "btagging_2017UL.json",
                                                year_ = "2017",
                                                SFmeasReg = "mujets")
btagEffFtree_2018 = lambda : btageffVars_TopRun2UL(wp_   = "M",
                                                algo_ = ['deepJet'],
                                                json_  = btagpath + "btagging_2018UL.json",
                                                year_ = "2018",
                                                SFmeasReg = "mujets")


### Friend tree 5
from CMGTools.TTHAnalysis.tools.particleAndPartonVars_tagBjets_TopRun2UL import particleAndPartonVars_tagBjets_TopRun2UL 
tagIsHF = lambda : particleAndPartonVars_tagBjets_TopRun2UL()

