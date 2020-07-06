from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection as Collection
from CMGTools.TTHAnalysis.tools.nanoAOD.friendVariableProducerTools import writeOutput
from CMGTools.TTHAnalysis.tools.mvaTool import *


def obtainIndex(event, selbranch, index):
    theval = getattr(event, selbranch)
    return int(theval[index])



class MVA_tWRun2(Module):
    def __init__(self):
        self.MVAs = {}
        self.vars_1j1t = [
            MVAVar("TnLooseCentral",        func = lambda ev : ev.nJetSel20_Recl),
            MVAVar("TnBTotal",              func = lambda ev : ev.nBJetSelMedium30_Recl),
            MVAVar("TDilepMETJetPt",        func = lambda ev : ev.Lep1Lep2Jet1MET_Pt),
            MVAVar("TTHTtot",               func = lambda ev : ev.HTtot),
            MVAVar("TTJet1_pt",             func = lambda ev : ev.Jet1_Pt),
            MVAVar("TTJetLooseCentralpt",   func = lambda ev : ev.JetLoose1_Pt),
            MVAVar("TDilepMETJetPt_THTtot", func = lambda ev : ev.Lep1Lep2Jet1MET_Pt / ev.HTtot),
            MVAVar("TMSys",                 func = lambda ev : ev.Lep1Lep2Jet1MET_M),
            MVAVar("TC_jll_v2",             func = lambda ev : ev.Lep1Lep2Jet1_C),
            MVAVar("THTLepOverHT",          func = lambda ev : (ev.LepGood_pt[0] + ev.LepGood_pt[1]) / ev.HTtot if ev.nLepGood >= 2 else 0),
            MVAVar("TDilepJetPt",           func = lambda ev : ev.Lep1Lep2Jet1_Pt),
        ]
        self.vars_1j1t_jesup = [
            MVAVar("TnLooseCentral",        func = lambda ev : ev.nJetSel20_jesTotalUp_Recl),
            MVAVar("TnBTotal",              func = lambda ev : ev.nBJetSelMedium30_jesTotalUp_Recl),
            MVAVar("TDilepMETJetPt",        func = lambda ev : ev.Lep1Lep2Jet1MET_Pt_jesTotalUp),
            MVAVar("TTHTtot",               func = lambda ev : ev.HTtot_jesTotalUp),
            MVAVar("TTJet1_pt",             func = lambda ev : ev.Jet1_Pt_jesTotalUp),
            MVAVar("TTJetLooseCentralpt",   func = lambda ev : ev.JetLoose1_Pt_jesTotalUp),
            MVAVar("TDilepMETJetPt_THTtot", func = lambda ev : ev.Lep1Lep2Jet1MET_Pt_jesTotalUp / ev.HTtot_jesTotalUp),
            MVAVar("TMSys",                 func = lambda ev : ev.Lep1Lep2Jet1MET_M_jesTotalUp),
            MVAVar("TC_jll_v2",             func = lambda ev : ev.Lep1Lep2Jet1_C_jesTotalUp),
            MVAVar("THTLepOverHT",          func = lambda ev : (ev.LepGood_pt[0] + ev.LepGood_pt[1]) / ev.HTtot_jesTotalUp if ev.nLepGood >= 2 else 0),
            MVAVar("TDilepJetPt",           func = lambda ev : ev.Lep1Lep2Jet1_Pt_jesTotalUp),
        ]
        self.vars_1j1t_jesdn = [
            MVAVar("TnLooseCentral",        func = lambda ev : ev.nJetSel20_jesTotalDown_Recl),
            MVAVar("TnBTotal",              func = lambda ev : ev.nBJetSelMedium30_jesTotalDown_Recl),
            MVAVar("TDilepMETJetPt",        func = lambda ev : ev.Lep1Lep2Jet1MET_Pt_jesTotalDown),
            MVAVar("TTHTtot",               func = lambda ev : ev.HTtot_jesTotalDown),
            MVAVar("TTJet1_pt",             func = lambda ev : ev.Jet1_Pt_jesTotalDown),
            MVAVar("TTJetLooseCentralpt",   func = lambda ev : ev.JetLoose1_Pt_jesTotalDown),
            MVAVar("TDilepMETJetPt_THTtot", func = lambda ev : ev.Lep1Lep2Jet1MET_Pt_jesTotalDown / ev.HTtot_jesTotalDown),
            MVAVar("TMSys",                 func = lambda ev : ev.Lep1Lep2Jet1MET_M_jesTotalDown),
            MVAVar("TC_jll_v2",             func = lambda ev : ev.Lep1Lep2Jet1_C_jesTotalDown),
            MVAVar("THTLepOverHT",          func = lambda ev : (ev.LepGood_pt[0] + ev.LepGood_pt[1]) / ev.HTtot_jesTotalDown if ev.nLepGood >= 2 else 0),
            MVAVar("TDilepJetPt",           func = lambda ev : ev.Lep1Lep2Jet1_Pt_jesTotalDown),
        ]
        self.vars_1j1t_jerup = [
            MVAVar("TnLooseCentral",        func = lambda ev : ev.nJetSel20_jerUp_Recl),
            MVAVar("TnBTotal",              func = lambda ev : ev.nBJetSelMedium30_jerUp_Recl),
            MVAVar("TDilepMETJetPt",        func = lambda ev : ev.Lep1Lep2Jet1MET_Pt_jerUp),
            MVAVar("TTHTtot",               func = lambda ev : ev.HTtot_jerUp),
            MVAVar("TTJet1_pt",             func = lambda ev : ev.Jet1_Pt_jerUp),
            MVAVar("TTJetLooseCentralpt",   func = lambda ev : ev.JetLoose1_Pt_jerUp),
            MVAVar("TDilepMETJetPt_THTtot", func = lambda ev : ev.Lep1Lep2Jet1MET_Pt_jerUp / ev.HTtot_jerUp),
            MVAVar("TMSys",                 func = lambda ev : ev.Lep1Lep2Jet1MET_M_jerUp),
            MVAVar("TC_jll_v2",             func = lambda ev : ev.Lep1Lep2Jet1_C_jerUp),
            MVAVar("THTLepOverHT",          func = lambda ev : (ev.LepGood_pt[0] + ev.LepGood_pt[1]) / ev.HTtot_jerUp if ev.nLepGood >= 2 else 0),
            MVAVar("TDilepJetPt",           func = lambda ev : ev.Lep1Lep2Jet1_Pt_jerUp),
        ]
        self.vars_1j1t_jerdn = [
            MVAVar("TnLooseCentral",        func = lambda ev : ev.nJetSel20_jerDown_Recl),
            MVAVar("TnBTotal",              func = lambda ev : ev.nBJetSelMedium30_jerDown_Recl),
            MVAVar("TDilepMETJetPt",        func = lambda ev : ev.Lep1Lep2Jet1MET_Pt_jerDown),
            MVAVar("TTHTtot",               func = lambda ev : ev.HTtot_jerDown),
            MVAVar("TTJet1_pt",             func = lambda ev : ev.Jet1_Pt_jerDown),
            MVAVar("TTJetLooseCentralpt",   func = lambda ev : ev.JetLoose1_Pt_jerDown),
            MVAVar("TDilepMETJetPt_THTtot", func = lambda ev : ev.Lep1Lep2Jet1MET_Pt_jerDown / ev.HTtot_jerDown),
            MVAVar("TMSys",                 func = lambda ev : ev.Lep1Lep2Jet1MET_M_jerDown),
            MVAVar("TC_jll_v2",             func = lambda ev : ev.Lep1Lep2Jet1_C_jerDown),
            MVAVar("THTLepOverHT",          func = lambda ev : (ev.LepGood_pt[0] + ev.LepGood_pt[1]) / ev.HTtot_jerDown if ev.nLepGood >= 2 else 0),
            MVAVar("TDilepJetPt",           func = lambda ev : ev.Lep1Lep2Jet1_Pt_jerDown),
        ]

        self.vars_2j1t = [
            #MVAVar("jetPtSubLeading_",      func = lambda ev : ev.Jet_pt[ev.iJetSel30_Recl[1]] if ev.nJetSel30_Recl >= 2 else 0),
            MVAVar("jetPtSubLeading_",      func = lambda ev : ev.Jet2_Pt),
            MVAVar("deltaRL1_J1_",          func = lambda ev : ev.Lep1Jet1_DR),
            MVAVar("deltaRL1L2_J1J2_",      func = lambda ev : ev.Lep12Jet12_DR),
            MVAVar("deltaRL1L2_metJ1J2_",   func = lambda ev : ev.Lep12Jet12MET_DR),
        ]
        self.vars_2j1t_jesup = [
            #MVAVar("jetPtSubLeading_",      func = lambda ev : ev.Jet_pt[ev.iJetSel30_Recl[1]] if ev.nJetSel30_Recl >= 2 else 0),
            MVAVar("jetPtSubLeading_",      func = lambda ev : ev.Jet2_Pt_jesTotalUp),
            MVAVar("deltaRL1_J1_",          func = lambda ev : ev.Lep1Jet1_DR_jesTotalUp),
            MVAVar("deltaRL1L2_J1J2_",      func = lambda ev : ev.Lep12Jet12_DR_jesTotalUp),
            MVAVar("deltaRL1L2_metJ1J2_",   func = lambda ev : ev.Lep12Jet12MET_DR_jesTotalUp),
        ]
        self.vars_2j1t_jesdn = [
            #MVAVar("jetPtSubLeading_",      func = lambda ev : ev.Jet_pt[ev.iJetSel30_Recl[1]] if ev.nJetSel30_Recl >= 2 else 0),
            MVAVar("jetPtSubLeading_",      func = lambda ev : ev.Jet2_Pt_jesTotalDown),
            MVAVar("deltaRL1_J1_",          func = lambda ev : ev.Lep1Jet1_DR_jesTotalDown),
            MVAVar("deltaRL1L2_J1J2_",      func = lambda ev : ev.Lep12Jet12_DR_jesTotalDown),
            MVAVar("deltaRL1L2_metJ1J2_",   func = lambda ev : ev.Lep12Jet12MET_DR_jesTotalDown),
        ]
        self.vars_2j1t_jerup = [
            #MVAVar("jetPtSubLeading_",      func = lambda ev : ev.Jet_pt[ev.iJetSel30_Recl[1]] if ev.nJetSel30_Recl >= 2 else 0),
            MVAVar("jetPtSubLeading_",      func = lambda ev : ev.Jet2_Pt_jerUp),
            MVAVar("deltaRL1_J1_",          func = lambda ev : ev.Lep1Jet1_DR_jerUp),
            MVAVar("deltaRL1L2_J1J2_",      func = lambda ev : ev.Lep12Jet12_DR_jerUp),
            MVAVar("deltaRL1L2_metJ1J2_",   func = lambda ev : ev.Lep12Jet12MET_DR_jerUp),
        ]
        self.vars_2j1t_jerdn = [
            #MVAVar("jetPtSubLeading_",      func = lambda ev : ev.Jet_pt[ev.iJetSel30_Recl[1]] if ev.nJetSel30_Recl >= 2 else 0),
            MVAVar("jetPtSubLeading_",      func = lambda ev : ev.Jet2_Pt_jerDown),
            MVAVar("deltaRL1_J1_",          func = lambda ev : ev.Lep1Jet1_DR_jerDown),
            MVAVar("deltaRL1L2_J1J2_",      func = lambda ev : ev.Lep12Jet12_DR_jerDown),
            MVAVar("deltaRL1L2_metJ1J2_",   func = lambda ev : ev.Lep12Jet12MET_DR_jerDown),
        ]

        mvas_path = "/pool/phedex/userstorage/vrbouza/proyectos/tw_run2/mvas"

        path_SergioBDT_1j1b = mvas_path + "/2018_00_00_SergioBDT/bdt1j1t/weights/TMVAClassification_GradBoost_2000_0.01.weights.xml"
        path_SergioBDT_2j1b = mvas_path + "/2018_00_00_SergioBDT/bdt2j1t/weights/TMVAClassification_GradBoost_200_005_4.weights.xml"

        self.MVAs["SergioBDT_1j1b"]       = MVATool("BDT",         path_SergioBDT_1j1b, self.vars_1j1t)
        self.MVAs["SergioBDT_1j1b_JESUp"] = MVATool("BDT_JESUp",   path_SergioBDT_1j1b, self.vars_1j1t_jesup)
        self.MVAs["SergioBDT_1j1b_JESDn"] = MVATool("BDT_JESDown", path_SergioBDT_1j1b, self.vars_1j1t_jesdn)
        self.MVAs["SergioBDT_1j1b_JERUp"] = MVATool("BDT_JER",     path_SergioBDT_1j1b, self.vars_1j1t_jerup)
        self.MVAs["SergioBDT_1j1b_JERDn"] = MVATool("BDT_JER",     path_SergioBDT_1j1b, self.vars_1j1t_jerdn)
        self.MVAs["SergioBDT_2j1b"]       = MVATool("2j1b",        path_SergioBDT_2j1b, self.vars_2j1t)
        self.MVAs["SergioBDT_2j1b_JESUp"] = MVATool("2j1b",        path_SergioBDT_2j1b, self.vars_2j1t_jesup)
        self.MVAs["SergioBDT_2j1b_JESDn"] = MVATool("2j1b",        path_SergioBDT_2j1b, self.vars_2j1t_jesdn)
        self.MVAs["SergioBDT_2j1b_JERUp"] = MVATool("2j1b",        path_SergioBDT_2j1b, self.vars_2j1t_jerup)
        self.MVAs["SergioBDT_2j1b_JERDn"] = MVATool("2j1b",        path_SergioBDT_2j1b, self.vars_2j1t_jerdn)

        return


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.wrappedOutputTree = wrappedOutputTree
        for out in self.MVAs.keys():
            self.wrappedOutputTree.branch(out, 'F')
        return


    def analyze(self, event):
        writeOutput(self, dict([ (name, mva(event)) for name, mva in self.MVAs.iteritems()]))
        return True
