from math import sqrt, cos
from copy import deepcopy
import struct as st
import warnings as wr
import ROOT as r

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection as NanoAODCollection 
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR,deltaPhi

from CMGTools.TTHAnalysis.treeReAnalyzer import Collection as CMGCollection
from CMGTools.TTHAnalysis.tools.nanoAOD.friendVariableProducerTools import declareOutput, writeOutput
from CMGTools.TTHAnalysis.tools.nanoAOD.TopRun2_modules import ch, tags


class createTrainingMiniTree_TopRun3(Module):
    def __init__(self):
        self.branches = [("njets",   "I"),
                         ("nbjets",  "I"),
                         ("channel", "I"),
                         ("train_nloosejets",  "I"),
                         ("train_nbloosejets", "I"),
                         "train_jet1_pt",
                         "train_loosejet1_pt",
                         "train_lep1lep2jet1_pt",
                         "train_lep1lep2jet1met_pt",
                         "train_lep1lep2jet1met_ptOVERhttot",
                         "train_lep1lep2jet1met_m",
                         "train_htlepOVERhttot",
                         "train_lep1lep2jet1_c",
                         "train_lep1lep2jet1_cscalar",
                         "train_httot",

                         "train_jet2_pt",
                         "train_lep1jet1_dr",
                         "train_lep12jet12_dr",
                         "train_lep12jet12met_dr",

                         #### NUEVAS
                         "train_lep1lep2jet1met_mt",
                         "train_lep1lep2jet1_pz",
                         "train_lep1lep2jet1_m",
                         "train_lep1lep2jet1_e",
                         "train_lep1jet1_pt",
                         "train_lep1jet1_m",
                         "train_lep1jet2_m",
                         "train_lep2jet1_pt",
                         "train_lep2jet1_m",
                         "train_lep1jet2_m",
                         "train_lepjet11lep2_dr",
                         "train_lepjet11lep2_dphi",
                         "train_lepjet11lep2_deta",
                         "train_lep1jetjet21_dr",
                         "train_lep1jetjet21_dphi",
                         "train_lep1jetjet21_deta",
                         "train_met_pt",
                         "train_met_phi",
                         ("train_nfwdjet",  "I"),
                         ("train_nloosefwdjet",  "I"),

                         "train_lep1lep2_pt",
                         "train_lep1lep2_ptsum",
                         "train_lep1lep2_dr",
                         "train_lep1lep2_dphi",
                         "train_lep1lep2_deta",
                         "train_lep1lep2_m",

                         "train_lep1_pt",
                         "train_lep1_eta",
                         "train_lep1_phi",
                         "train_lep1_m",
                         "train_lep2_pt",
                         "train_lep2_eta",
                         "train_lep2_phi",
                         "train_lep2_m",

                         #"traiW_jet1_pt", # esta antes
                         "train_jet1_eta",
                         "train_jet1_phi",
                         "train_jet1_m",
                         #"traiW_jet2_pt", # esta antes
                         "train_jet2_eta",
                         "train_jet2_phi",
                         "train_jet2_m",
                         "train_minimax",

                         #"traiW_loosejet1_pt", # esta antes
                         "train_loosejet1_eta",
                         "train_loosejet1_phi",
                         "train_loosejet1_m",

                         "train_fwdjet1_pt",
                         "train_fwdjet1_eta",
                         "train_fwdjet1_phi",
                         "train_fwdjet1_m",

                         "train_fwdloosejet1_pt",
                         "train_fwdloosejet1_eta",
                         "train_fwdloosejet1_phi",
                         "train_fwdloosejet1_m",

                         "allweights"]
        return


    # New interface (nanoAOD-tools)
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        declareOutput(self, wrappedOutputTree, self.branches)
        return


    def analyze(self, event):
        if event.pass_trigger:
            if event.Flag_goodVertices and event.Flag_globalSuperTightHalo2016Filter and event.Flag_HBHENoiseFilter and event.Flag_HBHENoiseIsoFilter and event.Flag_EcalDeadCellTriggerPrimitiveFilter and event.Flag_BadPFMuonFilter:
                if event.nLepGood >= 2 and event.isSS == 0 and event.nJetSel30_Recl > 0:
                    if event.LepGood_pt_corrAll[0] > 25 and event.LepGood_pt_corrAll[1] > 20 and event.minMllAFAS > 20:
                        #if event.nGenDressedLepton >= 2: ## Why this?
                        writeOutput(self, self.run(event, NanoAODCollection))
                        return True
        #if event.separationIndex == 0:
            #writeOutput(self, self.run(event, NanoAODCollection))
            #return True
        return False


    # Common processing
    def run(self, event, Collection):
        allret = {}

        allret["njets"]   = event.nJetSel30_Recl
        allret["nbjets"]  = event.nBJetSelMedium30_Recl
        allret["channel"] = event.channel

        allret["train_nloosejets"                 ] = event.nJetSel20_Recl
        allret["train_nbloosejets"                ] = event.nBJetSelMedium20_Recl
        allret["train_jet1_pt"                    ] = event.Jet1_Pt                     if event.Jet1_Pt                     >= 0 else 0
        allret["train_loosejet1_pt"               ] = event.JetLoose1_Pt                if event.JetLoose1_Pt                >= 0 else 0
        allret["train_lep1lep2jet1_pt"            ] = event.Lep1Lep2Jet1_Pt             if event.Lep1Lep2Jet1_Pt             >= 0 else 0
        allret["train_lep1lep2jet1met_pt"         ] = event.Lep1Lep2Jet1MET_Pt          if event.Lep1Lep2Jet1MET_Pt          >= 0 else 0
        allret["train_lep1lep2jet1met_ptOVERhttot"] = event.Lep1Lep2Jet1MET_PtOverHTtot if event.Lep1Lep2Jet1MET_PtOverHTtot >= 0 else 0
        allret["train_lep1lep2jet1met_m"          ] = event.Lep1Lep2Jet1MET_M           if event.Lep1Lep2Jet1MET_M           >= 0 else 0
        allret["train_htlepOVERhttot"             ] = event.Lep1_PtLep2_PtOverHTtot     if event.Lep1_PtLep2_PtOverHTtot     >= 0 else 0
        allret["train_lep1lep2jet1_c"             ] = event.Lep1Lep2Jet1_C              if event.Lep1Lep2Jet1_C              >= 0 else 0
        allret["train_lep1lep2jet1_cscalar"       ] = event.Lep1Lep2Jet1_Cscalar        if event.Lep1Lep2Jet1_Cscalar        >= 0 else 0
        allret["train_httot"                      ] = event.HTtot                       if event.HTtot                       >= 0 else 0

        allret["train_jet2_pt"                    ] = event.Jet2_Pt                     if event.Jet2_Pt                     >= 0 else 0
        allret["train_lep1jet1_dr"                ] = event.Lep1Jet1_DR                 if event.Lep1Jet1_DR                 >= 0 else 0
        allret["train_lep12jet12_dr"              ] = event.Lep12Jet12_DR               if event.Lep12Jet12_DR               >= 0 else 0
        allret["train_lep12jet12met_dr"           ] = event.Lep12Jet12MET_DR            if event.Lep12Jet12MET_DR            >= 0 else 0
        allret["train_minimax"                    ] = event.minimax                     if event.minimax                     >= 0 else 0

        #### NUEVAS
        allret["train_lep1lep2jet1met_mt"         ] = event.Lep1Lep2Jet1MET_Mt          if event.Lep1Lep2Jet1MET_Mt          >= 0 else 0
        allret["train_lep1lep2jet1_pz"            ] = event.Lep1Lep2Jet1_Pz             if event.Lep1Lep2Jet1_Pz             >= 0 else 0
        allret["train_lep1lep2jet1_m"             ] = event.Lep1Lep2Jet1_M              if event.Lep1Lep2Jet1_M              >= 0 else 0
        allret["train_lep1lep2jet1_e"             ] = event.Lep1Lep2Jet1_E              if event.Lep1Lep2Jet1_E              >= 0 else 0
        allret["train_lep1jet1_pt"                ] = event.Lep1Jet1_Pt                 if event.Lep1Jet1_Pt                 >= 0 else 0
        allret["train_lep1jet1_m"                 ] = event.Lep1Jet1_M                  if event.Lep1Jet1_M                  >= 0 else 0
        allret["train_lep1jet2_m"                 ] = event.Lep1Jet2_M                  if event.Lep1Jet2_M                  >= 0 else 0
        allret["train_lep2jet1_pt"                ] = event.Lep2Jet1_Pt                 if event.Lep2Jet1_Pt                 >= 0 else 0
        allret["train_lep2jet1_m"                 ] = event.Lep2Jet1_M                  if event.Lep2Jet1_M                  >= 0 else 0
        allret["train_lep1jet2_m"                 ] = event.Lep1Jet2_M                  if event.Lep1Jet2_M                  >= 0 else 0
        allret["train_lepjet11lep2_dr"            ] = event.LepJet11Lep2_DR             if event.LepJet11Lep2_DR             >= 0 else 0
        allret["train_lepjet11lep2_dphi"          ] = event.LepJet11Lep2_DPhi           if event.LepJet11Lep2_DPhi           >= 0 else 0
        allret["train_lepjet11lep2_deta"          ] = event.LepJet11Lep2_DEta           if event.LepJet11Lep2_DEta           >= 0 else 0
        allret["train_lep1jetjet21_dr"            ] = event.Lep1LepJet21_DR             if event.Lep1LepJet21_DR             >= 0 else 0
        allret["train_lep1jetjet21_dphi"          ] = event.Lep1LepJet21_DPhi           if event.Lep1LepJet21_DPhi           >= 0 else 0
        allret["train_lep1jetjet21_deta"          ] = event.Lep1LepJet21_DEta           if event.Lep1LepJet21_DEta           >= 0 else 0
        allret["train_met_pt"                     ] = event.METgood_pt
        allret["train_met_phi"                    ] = event.METgood_phi
        allret["train_nfwdjet"                    ] = event.nFwdJetSel30_Recl
        allret["train_nloosefwdjet"               ] = event.nFwdJetSel20_Recl

        allret["train_lep1lep2_pt"                ] = event.Lep1Lep2_Pt                 if event.Lep1Lep2_Pt                 >= 0 else 0
        allret["train_lep1lep2_ptsum"             ] = event.Lep1Lep2_PtSum              if event.Lep1Lep2_PtSum              >= 0 else 0
        allret["train_lep1lep2_dr"                ] = event.Lep1Lep2_DR                 if event.Lep1Lep2_DR                 >= 0 else 0
        allret["train_lep1lep2_dphi"              ] = event.Lep1Lep2_DPhi               if event.Lep1Lep2_DPhi               >= 0 else 0
        allret["train_lep1lep2_deta"              ] = event.Lep1Lep2_DEta               if event.Lep1Lep2_DEta               >= 0 else 0
        allret["train_lep1lep2_m"                 ] = event.Mll                         if event.Mll                         >= 0 else 0

        allret["train_lep1_pt"                    ] = event.LepGood_pt_corrAll[0]
        allret["train_lep1_eta"                   ] = event.LepGood_eta[0]
        allret["train_lep1_phi"                   ] = event.LepGood_phi[0]
        allret["train_lep1_m"                     ] = event.LepGood_mass[0]
        allret["train_lep2_pt"                    ] = event.LepGood_pt_corrAll[1]
        allret["train_lep2_eta"                   ] = event.LepGood_eta[1]
        allret["train_lep2_phi"                   ] = event.LepGood_phi[1]
        allret["train_lep2_m"                     ] = event.LepGood_mass[1]

        if event.nJetSel30_Recl > 0:
            tmpind = event.iJetSel30_Recl[0] #### WARNING: desconozco por que, pero parez que es importante entrar de antemano en esta variable, o en una de ellas de esta forma; si no, peta. PORSIACASO, replico esto en los siguientes condicionales, pero aparentemente funciona con este solo y dejando el resto con el ijet...[1]
            #allret["train_jet1_pt" # esta antes      ] = event.
            allret["train_jet1_eta"                   ] = event.Jet_eta[tmpind]
            allret["train_jet1_phi"                   ] = event.Jet_phi[tmpind]
            allret["train_jet1_m"                     ] = event.Jet_mass_nom[tmpind]

            if event.nJetSel30_Recl > 1:
                tmpind = event.iJetSel30_Recl[1]
                #allret["train_jet2_pt" # esta antes      ] = event.   if event.  >= 0 else 0
                allret["train_jet2_eta"                   ] = event.Jet_eta[tmpind]
                allret["train_jet2_phi"                   ] = event.Jet_phi[tmpind]
                allret["train_jet2_m"                     ] = event.Jet_mass_nom[tmpind]
            else:
                #allret["train_jet2_pt" # esta antes      ] = 0
                allret["train_jet2_eta"                   ] = 0
                allret["train_jet2_phi"                   ] = 0
                allret["train_jet2_m"                     ] = 0
        else:
            #allret["train_jet1_pt" # esta antes      ] = 0
            allret["train_jet1_eta"                   ] = 0
            allret["train_jet1_phi"                   ] = 0
            allret["train_jet1_m"                     ] = 0

            #allret["train_jet2_pt" # esta antes      ] = 0
            allret["train_jet2_eta"                   ] = 0
            allret["train_jet2_phi"                   ] = 0
            allret["train_jet2_m"                     ] = 0

        if event.nJetSel20_Recl > 0:
            tmpind = event.iJetSel20_Recl[0]
            #allret["train_loosejet1_pt" # esta antes ] = event.   if event.  >= 0 else 0
            allret["train_loosejet1_eta"              ] = event.Jet_eta[tmpind]
            allret["train_loosejet1_phi"              ] = event.Jet_phi[tmpind]
            allret["train_loosejet1_m"                ] = event.Jet_mass_nom[tmpind]
        else:
            #allret["train_loosejet1_pt" # esta antes ] = 0
            allret["train_loosejet1_eta"              ] = 0
            allret["train_loosejet1_phi"              ] = 0
            allret["train_loosejet1_m"                ] = 0

        if event.nFwdJetSel30_Recl > 0:
            tmpind = event.iFwdJetSel30_Recl[0]
            allret["train_fwdjet1_pt"                 ] = event.Jet_pt_nom[tmpind]
            allret["train_fwdjet1_eta"                ] = event.Jet_eta[tmpind]
            allret["train_fwdjet1_phi"                ] = event.Jet_phi[tmpind]
            allret["train_fwdjet1_m"                  ] = event.Jet_mass_nom[tmpind]
        else:
            allret["train_fwdjet1_pt"                 ] = 0
            allret["train_fwdjet1_eta"                ] = 0
            allret["train_fwdjet1_phi"                ] = 0
            allret["train_fwdjet1_m"                  ] = 0

        if event.nFwdJetSel20_Recl > 0:
            tmpind = event.iFwdJetSel30_Recl[1]
            allret["train_fwdloosejet1_pt"            ] = event.Jet_pt_nom[tmpind]
            allret["train_fwdloosejet1_eta"           ] = event.Jet_eta[tmpind]
            allret["train_fwdloosejet1_phi"           ] = event.Jet_phi[tmpind]
            allret["train_fwdloosejet1_m"             ] = event.Jet_mass_nom[tmpind]
        else:
            allret["train_fwdloosejet1_pt"            ] = 0
            allret["train_fwdloosejet1_eta"           ] = 0
            allret["train_fwdloosejet1_phi"           ] = 0
            allret["train_fwdloosejet1_m"             ] = 0

        # Do we need to add the weights?
        #allret["allweights"]                        = (event.genWeight * event.MuonIDSF * event.ElecIDSF *
        #                                               event.TrigSF * event.puWeight * event.bTagWeight)

        return allret
