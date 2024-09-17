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
from CMGTools.TTHAnalysis.tools.nanoAOD.TopRun2UL_modules import ch, tags, emass

mresta_ = 172.5 ** 2 - 80.379 ** 2


class EventVars_TopRun2UL(Module):
    def __init__(self, label = "", recllabel = 'Recl', isMC = True,
                 jecvars = ["jesTotal", "jer"], lepvars = ["mu"]):

        self.jecbranches = ["Lep1Lep2Jet1MET_Pt",
                            "Lep1Lep2Jet1MET_M",
                            "Lep1Lep2Jet1MET_Mt",
                            "Lep1Lep2Jet1MET_PtOverHTtot",
                            "Lep1_PtLep2_PtOverHTtot",
                            "Lep1Lep2Jet1_Pt",
                            "Lep1Lep2Jet1_M",
                            "Lep1Lep2Jet1_E",
                            "Lep1Lep2Jet1_Pz",
                            "Lep1Jet1_Pt",
                            "Lep1Jet1_M",
                            "Lep1Jet1_DPhi",
                            "Lep2Jet1_Pt",
                            "Lep2Jet1_M",
                            "Lep2Jet1_DPhi",
                            "cosThetaStar1",
                            "cosThetaStar2",
                            "Jet1_E",
                            "Jet1_Pt",
                            "Jet2_Pt",
                            "JetLoose1_Pt",

                            "Lep1Jet1_DR",
                            "Lep12Jet12_DR",
                            "Lep12Jet12MET_DR",
                            "LepJet11Lep2_DR",
                            "LepJet11Lep2_DPhi",
                            "LepJet11Lep2_DEta",
                            "Lep1LepJet21_DR",
                            "Lep1LepJet21_DPhi",
                            "Lep1LepJet21_DEta",
                            "Lep1Lep2Jet1_C",
                            "Lep1Lep2Jet1_Cscalar",
                            "Lep1Lep2Jet1_CwithP",
                            "HTtot",

                            "minimax",
                            "Lep1Jet2_M",
                            "Lep2Jet2_M",

                            "METgood_pt",
                            "METgood_phi",]

        self.lepenergyvars = [("isSS", "I"),
                              ("channel", "I"),
                              "Lep1_Pt",
                              "Lep2_Pt",
                              "Lep1Lep2_Pt",
                              "Lep1Lep2_PtSum",
                              "Lep1Lep2_DR",
                              "Lep1Lep2_DPhi",
                              "Lep1Lep2_DEta",
                              "Mll",

                              "Lep1Lep2Jet1MET_Pt",
                              "Lep1Lep2Jet1MET_M",
                              "Lep1Lep2Jet1MET_Mt",
                              "Lep1Lep2Jet1MET_PtOverHTtot",
                              "Lep1_PtLep2_PtOverHTtot",
                              "Lep1Lep2Jet1_Pt",
                              "Lep1Lep2Jet1_M",
                              "Lep1Lep2Jet1_E",
                              "Lep1Lep2Jet1_Pz",
                              "Lep1Jet1_Pt",
                              "Lep1Jet1_M",
                              "Lep1Jet1_DPhi",
                              "Lep2Jet1_Pt",
                              "Lep2Jet1_M",
                              "Lep2Jet1_DPhi",
                              "cosThetaStar1",
                              "cosThetaStar2",
                              "Jet1_E",
                              "Jet1_Pt",
                              "Jet2_Pt",
                              "JetLoose1_Pt",

                              "Lep1Jet1_DR",
                              "Lep12Jet12_DR",
                              "Lep12Jet12MET_DR",
                              "LepJet11Lep2_DR",
                              "LepJet11Lep2_DPhi",
                              "LepJet11Lep2_DEta",
                              "Lep1LepJet21_DR",
                              "Lep1LepJet21_DPhi",
                              "Lep1LepJet21_DEta",
                              "Lep1Lep2Jet1_C",
                              "Lep1Lep2Jet1_Cscalar",
                              "Lep1Lep2Jet1_CwithP",
                              "HTtot",

                              "minimax",
                              "Lep1Jet2_M",
                              "Lep2Jet2_M",
                              "minMllAFAS",]

        self.otherbranches = [("isSS", "I"),
                              #("isSS", "O"),
                              #("channel", "B"),
                              ("channel", "I"),
                              ("datatag", "I"),
                              ("isTop", "I"),

                              #### Variables that are not susceptible to JEC variations that need a nominal value
                              # (note that being in the lepenergyvars list does not give this!)
                              "Lep1_Pt",
                              "Lep2_Pt",
                              "Lep1Lep2_Pt",
                              "Lep1Lep2_PtSum",
                              "Lep1Lep2_DR",
                              "Lep1Lep2_DPhi",
                              "Lep1Lep2_DEta",
                              "Mll",
                              "minMllAFAS",
                             ]


        self.label    = "" if (label in ["", None]) else ("_" + label)

        self.systsJEC   = {0: ""}
        self.systsLepEn = {}
        self.isMC       = isMC
        self.nominaljecscaff = "_nom"
        if not self.isMC:
            jecvars = []
            self.nominaljecscaff = ""
            if "elscale" in lepvars:
                lepvars = ["elscale"]
            else:
                lepvars = []

        if   len(jecvars):
            for i, var in enumerate(jecvars):
                self.systsJEC[i+1]    = "_%sUp"%var
                self.systsJEC[-(i+1)] = "_%sDown"%var

        if   len(lepvars):
            for i, var in enumerate(lepvars):
                self.systsLepEn[i+1]    = "_%sUp"%var
                self.systsLepEn[-(i+1)] = "_%sDown"%var

        self.inputlabel = '_' + recllabel
        self.branches   = []

        for delta,var in self.systsJEC.items():
            self.branches.extend([br + self.label + var for br in self.jecbranches])
        for delta,var in self.systsLepEn.items():
            self.branches.extend([el     + self.label + var for el in self.lepenergyvars if type(el) != tuple])
            self.branches.extend([(el[0] + self.label + var, el[1]) for el in self.lepenergyvars if type(el) == tuple])

        self.branches.extend(self.otherbranches)
        return


    # Old interface (CMG)
    def listBranches(self):
        return self.branches[:]


    def __call__(self, event):
        return self.run(event, CMGCollection, "met")


    # New interface (nanoAOD-tools)
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        declareOutput(self, wrappedOutputTree, self.branches)


    def analyze(self, event):
        writeOutput(self, self.run(event, NanoAODCollection))
        return True


    # Common processing
    def run(self, event, Collection):
        # ============================ Definitions and declarations
        allret = {}

        leps    = [l for l in Collection(event, "LepGood")]
        leps_4m = [l.p4() for l in leps]

        for i in range(len(leps_4m)):
            leps_4m[i].SetPtEtaPhiM(leps[i].pt_corrAll, leps_4m[i].Eta(),
                                    leps_4m[i].Phi(),   leps_4m[i].M() if abs(leps[i].pdgId) != 11 else emass)

        coldirs    = {} #### NOTE: the keys begin always with "_"!
        coldirs_4m = {}
        for delta,var in self.systsLepEn.items():
            coldirs[var]    = [l      for l in Collection(event, "LepGood" + var[1:])]
            coldirs_4m[var] = [l.p4() for l in coldirs[var]]
            for i in range(len(coldirs_4m[var])):
                coldirs_4m[var][i].SetPtEtaPhiM(getattr(coldirs[var][i], "pt" + var), coldirs_4m[var][i].Eta(),
                                                coldirs_4m[var][i].Phi(), coldirs_4m[var][i].M() if abs(coldirs[var][i].pdgId) != 11 else emass)

        all_jets = [j for j in Collection(event, "Jet")]
        jets     = []
        jets_4m  = []


        # ============================ Initialisations
        for var in self.jecbranches:
            for delta,sys in self.systsJEC.items():
                allret[var + sys] = -99
                allret[var + sys] = -99

        for prevar in self.lepenergyvars:
            var = prevar
            if type(prevar) ==  tuple:
                var = prevar[0]
            for delta,sys in self.systsLepEn.items():
                allret[var + sys] = -99
                allret[var + sys] = -99

        for var in self.otherbranches:
            if type(var) == tuple:
                if var[0] == "channel":
                    allret["channel"] = ch.NoChan
                elif var[0] == "isTop":
                    allret["isTop"] = 0
                elif var[0] == "datatag":
                    allret["datatag"] = tags.NoTag
                else:
                    allret[var[0]] = -99
            else:
                allret[var] = -99
                allret[var] = -99

        # ============================ Calculations
        thesample = ""
        for i in range(event.nDatasetName):
            thesample += str(event.DatasetName_name[i])

        if event.isData:
            if   "muoneg"     in thesample.lower():
                allret["datatag"] = tags.muoneg
            elif "singlemuon" in thesample.lower():
                allret["datatag"] = tags.singlemuon
            elif "singleelec" in thesample.lower() or "egamma" in thesample.lower():
                allret["datatag"] = tags.singleelec
            elif "doublemuon" in thesample.lower():
                allret["datatag"] = tags.doublemuon
            elif "doubleeg" in thesample.lower():
                allret["datatag"] = tags.doubleeg
        elif "ttto2l2nu" in thesample.lower():
            allret["datatag"] = tags.mc
            allret["isTop"]   = 1
        else:
            allret["datatag"] = tags.mc


        # Nominal variables and those susceptible to JEC enter in this if-statement.
        for iL in range(len(leps)):
            for jL in range(iL + 1, len(leps)):
                tmpM = (leps_4m[iL] + leps_4m[jL]).M()
                if tmpM < allret["minMllAFAS"] or allret["minMllAFAS"] == -99: allret["minMllAFAS"] = tmpM

        if event.nLepGood >= 2:
            # ============================ Nominal variables
            allret["isSS"] = int(leps[0].charge == leps[1].charge)
            if   ((abs(leps[0].pdgId) == 13 and abs(leps[1].pdgId) == 11) or
                 (abs(leps[0].pdgId) == 11 and abs(leps[1].pdgId) == 13)):
                allret["channel"] = ch.ElMu
            elif (abs(leps[0].pdgId) == 13 and abs(leps[1].pdgId) == 13):
                allret["channel"] = ch.Muon
            elif (abs(leps[0].pdgId) == 11 and abs(leps[1].pdgId) == 11):
                allret["channel"] = ch.Elec
            else:
                allret["channel"] = ch.NoChan

            allret["Lep1_Pt"]        = leps_4m[0].Pt()
            allret["Lep2_Pt"]        = leps_4m[1].Pt()
            allret["Lep1Lep2_Pt"]    = (leps_4m[0] + leps_4m[1]).Pt()
            allret["Lep1Lep2_PtSum"] = leps_4m[0].Pt() + leps_4m[1].Pt()
            allret["Lep1Lep2_DR"]    = leps_4m[0].DeltaR(leps_4m[1])
            allret["Lep1Lep2_DPhi"]  = abs(deltaPhi(leps[0], leps[1]))/r.TMath.Pi()
            allret["Lep1Lep2_DEta"]  = abs(leps_4m[0].Eta() - leps_4m[1].Eta())
            allret["Mll"]            = (leps_4m[0] + leps_4m[1]).M()


            # ============================ Variables susceptible to JEC (nominal and their variations)
            met_4m = r.TLorentzVector()
            for delta,sys in self.systsJEC.items():
                metpt  = -99.
                metphi = -99.
                if event.isData:
                    metpt  = event.MET_pt
                    metphi = event.MET_phi
                else:
                    isJER = False
                    for JERSource in ["jerUp","jer0Up","jer1Up","jer2Up","jer3Up","jer4Up","jer5Up","jerDown","jer0Down","jer1Down","jer2Down","jer3Down","jer4Down","jer5Down"]:
                        if JERSource in sys:
                            isJER = True

                    if isJER: ###### Type I corrected MET should not have JER variations propagated
                        metpt  = getattr(event, 'MET_T1_pt')
                        metphi = getattr(event, 'MET_T1_phi')
                    else:
                        metpt  = getattr(event, 'MET_T1_pt{v}'.format( v = sys if sys != "" else ""))
                        metphi = getattr(event, 'MET_T1_phi{v}'.format(v = sys if sys != "" else ""))

                met_4m.SetPtEtaPhiM(metpt, 0, metphi, 0)

                ### jets
                jets    = [all_jets[getattr(event, 'iJetSel30{v}_Recl'.format(v = sys if "unclustEn" not in sys else ""))[j]]
                           for j in range(min([getattr(event, 'nJetSel30{v}_Recl'.format(v = sys if "unclustEn" not in sys else "")), 5]))]
                jets_4m = [j.p4() for j in jets]

                jetjecsysscaff = (sys if (sys != "" and "unclustEn" not in sys) else self.nominaljecscaff)
                if event.isData: jetjecsysscaff = ""

                for i in range(len(jets_4m)):
                    jets_4m[i].SetPtEtaPhiM(getattr(jets[i], "pt" + jetjecsysscaff), jets_4m[i].Eta(),
                                            jets_4m[i].Phi(), getattr(jets[i], "mass" + jetjecsysscaff))

                ### loose jets
                loosejets = [all_jets[getattr(event, 'iJetSel20{v}_Recl'.format(v = sys if "unclustEn" not in sys else ""))[j]]
                           for j in range(min([getattr(event, 'nJetSel20{v}_Recl'.format(v = sys if "unclustEn" not in sys else "")), 5]))]
                loosejets_4m = [j.p4() for j in loosejets]

                for i in range(len(loosejets_4m)):
                    loosejets_4m[i].SetPtEtaPhiM(getattr(loosejets[i], "pt" + jetjecsysscaff), loosejets_4m[i].Eta(),
                                            loosejets_4m[i].Phi(), getattr(loosejets[i], "mass" + jetjecsysscaff))


                if getattr(event, 'nJetSel20{v}_Recl'.format(v = sys if "unclustEn" not in sys else "")) > 0:
                    allret["JetLoose1_Pt" + sys] = loosejets_4m[0].Pt()

                allret["METgood_pt"  + sys] = met_4m.Pt()
                allret["METgood_phi" + sys] = met_4m.Phi()

                if getattr(event, 'nJetSel30{v}_Recl'.format(v = sys if "unclustEn" not in sys else "")) > 0:
                    allret["Lep1Lep2Jet1MET_Pt"          + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0] + met_4m).Pt()
                    allret["Lep1Lep2Jet1MET_M"           + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0] + met_4m).M()
                    allret["Lep1Lep2Jet1MET_Mt"          + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0] + met_4m).Mt()
                    allret["Lep1Lep2Jet1_Pt"             + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0]).Pt()
                    allret["Lep1Lep2Jet1_M"              + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0]).M()
                    allret["Lep1Lep2Jet1_E"              + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0]).E()
                    allret["Lep1Lep2Jet1_Pz"             + sys] = abs((leps_4m[0] + leps_4m[1] + jets_4m[0]).Pz())
                    allret["Lep1Jet1_Pt"                 + sys] = (leps_4m[0] + jets_4m[0]).Pt()
                    allret["Lep1Jet1_M"                  + sys] = (leps_4m[0] + jets_4m[0]).M()
                    allret["Lep1Jet1_DPhi"               + sys] = abs(leps_4m[0].Phi() - jets_4m[0].Phi()) / r.TMath.Pi()
                    allret["cosThetaStar1"               + sys] = 2 * (allret["Lep1Jet1_M" + sys] ** (2) / mresta_) - 1
                    allret["Lep2Jet1_Pt"                 + sys] = (leps_4m[1] + jets_4m[0]).Pt()
                    allret["Lep2Jet1_M"                  + sys] = (leps_4m[1] + jets_4m[0]).M()
                    allret["Lep2Jet1_DPhi"               + sys] = abs(leps_4m[1].Phi() - jets_4m[0].Phi()) / r.TMath.Pi()
                    allret["cosThetaStar2"               + sys] = 2 * (allret["Lep2Jet1_M" + sys] ** (2) / mresta_) - 1
                    allret["Jet1_Pt"                     + sys] = jets_4m[0].Pt()
                    allret["Jet1_E"                      + sys] = jets_4m[0].E()

                    allret["Lep1Jet1_DR"                 + sys] = leps_4m[0].DeltaR(jets_4m[0])
                    allret["Lep1Lep2Jet1_C"              + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0]).Et() / (leps_4m[0] + leps_4m[1] + jets_4m[0]).E()
                    allret["Lep1Lep2Jet1_Cscalar"        + sys] = (leps_4m[0].Et() + leps_4m[1].Et() + jets_4m[0].Et()) / (leps_4m[0].E() + leps_4m[1].E() + jets_4m[0].E())
                    allret["Lep1Lep2Jet1_CwithP"         + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0]).Pt() / (leps_4m[0] + leps_4m[1] + jets_4m[0]).P()
                    allret["HTtot"                       + sys] = leps_4m[0].Pt() + leps_4m[1].Pt() + jets_4m[0].Pt() + met_4m.Pt()
                    allret["Lep1Lep2Jet1MET_PtOverHTtot" + sys] = allret["Lep1Lep2Jet1MET_Pt" + sys] / allret["HTtot" + sys]
                    allret["Lep1_PtLep2_PtOverHTtot"     + sys] = (leps_4m[0].Pt() + leps_4m[1].Pt()) / allret["HTtot" + sys]

                    allret["LepJet11Lep2_DR"            + sys] = (leps_4m[0] + jets_4m[0]).DeltaR(leps_4m[1])
                    allret["LepJet11Lep2_DPhi"          + sys] = abs((leps_4m[0] + jets_4m[0]).Phi() - leps_4m[1].Phi()) / r.TMath.Pi()
                    allret["LepJet11Lep2_DEta"          + sys] = abs((leps_4m[0] + jets_4m[0]).Eta() - leps_4m[1].Eta())
                    allret["Lep1LepJet21_DR"            + sys] = leps_4m[0].DeltaR(leps_4m[1] + jets_4m[0])
                    allret["Lep1LepJet21_DPhi"          + sys] = abs(leps_4m[0].Phi() - (jets_4m[0] + leps_4m[1]).Phi()) / r.TMath.Pi()
                    allret["Lep1LepJet21_DEta"          + sys] = abs(leps_4m[0].Phi() - (jets_4m[0] + leps_4m[1]).Eta())


                    if getattr(event, 'nJetSel30{v}_Recl'.format(v = sys if "unclustEn" not in sys else "")) > 1:
                        allret["Lep12Jet12_DR"    + sys] = (leps_4m[0] + leps_4m[1]).DeltaR(jets_4m[0] + jets_4m[1])
                        allret["Lep12Jet12MET_DR" + sys] = (leps_4m[0] + leps_4m[1]).DeltaR(jets_4m[0] + jets_4m[1] + met_4m)
                        allret["Lep1Jet2_M"       + sys] = (leps_4m[0] + jets_4m[1]).M()
                        allret["Lep2Jet2_M"       + sys] = (leps_4m[1] + jets_4m[1]).M()
                        allret["Jet2_Pt"          + sys] = jets_4m[1].Pt()
                        #### WARNING: this minimax variable is only equal to that of ATLAS' PRL when the signal region is njets == 2, nbjets == 2.
                        allret["minimax"          + sys] = min([max([allret["Lep1Jet1_M" + sys],
                                                                     allret["Lep2Jet2_M" + sys]]),
                                                                max( allret["Lep2Jet1_M" + sys],
                                                                     allret["Lep1Jet2_M" + sys])])



        # ============================ Variations of lepton energy corrections
        metpt  = -99.
        metphi = -99.
        if event.isData:
            metpt  = event.MET_pt
            metphi = event.MET_phi
        else:
            metpt  = getattr(event, 'MET_T1_pt')
            metphi = getattr(event, 'MET_T1_phi')

        met_4m = r.TLorentzVector()
        met_4m.SetPtEtaPhiM(metpt,  0, metphi, 0)

        for delta,sys in self.systsLepEn.items():
            leps    = [l for l in Collection(event, "LepGood" + sys[1:])]
            leps_4m = [l.p4() for l in leps]

            for i in range(len(leps_4m)):
                leps_4m[i].SetPtEtaPhiM(getattr(leps[i], "pt" + sys), leps_4m[i].Eta(),
                                        leps_4m[i].Phi(),             leps_4m[i].M() if abs(leps[i].pdgId) != 11 else emass)

            for iL in range(len(leps)):
                for jL in range(iL + 1, len(leps)):
                    tmpM = (leps_4m[iL] + leps_4m[jL]).M()
                    if tmpM < allret["minMllAFAS" + sys] or allret["minMllAFAS" + sys] == -99: allret["minMllAFAS" + sys] = tmpM

            if getattr(event, "nLepGood" + sys[1:]) >= 2:
                jets    = [all_jets[getattr(event, 'iJetSel30{v}_Recl'.format(v = sys))[j]]
                           for j in range(min([getattr(event, 'nJetSel30{v}_Recl'.format(v = sys)), 5]))]
                jets_4m = [j.p4() for j in jets]

                jetjecsysscaff = self.nominaljecscaff
                if event.isData: jetjecsysscaff = ""

                for i in range(len(jets_4m)):
                    jets_4m[i].SetPtEtaPhiM(getattr(jets[i], "pt" + jetjecsysscaff), jets_4m[i].Eta(),
                                            jets_4m[i].Phi(), getattr(jets[i], "mass" + jetjecsysscaff))

                ### loose jets
                loosejets = [all_jets[getattr(event, 'iJetSel20{v}_Recl'.format(v = sys))[j]]
                           for j in range(min([getattr(event, 'nJetSel20{v}_Recl'.format(v = sys)), 5]))]
                loosejets_4m = [j.p4() for j in loosejets]

                for i in range(len(loosejets_4m)):
                    loosejets_4m[i].SetPtEtaPhiM(getattr(loosejets[i], "pt" + jetjecsysscaff), loosejets_4m[i].Eta(),
                                            loosejets_4m[i].Phi(), getattr(loosejets[i], "mass" + jetjecsysscaff))

                # ============================ Nominal variables
                allret["isSS" + sys] = int(leps[0].charge == leps[1].charge)
                if   ((abs(leps[0].pdgId) == 13 and abs(leps[1].pdgId) == 11) or
                    (abs(leps[0].pdgId) == 11 and abs(leps[1].pdgId) == 13)):
                    allret["channel" + sys] = ch.ElMu
                elif (abs(leps[0].pdgId) == 13 and abs(leps[1].pdgId) == 13):
                    allret["channel" + sys] = ch.Muon
                elif (abs(leps[0].pdgId) == 11 and abs(leps[1].pdgId) == 11):
                    allret["channel" + sys] = ch.Elec
                else:
                    allret["channel" + sys] = ch.NoChan

                allret["Lep1_Pt"        + sys] = leps_4m[0].Pt()
                allret["Lep2_Pt"        + sys] = leps_4m[1].Pt()
                allret["Lep1Lep2_Pt"    + sys] = (leps_4m[0] + leps_4m[1]).Pt()
                allret["Lep1Lep2_PtSum" + sys] = leps_4m[0].Pt() + leps_4m[1].Pt()
                allret["Lep1Lep2_DR"    + sys] = leps_4m[0].DeltaR(leps_4m[1])
                allret["Lep1Lep2_DPhi"  + sys] = abs(deltaPhi(leps[0], leps[1]))/r.TMath.Pi()
                allret["Lep1Lep2_DEta"  + sys] = abs(leps_4m[0].Eta() - leps_4m[1].Eta())
                allret["Mll"            + sys] = (leps_4m[0] + leps_4m[1]).M()

                # ============================ Variables susceptible also to JEC
                if getattr(event, 'nJetSel20{v}_Recl'.format(v = sys)) > 0:
                    allret["JetLoose1_Pt" + sys] = loosejets_4m[0].Pt()

                if getattr(event, 'nJetSel30{v}_Recl'.format(v = sys)) > 0:
                    allret["Lep1Lep2Jet1MET_Pt"          + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0] + met_4m).Pt()
                    allret["Lep1Lep2Jet1MET_M"           + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0] + met_4m).M()
                    allret["Lep1Lep2Jet1MET_Mt"          + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0] + met_4m).Mt()
                    allret["Lep1Lep2Jet1_Pt"             + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0]).Pt()
                    allret["Lep1Lep2Jet1_M"              + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0]).M()
                    allret["Lep1Lep2Jet1_E"              + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0]).E()
                    allret["Lep1Lep2Jet1_Pz"             + sys] = abs((leps_4m[0] + leps_4m[1] + jets_4m[0]).Pz())
                    allret["Lep1Jet1_Pt"                 + sys] = (leps_4m[0] + jets_4m[0]).Pt()
                    allret["Lep1Jet1_M"                  + sys] = (leps_4m[0] + jets_4m[0]).M()
                    allret["Lep1Jet1_DPhi"               + sys] = abs(leps_4m[0].Phi() - jets_4m[0].Phi()) / r.TMath.Pi()
                    allret["cosThetaStar1"               + sys] = 2 * (allret["Lep1Jet1_M" + sys] ** (2) / mresta_) - 1
                    allret["Lep2Jet1_Pt"                 + sys] = (leps_4m[1] + jets_4m[0]).Pt()
                    allret["Lep2Jet1_M"                  + sys] = (leps_4m[1] + jets_4m[0]).M()
                    allret["Lep2Jet1_DPhi"               + sys] = abs(leps_4m[1].Phi() - jets_4m[0].Phi()) / r.TMath.Pi()
                    allret["cosThetaStar2"               + sys] = 2 * (allret["Lep2Jet1_M" + sys] ** (2) / mresta_) - 1
                    allret["Jet1_Pt"                     + sys] = jets_4m[0].Pt()
                    allret["Jet1_E"                      + sys] = jets_4m[0].E()

                    allret["Lep1Jet1_DR"                 + sys] = leps_4m[0].DeltaR(jets_4m[0])
                    allret["Lep1Lep2Jet1_C"              + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0]).Et() / (leps_4m[0] + leps_4m[1] + jets_4m[0]).E()
                    allret["Lep1Lep2Jet1_Cscalar"        + sys] = (leps_4m[0].Et() + leps_4m[1].Et() + jets_4m[0].Et()) / (leps_4m[0].E() + leps_4m[1].E() + jets_4m[0].E())
                    allret["Lep1Lep2Jet1_CwithP"         + sys] = (leps_4m[0] + leps_4m[1] + jets_4m[0]).Pt() / (leps_4m[0] + leps_4m[1] + jets_4m[0]).P()

                    allret["LepJet11Lep2_DR"             + sys] = (leps_4m[0] + jets_4m[0]).DeltaR(leps_4m[1])
                    allret["LepJet11Lep2_DPhi"           + sys] = abs((leps_4m[0] + jets_4m[0]).Phi() - leps_4m[1].Phi()) / r.TMath.Pi()
                    allret["LepJet11Lep2_DEta"           + sys] = abs((leps_4m[0] + jets_4m[0]).Eta() - leps_4m[1].Eta())
                    allret["Lep1LepJet21_DR"             + sys] = leps_4m[0].DeltaR(leps_4m[1] + jets_4m[0])
                    allret["Lep1LepJet21_DPhi"           + sys] = abs(leps_4m[0].Phi() - (jets_4m[0] + leps_4m[1]).Phi()) / r.TMath.Pi()
                    allret["Lep1LepJet21_DEta"           + sys] = abs(leps_4m[0].Phi() - (jets_4m[0] + leps_4m[1]).Eta())

                    allret["HTtot"                       + sys] = leps_4m[0].Pt() + leps_4m[1].Pt() + jets_4m[0].Pt() + met_4m.Pt()
                    allret["Lep1Lep2Jet1MET_PtOverHTtot" + sys] = allret["Lep1Lep2Jet1MET_Pt" + sys] / allret["HTtot" + sys]
                    allret["Lep1_PtLep2_PtOverHTtot"     + sys] = (leps_4m[0].Pt() + leps_4m[1].Pt()) / allret["HTtot" + sys]

                    if getattr(event, 'nJetSel30{v}_Recl'.format(v = sys)) > 1:
                        allret["Lep12Jet12_DR"    + sys] = (leps_4m[0] + leps_4m[1]).DeltaR(jets_4m[0] + jets_4m[1])
                        allret["Lep12Jet12MET_DR" + sys] = (leps_4m[0] + leps_4m[1]).DeltaR(jets_4m[0] + jets_4m[1] + met_4m)
                        allret["Lep1Jet2_M"       + sys] = (leps_4m[0] + jets_4m[1]).M()
                        allret["Lep2Jet2_M"       + sys] = (leps_4m[1] + jets_4m[1]).M()
                        allret["Jet2_Pt"          + sys] = jets_4m[1].Pt()
                        #### WARNING: this minimax variable is only equal to that of ATLAS' PRL when the signal region is njets == 2, nbjets == 2.
                        allret["minimax"          + sys] = min([max([allret["Lep1Jet1_M" + sys],
                                                                     allret["Lep2Jet2_M" + sys]]),
                                                                max( allret["Lep2Jet1_M" + sys],
                                                                     allret["Lep1Jet2_M" + sys])])
        return allret
