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
from CMGTools.TTHAnalysis.tools.nanoAOD.TopRun2UL_modules import ch, tags


class particleAndPartonVars_TopRun2UL(Module):
    def __init__(self):
        self.wmass = 80.379
        self.bmass = 4.18

        self.branches = [("DressisSS", "I"),
                         #("isSS", "O"),
                         #("channel", "B"),
                         ("Dresschannel", "I"),
                         ("DressEXTchannel", "I"),
                         ("Origchannel", "I"),
                         "DressLep1_Pt",
                         "DressLep2_Pt",
                         "DressLep1Lep2_Pt",
                         "DressLep1Lep2_DPhi",
                         "DressMll",
                         "DressminMllAFAS",
                         "DressLep1Lep2Jet1MET_Pz",
                         "DressLep1Lep2Jet1MET_Pt",
                         "DressLep1Lep2Jet1MET_M",
                         "DressLep1Lep2Jet1MET_Mt",
                         "DressLep1Lep2BJet1BJet2_Pt",
                         "DressLep1Lep2BJet1BJet2_M",
                         "DressLep1Lep2Jet1_Pt",
                         "DressLep1Lep2Jet1_M",
                         "DressLep1Lep2Jet1_E",
                         "DressLep1Jet1_Pt",
                         "DressLep1Jet1_M",
                         "DressLep2Jet1_Pt",
                         "DressLep2Jet1_M",
                         "DressJet1_E",
                         "DressJet1_Pt",
                         "DressJet2_Pt",
                         "DressJetLoose1_Pt",
                         "DressJetB1_Pt",
                         "DressJetB2_Pt",
                         "DressJetB1_M",
                         "DressJetB2_M",

                         "DressLep1Jet1_DR",
                         "DressLep12Jet12_DR",
                         "DressLep12Jet12MET_DR",
                         "DressLep1Lep2Jet1_C",
                         "DressHTtot",

                         "Dressminimax",
                         "DressLep1Lep2Jet1Jet2MET_M",
                         "DressLep1Jet2_M",
                         "DressLep2Jet2_M",

                         "GenWWbb_M",
                         "GenWWbb_Pt",]

        return


    # New interface (nanoAOD-tools)
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        declareOutput(self, wrappedOutputTree, self.branches)


    def analyze(self, event):
        writeOutput(self, self.run(event, NanoAODCollection))
        return True


    # Common processing
    def run(self, event, Collection):
        allret = {}

        genparts    = [l for l in Collection(event, "GenPart")]

        # ============================ Variables not susceptible to JEC
        all_leps = [l for l in Collection(event, "GenDressedLepton")]

        leps     = [all_leps[getattr(event, 'iDressSelLep')[l]]
                    for l in range(min([getattr(event, 'nDressSelLep'), 5]))]
        leps_4m  = [l.p4() for l in leps]

        met_4m = r.TLorentzVector()
        met_4m.SetPtEtaPhiM(event.MET_fiducialGenPt, 0, event.MET_fiducialGenPhi, 0)


        all_jets = [j for j in Collection(event, "GenJet")]

        jets     = [all_jets[getattr(event, 'iDressSelJet')[j]]
                    for j in range(min([getattr(event, 'nDressSelJet'), 5]))]
        jets_4m  = [j.p4() for j in jets]

        loosejets    = [all_jets[getattr(event, 'iDressSelLooseJet')[j]]
                        for j in range(min([getattr(event, 'nDressSelLooseJet'), 5]))]
        loosejets_4m = [j.p4() for j in loosejets]

        ### b jets
        isbjet = [j.hadronFlavour == 5 for j in jets]
        jetsb_4m = []
        for i,ijet_4m in enumerate(jets_4m):
            if isbjet[i]:
                jetsb_4m.append(ijet_4m)


        allret["Dresschannel"]            = ch.NoChan
        allret["DressEXTchannel"]         = ch.NoChan
        allret["Origchannel"]             = ch.NoChan
        allret["DressisSS"]               = -99
        allret["DressminMllAFAS"]         = -99
        allret["DressLep1_Pt"]            = -99
        allret["DressLep2_Pt"]            = -99
        allret["DressLep1Lep2_Pt"]        = -99
        allret["DressLep1Lep2_DPhi"]      = -99
        allret["DressMll"]                = -99
        allret["DressLep1Lep2Jet1MET_Pz"] = -99
        allret["DressLep1Lep2Jet1MET_Pt"] = -99
        allret["DressLep1Lep2Jet1MET_M"]  = -99
        allret["DressLep1Lep2Jet1MET_Mt"] = -99
        allret["DressLep1Lep2Jet1_Pt"]    = -99
        allret["DressLep1Lep2Jet1_E"]     = -99
        allret["DressLep1Jet1_Pt"]        = -99
        allret["DressLep1Jet1_M"]         = -99
        allret["DressLep2Jet1_Pt"]        = -99
        allret["DressLep2Jet1_M"]         = -99
        allret["DressJet1_Pt"]            = -99
        allret["DressJet1_E"]             = -99
        allret["DressJet2_Pt"]            = -99
        allret["DressJetLoose1_Pt"]       = -99
        allret["DressJetB1_Pt"]           = -99
        allret["DressJetB2_Pt"]           = -99
        allret["DressJetB1_M"]            = -99
        allret["DressJetB2_M"]            = -99
        allret["DressLep1Jet1_DR"]        = -99
        allret["DressLep12Jet12_DR"]      = -99
        allret["DressLep12Jet12MET_DR"]   = -99
        allret["DressLep1Lep2Jet1_C"]     = -99
        allret["DressHTtot"]              = -99
        allret["DressLep1Jet2_M"]         = -99
        allret["DressLep2Jet2_M"]         = -99
        allret["Dressminimax"]            = -99
        allret["DressLep1Lep2BJet1BJet2_Pt"] = -99
        allret["DressLep1Lep2BJet1BJet2_M"]  = -99
        allret["DressLep1Lep2Jet1Jet2MET_M"] = -99
        allret["GenWWbb_M"]               = -99

        if   event.nGenDressedLepton == 0:
            allret["Origchannel"] = 0 # No leps

        elif event.nGenDressedLepton == 1:
            if (abs(all_leps[0].pdgId) == 13 and not all_leps[0].hasTauAnc):
                allret["Origchannel"] = 1 # one muon
            elif (abs(all_leps[0].pdgId) == 11 and not all_leps[0].hasTauAnc):
                allret["Origchannel"] = 2 # one electron
            elif ((abs(all_leps[0].pdgId == 11) or abs(all_leps[0].pdgId == 13)) and all_leps[0].hasTauAnc):
                allret["Origchannel"] = 3 # one tau that decays leptonically

        elif event.nGenDressedLepton == 2:
            if ((abs(all_leps[0].pdgId) == 13 and abs(all_leps[1].pdgId) == 13) and (not all_leps[0].hasTauAnc and not all_leps[1].hasTauAnc)):
                allret["Origchannel"] = 4 # two muons
            elif ((abs(all_leps[0].pdgId) == 11 and abs(all_leps[1].pdgId) == 11) and (not all_leps[0].hasTauAnc and not all_leps[1].hasTauAnc)):
                allret["Origchannel"] = 5 # two electrons
            elif (((abs(all_leps[0].pdgId) == 13 and abs(all_leps[1].pdgId) == 11) or
                   (abs(all_leps[0].pdgId) == 11 and abs(all_leps[1].pdgId) == 13)) and
                  (not all_leps[0].hasTauAnc and not all_leps[1].hasTauAnc)):
                allret["Origchannel"] = 6 # one electron and one muon
            elif (((abs(all_leps[0].pdgId) == 13 and not all_leps[0].hasTauAnc) and ((abs(all_leps[1].pdgId) == 11 or abs(all_leps[1].pdgId) == 13) and all_leps[1].hasTauAnc)) or
                  ((abs(all_leps[1].pdgId) == 13 and not all_leps[1].hasTauAnc) and ((abs(all_leps[0].pdgId) == 11 or abs(all_leps[0].pdgId) == 13) and all_leps[0].hasTauAnc))):
                allret["Origchannel"] = 7 # one muon and a tau that decays leptonically
            elif (((abs(all_leps[0].pdgId) == 11 and not all_leps[0].hasTauAnc) and ((abs(all_leps[1].pdgId) == 11 or abs(all_leps[1].pdgId) == 13) and all_leps[1].hasTauAnc)) or
                  ((abs(all_leps[1].pdgId) == 11 and not all_leps[1].hasTauAnc) and ((abs(all_leps[0].pdgId) == 11 or abs(all_leps[0].pdgId) == 13) and all_leps[0].hasTauAnc))):
                allret["Origchannel"] = 8 # one electron and a tau that decays leptonically
            elif (all_leps[0].hasTauAnc and all_leps[1].hasTauAnc):
                allret["Origchannel"] = 9 # two taus that decay leptonically

        elif event.nGenDressedLepton >= 3:
            allret["Origchannel"] = 10    # three or more dressed leptons

        if event.nDressSelLep >= 1:
            allret["DressLep1_Pt"] = leps_4m[0].Pt()
            if event.nDressSelLep >= 2:
                allret["DressLep2_Pt"] = leps_4m[1].Pt()
                allret["DressisSS"] = int((leps[0].pdgId > 0 and leps[1].pdgId > 0) or (leps[0].pdgId < 0 and leps[1].pdgId < 0))
                if   ((abs(leps[0].pdgId) == 13 and abs(leps[1].pdgId) == 11) or
                    (abs(leps[0].pdgId) == 11 and abs(leps[1].pdgId) == 13)):
                    allret["Dresschannel"] = ch.ElMu
                elif (abs(leps[0].pdgId) == 13 and abs(leps[1].pdgId) == 13):
                    allret["Dresschannel"] = ch.Muon
                elif (abs(leps[0].pdgId) == 11 and abs(leps[1].pdgId) == 11):
                    allret["Dresschannel"] = ch.Elec
                else:
                    allret["Dresschannel"] = ch.NoChan


                if   (((abs(leps[0].pdgId) == 13 and abs(leps[1].pdgId) == 11) or
                    (abs(leps[0].pdgId) == 11 and abs(leps[1].pdgId) == 13)) and
                    (leps[0].hasTauAnc == 0) and (leps[1].hasTauAnc == 0)):
                    allret["DressEXTchannel"] = ch.ElMu
                elif ((abs(leps[0].pdgId) == 13 and abs(leps[1].pdgId) == 13) and
                    (leps[0].hasTauAnc == 0) and (leps[1].hasTauAnc == 0)):
                    allret["DressEXTchannel"] = ch.Muon
                elif ((abs(leps[0].pdgId) == 11 and abs(leps[1].pdgId) == 11) and
                    (leps[0].hasTauAnc == 0) and (leps[1].hasTauAnc == 0)):
                    allret["DressEXTchannel"] = ch.Elec
                elif (((abs(leps[0].pdgId) == 13 and abs(leps[1].pdgId) == 11) or
                    (abs(leps[0].pdgId) == 11 and abs(leps[1].pdgId) == 13)) and
                    (leps[0].hasTauAnc == 1) and (leps[1].hasTauAnc == 1)):
                    allret["DressEXTchannel"] = ch.ElMuFromTaus
                elif ((abs(leps[0].pdgId) == 13 and abs(leps[1].pdgId) == 13) and
                    (leps[0].hasTauAnc == 1) and (leps[1].hasTauAnc == 1)):
                    allret["DressEXTchannel"] = ch.MuonFromTaus
                elif ((abs(leps[0].pdgId) == 11 and abs(leps[1].pdgId) == 11) and
                    (leps[0].hasTauAnc == 1) and (leps[1].hasTauAnc == 1)):
                    allret["DressEXTchannel"] = ch.ElecFromTaus
                elif (((abs(leps[0].pdgId) == 13 and abs(leps[1].pdgId) == 11) or
                    (abs(leps[0].pdgId) == 11 and abs(leps[1].pdgId) == 13)) and
                    (int(leps[0].hasTauAnc) + int(leps[1].hasTauAnc) == 1)):
                    allret["DressEXTchannel"] = ch.ElMuMixedFromTaus
                elif ((abs(leps[0].pdgId) == 13 and abs(leps[1].pdgId) == 13) and
                    (int(leps[0].hasTauAnc) + int(leps[1].hasTauAnc) == 1)):
                    allret["DressEXTchannel"] = ch.MuonMixedFromTaus
                elif ((abs(leps[0].pdgId) == 11 and abs(leps[1].pdgId) == 11) and
                    (int(leps[0].hasTauAnc) + int(leps[1].hasTauAnc) == 1)):
                    allret["DressEXTchannel"] = ch.ElecMixedFromTaus
                else:
                    allret["DressEXTchannel"] = ch.NoChan

                allret["DressLep1Lep2_Pt"]   = (leps_4m[0] + leps_4m[1]).Pt()
                allret["DressLep1Lep2_DPhi"] = abs(deltaPhi(leps[0], leps[1]))/r.TMath.Pi()
                allret["DressMll"]           = (leps_4m[0] + leps_4m[1]).M()

                for iL in range(len(leps)):
                    for jL in range(iL + 1, len(leps)):
                        tmpM = (leps_4m[iL] + leps_4m[jL]).M()
                        if tmpM < allret["DressminMllAFAS"] or allret["DressminMllAFAS"] == -99: allret["DressminMllAFAS"] = tmpM

                if event.nDressSelLooseJet > 0:
                    allret["DressJetLoose1_Pt"] = loosejets_4m[0].Pt()


                if event.nDressSelJet > 0:
                    allret["DressLep1Lep2Jet1MET_Pz"] = (leps_4m[0] + leps_4m[1] + jets_4m[0] + met_4m).Pz()
                    allret["DressLep1Lep2Jet1MET_Pt"] = (leps_4m[0] + leps_4m[1] + jets_4m[0] + met_4m).Pt()
                    allret["DressLep1Lep2Jet1MET_M"]  = (leps_4m[0] + leps_4m[1] + jets_4m[0] + met_4m).M()
                    allret["DressLep1Lep2Jet1MET_Mt"] = (leps_4m[0] + leps_4m[1] + jets_4m[0] + met_4m).Mt()
                    allret["DressLep1Lep2Jet1_Pt"]    = (leps_4m[0] + leps_4m[1] + jets_4m[0]).Pt()
                    allret["DressLep1Lep2Jet1_M"]     = (leps_4m[0] + leps_4m[1] + jets_4m[0]).M()
                    allret["DressLep1Lep2Jet1_E"]     = (leps_4m[0] + leps_4m[1] + jets_4m[0]).E()
                    allret["DressLep1Jet1_Pt"]        = (leps_4m[0] + jets_4m[0]).Pt()
                    allret["DressLep1Jet1_M"]         = (leps_4m[0] + jets_4m[0]).M()
                    allret["DressLep2Jet1_Pt"]        = (leps_4m[1] + jets_4m[0]).Pt()
                    allret["DressLep2Jet1_M"]         = (leps_4m[1] + jets_4m[0]).M()
                    allret["DressJet1_Pt"]            = jets_4m[0].Pt()
                    allret["DressJet1_E"]             = jets_4m[0].E()
                    allret["DressLep1Jet1_DR"]        = leps_4m[0].DeltaR(jets_4m[0])
                    allret["DressLep1Lep2Jet1_C"]     = (leps_4m[0] + leps_4m[1] + jets_4m[0]).Et() / (leps_4m[0] + leps_4m[1] + jets_4m[0]).E()
                    allret["DressHTtot"]              = leps_4m[0].Pt() + leps_4m[1].Pt() + jets_4m[0].Pt() + met_4m.Pt()

                    if event.nDressSelJet > 1:
                        allret["DressLep12Jet12_DR"]    = (leps_4m[0] + leps_4m[1]).DeltaR(jets_4m[0] + jets_4m[1])
                        allret["DressLep12Jet12MET_DR"] = (leps_4m[0] + leps_4m[1]).DeltaR(jets_4m[0] + jets_4m[1] + met_4m)
                        allret["DressLep1Jet2_M"]       = (leps_4m[0] + jets_4m[1]).M()
                        allret["DressLep2Jet2_M"]       = (leps_4m[1] + jets_4m[1]).M()
                        allret["DressJet2_Pt"]          = jets_4m[1].Pt()
                        #### WARNING: this minimax variable is only equal to that of ATLAS' PRL when the signal region is njets == 2, nbjets == 2.
                        allret["Dressminimax"] = min([max([allret["DressLep1Jet1_M"],
                                                        allret["DressLep2Jet2_M"]]),
                                                      max(allret["DressLep2Jet1_M"],
                                                        allret["DressLep1Jet2_M"])])
                        allret["DressLep1Lep2Jet1Jet2MET_M"] = (leps_4m[0] + leps_4m[1] + jets_4m[0] + jets_4m[1] + met_4m).M()
                        if len(jetsb_4m) > 1:
                            allret["DressLep1Lep2BJet1BJet2_Pt"] = (leps_4m[0] + leps_4m[1] + jetsb_4m[0] + jetsb_4m[1]).Pt()
                            allret["DressLep1Lep2BJet1BJet2_M"]  = (leps_4m[0] + leps_4m[1] + jetsb_4m[0] + jetsb_4m[1]).M() 

        thesample = ""
        for i in range(event.nDatasetName):
            thesample += str(event.DatasetName_name[i])

        if "b_bbar_4l" in thesample.lower() or "ttto2l2nu" in thesample.lower():
            the4mom = None
            for p in genparts:
                if abs(genparts[p.genPartIdxMother].pdgId) == 6:  # top quark
                    if   abs(p.pdgId) == 24:  # W boson
                        tmpvec = r.TLorentzVector()
                        tmpvec.SetPtEtaPhiM(p.pt, p.eta, p.phi, self.wmass)

                        if not the4mom:
                            the4mom = tmpvec
                        else:
                            the4mom += tmpvec

                    elif abs(p.pdgId) == 5: # b quark
                        tmpvec = r.TLorentzVector()
                        tmpvec.SetPtEtaPhiM(p.pt, p.eta, p.phi, self.bmass)

                        if not the4mom:
                            the4mom = tmpvec
                        else:
                            the4mom += tmpvec
            if the4mom:
                allret["GenWWbb_M"]  = the4mom.M()
                allret["GenWWbb_Pt"] = the4mom.Pt()
            else:
                allret["GenWWbb_M"]  = 0
                allret["GenWWbb_Pt"] = 0

            for iJ in all_jets:
                if iJ.hadronFlavour == 5:
                    if allret["DressJetB1_Pt"] == -99:
                        allret["DressJetB1_Pt"] = iJ.pt
                        allret["DressJetB1_M"]  = iJ.mass
                    else:
                        allret["DressJetB2_Pt"] = iJ.pt
                        allret["DressJetB2_M"]  = iJ.mass
                        break
        return allret
