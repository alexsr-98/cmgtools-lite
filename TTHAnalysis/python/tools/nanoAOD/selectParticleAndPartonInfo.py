from math import sqrt, cos
from copy import deepcopy
import struct as st
import warnings as wr
import ROOT as r

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection as NanoAODCollection 
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR,deltaPhi
from CMGTools.TTHAnalysis.tools.collectionSkimmer import CollectionSkimmer

from CMGTools.TTHAnalysis.treeReAnalyzer import Collection as CMGCollection
from CMGTools.TTHAnalysis.tools.nanoAOD.friendVariableProducerTools import declareOutput, writeOutput
from array import array

class selectParticleAndPartonInfo(Module):
    def __init__(self, dresslepSel_         = lambda l: True,
                       dressjetSel_         = lambda j: True,
                       dressloosejetSel_    = lambda j: True,
                       dressfwdjetSel_      = lambda j: True,
                       dressfwdloosejetSel_ = lambda j: True):

        self.branches = ["Top1_pt", "Top1_eta", "Top1_phi", ("Top1_charge", "I"),
                         "Top2_pt", "Top2_eta", "Top2_phi", ("Top2_charge", "I"),
                         ("nDressBSelJet", "I"), ("nDressBSelLooseJet", "I"),
                         ("nDressPartBSelJet", "I"), ("nDressPartBSelLooseJet", "I")]

        self.vars_common       = ["pt", "eta", "phi", "mass"]
        self.vars_dressleptons = ["pdgId", "hasTauAnc"]
        self.vars_dressjets    = ["partonFlavour", "hadronFlavour"]

        self.dresslepSel         = dresslepSel_
        self.dressjetSel         = dressjetSel_
        self.dressloosejetSel    = dressloosejetSel_
        self.dressfwdjetSel      = dressfwdjetSel_
        self.dressfwdloosejetSel = dressfwdloosejetSel_

        self.listdresslep         = []
        self.listdressjet         = []
        self.listdressloosejet    = []
        self.listdressfwdjet      = []
        self.listdressfwdloosejet = []

        return

    # New interface (nanoAOD-tools)
    def beginJob(self, histFile = None, histDirName = None):
        self.colls = {}
        self.colls["DressSelLep"        ] = CollectionSkimmer("DressSelLep",         "GenDressedLepton", floats = [], saveSelectedIndices = True, maxSize = 5)
        self.colls["DressSelJet"        ] = CollectionSkimmer("DressSelJet",         "GenJet",           floats = [], saveSelectedIndices = True, maxSize = 5)
        self.colls["DressSelLooseJet"   ] = CollectionSkimmer("DressSelLooseJet",    "GenJet",           floats = [], saveSelectedIndices = True, maxSize = 5)
        self.colls["DressSelFwdJet"     ] = CollectionSkimmer("DressSelFwdJet",      "GenJet",           floats = [], saveSelectedIndices = True, maxSize = 5)
        self.colls["DressSelFwdLooseJet"] = CollectionSkimmer("DressSelFwdLooseJet", "GenJet",           floats = [], saveSelectedIndices = True, maxSize = 5)
        return


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        for col in self.colls:
            self.colls[col].initInputTree(inputTree)

        self.initReaders(inputTree)
        declareOutput(self, wrappedOutputTree, self.branches)

        for col in self.colls:
            self.colls[col].initOutputTree(wrappedOutputTree.tree(), True);
        return


    def analyze(self, event):
        from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

        #### Obtain particle-level dressed leptons and jets colls
        dressleps = [l for l in Collection(event, "GenDressedLepton")]
        dressjets = [j for j in Collection(event, "GenJet")]

        self.listdresslep         = []
        self.listdressjet         = []
        self.listdressloosejet    = []
        self.listdressfwdjet      = []
        self.listdressfwdloosejet = []

        self.vetoedjets           = []

        otherVarsDict = {}
        for b in self.branches:
            if not isinstance(b, tuple):
                otherVarsDict[b] = -99
                if "nDress" in b: otherVarsDict[b] = 0
            else:
                otherVarsDict[b[0]] = -99
                if "nDress" in b[0]: otherVarsDict[b[0]] = 0

        for i, lep in enumerate(dressleps):
            if self.dresslepSel(lep):
                self.vetoedjets.extend(self.vetoJets(event, lep))
                self.listdresslep.append(i)

        for i, jet in enumerate(dressjets):
            if   self.dressjetSel(jet):
                #if self.isClean(event, jet):
                if i not in self.vetoedjets:
                    self.listdressjet.append(i)
                    if jet.hadronFlavour == 5: otherVarsDict["nDressBSelJet"] += 1
                    if abs(jet.partonFlavour) == 5: otherVarsDict["nDressPartBSelJet"] += 1
            elif self.dressloosejetSel(jet):
                #if self.isClean(event, jet):
                if i not in self.vetoedjets:
                    self.listdressloosejet.append(i)
                    if jet.hadronFlavour == 5: otherVarsDict["nDressBSelLooseJet"] += 1
                    if abs(jet.partonFlavour) == 5: otherVarsDict["nDressPartBSelLooseJet"] += 1
            elif self.dressfwdjetSel(jet):
                #if self.isClean(event, jet):
                if i not in self.vetoedjets:
                    self.listdressfwdjet.append(i)
            elif self.dressfwdloosejetSel(jet):
                #if self.isClean(event, jet):
                if i not in self.vetoedjets:
                    self.listdressfwdloosejet.append(i)


        if self._ttreereaderversion != event._tree._ttreereaderversion:
            for col in self.colls: self.colls[col].initInputTree(event._tree)
            self.initReaders(event._tree)

        for col in self.colls: self.colls[col].clear()

        self.colls["DressSelLep"].push_back_all(self.listdresslep)
        self.colls["DressSelJet"].push_back_all(self.listdressjet)
        self.colls["DressSelLooseJet"].push_back_all(self.listdressloosejet)
        self.colls["DressSelFwdJet"].push_back_all(self.listdressfwdjet)
        self.colls["DressSelFwdLooseJet"].push_back_all(self.listdressfwdloosejet)


        #### Obtain information about the generated tops
        partobjs = [p for p in Collection(event, "GenPart")]
        candtops = []
        for i, part in enumerate(partobjs):
            #if abs(part.pdgId) == 6 and part.status == 22:
            if abs(part.pdgId) == 6 and part.status == 62:
                candtops.append(i)

        highptind = 0; lowptind = 1;

        if len(candtops) > 2 and event.isTop:
            raise RuntimeError("FATAL: more than two tops recognised in a ttbar sample.")
        elif len(candtops) == 2:
            if partobjs[candtops[0]].pt < partobjs[candtops[1]].pt:
                highptind = 1
                lowptind  = 0

        if len(candtops) > 0:
            otherVarsDict["Top1_pt"]     = partobjs[candtops[highptind]].pt
            otherVarsDict["Top1_eta"]    = partobjs[candtops[highptind]].eta
            otherVarsDict["Top1_phi"]    = partobjs[candtops[highptind]].phi
            otherVarsDict["Top1_charge"] = 1 if (partobjs[candtops[highptind]].pdgId > 0) else -1

            if len(candtops) > 1:
                otherVarsDict["Top2_pt"]     = partobjs[candtops[lowptind]].pt
                otherVarsDict["Top2_eta"]    = partobjs[candtops[lowptind]].eta
                otherVarsDict["Top2_phi"]    = partobjs[candtops[lowptind]].phi
                otherVarsDict["Top2_charge"] = 1 if (partobjs[candtops[lowptind]].pdgId > 1) else -1


        #### Write branches not in colls (and their counts).
        for b in self.branches:
            if not isinstance(b, tuple):
                self.wrappedOutputTree.fillBranch(b, otherVarsDict[b])
            else:
                self.wrappedOutputTree.fillBranch(b[0], otherVarsDict[b[0]])

        return True


    def initReaders(self, tree):
        self._ttreereaderversion = tree._ttreereaderversion
        for col in ["GenDressedLepton", "GenJet"]:
            setattr(self, 'n' + col, tree.valueReader('n' + col))
            _vars = self.vars_common[:]
            if col == 'GenDressedLepton': _vars.extend(self.vars_dressleptons)
            if col == "GenJet":           _vars.extend(self.vars_dressjets)

            for B in _vars:
                if type(B) == tuple:
                    setattr(self, "%s_%s"%(col, B[0]), tree.arrayReader("%s_%s"%(col, B[1])))
                else:
                    setattr(self, "%s_%s"%(col, B),    tree.arrayReader("%s_%s"%(col, B)))
        return True


    #### Old, standard cleaning
    def isClean(self, ev, j, dR = 0.4):
        tmpleps = [l for l in Collection(ev, "GenDressedLepton")]

        for iL in self.listdresslep:
            if abs(j.p4().DeltaR(tmpleps[iL].p4())) < dR:
                return False
        return True


    #### ttH-like cleaning
    def vetoJets(self, ev, l, dR = 0.4):
        tmpjets = [j for j in Collection(ev, "GenJet")]
        vetoed  = []

        for i, j in enumerate(tmpjets):
            if abs(j.p4().DeltaR(l.p4())) < dR:
                vetoed.append(i)
        return vetoed
