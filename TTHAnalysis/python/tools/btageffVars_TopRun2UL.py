from math import sqrt, cos
from copy import deepcopy
import struct as st
import warnings as wr
import ROOT as r
import correctionlib._core as core

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection as NanoAODCollection 
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR,deltaPhi

from CMGTools.TTHAnalysis.treeReAnalyzer import Collection as CMGCollection
from CMGTools.TTHAnalysis.tools.nanoAOD.friendVariableProducerTools import declareOutput, writeOutput
from CMGTools.TTHAnalysis.tools.nanoAOD.TopRun2UL_modules import ch, tags, emass


class btageffVars_TopRun2UL(Module):
    def __init__(self, year_ = "2016", lepCollection = "LepGood", algo_ = 'deepJet', wp_ = "M", json_ = None, SFmeasReg = "mujets"):
        self.selecsdict = {}
        self.selecsdict["2016apv"] = lambda jet: (jet.jetId > 1 and abs(jet.eta) < 2.4 and (jet.puId >= 1 if jet.pt_nom <= 50 else 1))
        self.selecsdict["2016"]    = lambda jet: (jet.jetId > 1 and abs(jet.eta) < 2.4 and (jet.puId >= 1 if jet.pt_nom <= 50 else 1))
        self.selecsdict["2017"]    = lambda jet: (jet.jetId > 1 and abs(jet.eta) < 2.4 and (jet.puId >= 4 if jet.pt_nom <= 50 else 1))
        self.selecsdict["2018"]    = lambda jet: (jet.jetId > 1 and abs(jet.eta) < 2.4 and (jet.puId >= 4 if jet.pt_nom <= 50 else 1))
        self.lc        = lepCollection
        self.isSet     = False
        self.year      = year_
        self.selection = lambda j : True
        self.deltaRcut = 0.4
        self.nominaljecscaff = "_nom"
        self.debug = False
        self.SFmeasReg = SFmeasReg
        
        if isinstance(algo_, list):
            self.algo = algo_
        else:
            self.algo = [algo_]
            
        self.algodict = {"deepJet" : "DeepFlav"}        
        self.branchbtag = [("btag" + self.algodict[alg] + "B") for alg in self.algo]
        self.wp         = wp_
        
        if json_ == None: raise RuntimeError("FATAL: no json given.")
        self.btvjson = core.CorrectionSet.from_file(json_)
        
        self.btaggingWPsEvaluators = {}
        for algo in self.algodict:
            self.btaggingWPsEvaluators[algo]    = self.btvjson[algo + "_wp_values"]
    
        self.cutdict = {}
        for algo in self.algodict:
            for _wp_ in ["L", "M", "T"]:
                self.cutdict[algo + "_{y}_{wp}".format(y = self.year, wp = _wp_)] = self.btaggingWPsEvaluators[algo].evaluate(_wp_)

        #print(self.cutdict)

        #self.cutdict["deepCSV"] = {}; self.cutdict["deepJet"] = {}
        #self.cutdict["deepCSV"][2022] = {"L" : 0.1208, "M" : 0.4168, "T": 0.7665}
        #self.cutdict["deepJet"][2022] = {"L" : 0.0490, "M" : 0.2783, "T": 0.7100}
        
        self.flavdict = {5 : "B",
                         4 : "C",
                         0 : "L"}
                
        
        return


    # New interface (nanoAOD-tools)
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.wrappedOutputTree = wrappedOutputTree
        
        for el in ["B", "C", "L"]:
            self.wrappedOutputTree.branch("nEffSFJet{f}".format(   f = el), 'I')
            self.wrappedOutputTree.branch("EffSFJet{f}_Pt".format( f = el), 'F', lenVar = "nEffSFJet{f}".format(f = el))
            self.wrappedOutputTree.branch("EffSFJet{f}_Eta".format(f = el), 'F', lenVar = "nEffSFJet{f}".format(f = el))
            for alg in self.algo:
                self.wrappedOutputTree.branch("EffSFJet{f}_{a}Istag".format(f = el, a = alg), 'F', lenVar = "nEffSFJet{f}".format(f = el))
                self.wrappedOutputTree.branch("EffSFJet{f}_{a}SF".format(   f = el, a = alg), 'F', lenVar = "nEffSFJet{f}".format(f = el))
                self.wrappedOutputTree.branch("EffSFJet{f}_{a}SFup".format( f = el, a = alg), 'F', lenVar = "nEffSFJet{f}".format(f = el))
                self.wrappedOutputTree.branch("EffSFJet{f}_{a}SFdn".format( f = el, a = alg), 'F', lenVar = "nEffSFJet{f}".format(f = el))
        return

    def areMyJetsCleanAndGood(self, thejets, theleps):
        clist = [self.selection(jet) for jet in thejets]
        if self.debug: print(("[pythonCleaningTopRun2::areMyJetsCleanAndGood] bad/good clist:", clist))

        for iL in range(len(theleps)):
            mindr = -1; best = -1;
            for iJ in range(len(thejets)):
                tmpdr = abs(theleps[iL].p4().DeltaR(thejets[iJ].p4()))
                if (mindr < 0 or tmpdr < mindr):
                    mindr = tmpdr
                    best = iJ
            if (best > -1 and mindr < self.deltaRcut):
                clist[best] = False
        if self.debug: print(("[pythonCleaningTopRun2::areMyJetsCleanAndGood] final clist:", clist))
        return clist

    def configureCleaning(self, ev):
        if self.year != None:
            self.selection = self.selecsdict[self.year]
        elif   hasattr(ev, "year"):
            self.selection = self.selecsdict[ev.year]
            self.year_     = ev.year
        return


    def analyze(self, event):
        if not self.isSet:
            if self.debug: print("[pythonCleaningTopRun2::analyze] Configuring settings")
            self.configureCleaning(event)
            self.isSet = True
        
        self.jets    = [j for j in NanoAODCollection(event, "Jet")]
        self.jets_4m = [j.p4() for j in self.jets]
        self.leps    = [l for l in NanoAODCollection(event, self.lc)]
        writeOutput(self, self.run(event, NanoAODCollection))
        return True
   
    
    def getSF(self, pt, eta, flavour, algorithm):
        #pt_cutoff  = max(20. , min(999., pt))
        #eta_cutoff = min(2.39 if flavour == 2 else 2.49, abs(eta))

        eta = abs(eta)

        if flavour != 0:
            measurementRegion = self.SFmeasReg
        else:
            measurementRegion = "incl"

        theReader  = self.btvjson[algorithm + "_" + measurementRegion]

        SF = theReader.evaluate("central", self.wp, flavour, eta, pt)
        SFup = theReader.evaluate("up", self.wp, flavour, eta, pt)
        SFdn = theReader.evaluate("down", self.wp, flavour, eta, pt)

        return [SF, SFup, SFdn]


    # Common processing
    def run(self, event, Collection):
        # ============================ Definitions and declarations
        allret  = {}

        # ============================ Initialisations
        for el in ["B", "C", "L"]:
            allret["nEffSFJet{f}".format(   f = el)] = 0
            allret["EffSFJet{f}_Pt".format( f = el)] = []
            allret["EffSFJet{f}_Eta".format(f = el)] = []
            for alg in self.algo:
                allret["EffSFJet{f}_{a}Istag".format(f = el, a = alg)] = []
                allret["EffSFJet{f}_{a}SF".format(   f = el, a = alg)] = []
                allret["EffSFJet{f}_{a}SFup".format( f = el, a = alg)] = []
                allret["EffSFJet{f}_{a}SFdn".format( f = el, a = alg)] = []

        
        # ============================ Execution
        cleanedandgoodjets = {}
        cleanedandgoodjets[""] = self.areMyJetsCleanAndGood(self.jets, self.leps)
            
        for iJ in range(len(self.jets)):
            if cleanedandgoodjets[""][iJ]:
                thept = getattr(self.jets[iJ], "pt" + self.nominaljecscaff)
                if thept <= 20: continue
                thefl = self.jets[iJ].hadronFlavour
                allret["nEffSFJet{f}".format(   f = self.flavdict[thefl])] += 1
                allret["EffSFJet{f}_Pt".format( f = self.flavdict[thefl])].append(thept)
                allret["EffSFJet{f}_Eta".format(f = self.flavdict[thefl])].append(self.jets[iJ].eta)
                
                for iA in range(len(self.algo)):
                    SF = self.getSF(thept, self.jets[iJ].eta, thefl, self.algo[iA])
                    allret["EffSFJet{f}_{a}Istag".format(f = self.flavdict[thefl], a = self.algo[iA])].append(int(getattr(self.jets[iJ], self.branchbtag[iA]) > self.cutdict[self.algo[iA] + "_{y}_{wp}".format(y=self.year, wp=self.wp)]))
                    allret["EffSFJet{f}_{a}SF".format(   f = self.flavdict[thefl], a = self.algo[iA])].append(SF[0])
                    allret["EffSFJet{f}_{a}SFup".format( f = self.flavdict[thefl], a = self.algo[iA])].append(SF[1])
                    allret["EffSFJet{f}_{a}SFdn".format( f = self.flavdict[thefl], a = self.algo[iA])].append(SF[2])
                
        return allret
