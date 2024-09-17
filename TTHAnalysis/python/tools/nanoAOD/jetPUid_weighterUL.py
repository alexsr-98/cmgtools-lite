from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop   import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel   import Collection
from CMGTools.TTHAnalysis.tools.nanoAOD.friendVariableProducerTools import writeOutput
import ROOT as r
from copy import deepcopy
import correctionlib._core as core


class jetPUid_weighterUL(Module):
    def __init__(self, fich, wp = "L", label = "", year = 2016, lepCollection = "LepGood",
                 jecvars = ["jesTotal", "jer"], lepenvars = ["mu"], debug = False):

        self.wp    = wp
        self.label = label
        self.year  = year
        self.debug = debug

        #### HERE W/O JETPUID AND ALL OTHER CUTS except pt
        self.selecsdict = {}
        self.selecsdict[2016] = lambda jet: (jet.jetId > 1 and (jet.puId >= 1 if jet.pt_nom <= 50 else 1))
        self.selecsdict[2017] = lambda jet: (jet.jetId > 1 and (jet.puId >= 4 if jet.pt_nom <= 50 else 1))
        self.selecsdict[2018] = lambda jet: (jet.jetId > 1 and (jet.puId >= 4 if jet.pt_nom <= 50 else 1))
        self.lc     = lepCollection
        self.isSet     = False
        self.selection = lambda j : True
        self.deltaRcut = 0.4

        self.systsJEC   = {0: ""}
        self.systsLepEn = {}
        self.nominaljecscaff = "_nom"

        if   len(jecvars):
            for i, var in enumerate(jecvars):
                self.systsJEC[i+1]    = "_%sUp"%var
                self.systsJEC[-(i+1)] = "_%sDown"%var
        if   len(lepenvars):
            for i, var in enumerate(lepenvars):
                self.systsLepEn[i+1]    = "_%sUp"%var
                self.systsLepEn[-(i+1)] = "_%sDown"%var

        self.ret = {}

        self.jmarjson = core.CorrectionSet.from_file(fich)
        map_name = "PUJetID_eff"
        self.evalPUJetID = self.jmarjson[map_name]
        return


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.wrappedOutputTree = wrappedOutputTree

        self.wrappedOutputTree.branch("jetPUidWeight" + self.label + "_tag_Up"   , "F")
        self.wrappedOutputTree.branch("jetPUidWeight" + self.label + "_tag_Dn"   , "F")

        for delta,jecVar in self.systsJEC.items():
            self.wrappedOutputTree.branch("jetPUidWeight" + self.label + jecVar , "F")

        for delta,lepVar in self.systsLepEn.items():
            self.wrappedOutputTree.branch("jetPUidWeight" + self.label + lepVar , "F")

        return


    def analyze(self, event):
        if not self.isSet:
            if self.debug: print("[pythonCleaningTopRun2::analyze] Configuring settings")
            self.configureCleaning(event)
            self.isSet = True
        self.event = event
        self.ret.clear()
        self.all_jets  = [j for j in Collection(event, "Jet")]
        self.leps = [l for l in Collection(event, self.lc)]
        self.computeWeights(event)

        writeOutput(self, self.ret)
        return True


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
        return

    def computeWeights(self, event):
        cleanedandgoodjets = {}
        cleanedandgoodjets[""] = self.areMyJetsCleanAndGood(self.all_jets, self.leps)

        for delta,lvar in self.systsLepEn.items():
            tmpleps = [l for l in Collection(event, self.lc + lvar[1:])]
            cleanedandgoodjets[lvar] = self.areMyJetsCleanAndGood(self.all_jets, tmpleps)

        for delta,jecVar in self.systsJEC.items():
            sftag       = 1. 
            unctagup    = 1.
            unctagdn    = 1.

            jetjecsysscaff = (jecVar if jecVar != "" else self.nominaljecscaff)

            for iJ in range(len(self.all_jets)):
                if cleanedandgoodjets[""][iJ]:
                    thept   = getattr(self.all_jets[iJ], "pt" + jetjecsysscaff)
                    if thept <= 20 or thept >= 50: continue
                    isPUjet = getattr(self.all_jets[iJ], "genJetIdx") == -1
                    if isPUjet: continue


                    sftag *= self.evalPUJetID.evaluate(self.all_jets[iJ].eta, thept, "nom", self.wp)

                    if jecVar == "":              
                        unctagup    *= self.evalPUJetID.evaluate(self.all_jets[iJ].eta, thept, "up", self.wp)
                        unctagdn    *= self.evalPUJetID.evaluate(self.all_jets[iJ].eta, thept, "down", self.wp)

            self.ret["jetPUidWeight" + self.label + jecVar] = sftag

            if jecVar == "":
                self.ret["jetPUidWeight" + self.label + "_tag_Up"]    = unctagup
                self.ret["jetPUidWeight" + self.label + "_tag_Dn"]    = unctagdn

                for ldelta,lepVar in self.systsLepEn.items():
                    sftaglep       = 1. 
                    jetjecsysscaff = self.nominaljecscaff

                    for iJ in range(len(self.all_jets)):
                        if cleanedandgoodjets[lepVar][iJ]:
                            thept   = getattr(self.all_jets[iJ],  "pt" + jetjecsysscaff)
                            if thept <= 20 or thept >= 50: continue
                            isPUjet = getattr(self.all_jets[iJ], "genJetIdx") == -1
                            if isPUjet: continue

                            sftaglep *= self.evalPUJetID.evaluate(self.all_jets[iJ].eta, thept, "nom", self.wp)

                    self.ret["jetPUidWeight" + self.label + lepVar] = sftaglep

        return




