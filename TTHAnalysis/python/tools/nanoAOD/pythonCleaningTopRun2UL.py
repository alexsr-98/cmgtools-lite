import ROOT as r
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from CMGTools.TTHAnalysis.tools.collectionSkimmer import CollectionSkimmer
from CMGTools.TTHAnalysis.tools.nanoAOD.friendVariableProducerTools import declareOutput
import warnings as wr
import os, sys

class pythonCleaningTopRun2UL(Module):
    def __init__(self, label, jetPts, jetPtNoisyFwd, jetCollection = 'Jet', lepCollection = "LepGood",
                 btagWP_ = "M", isMC = True, deltaRcut = 0.4,
                 jecvars = ["jesTotalCorr", "jesTotalUnCorr", "jer"],
                 lepenvars = ["mu"], algo = "DeepFlav",
                 year_ = None, debug = False):

        self.label = "" if (label in ["", None]) else ("_" + label)

        #### NOTE: here the tight selections for jet (excluding eta and pt requirements) are hardcoded per year
        self.selecsdict = {}

        # Con JetPUID (revisar con cuidao https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDUL)
        """
        if isMC:
            self.selecsdict[2016] = lambda jet: (jet.jetId > 1 and (jet.puId >= 3 if jet.pt_nom <= 50 else 1))
            self.selecsdict[2017] = lambda jet: (jet.jetId > 1 and (jet.puId >= 6 if jet.pt_nom <= 50 else 1))
            self.selecsdict[2018] = lambda jet: (jet.jetId > 1 and (jet.puId >= 6 if jet.pt_nom <= 50 else 1))
        else:
            self.selecsdict[2016] = lambda jet: (jet.jetId > 1 and (jet.puId >= 3 if jet.pt <= 50 else 1))
            self.selecsdict[2017] = lambda jet: (jet.jetId > 1 and (jet.puId >= 6 if jet.pt <= 50 else 1))
            self.selecsdict[2018] = lambda jet: (jet.jetId > 1 and (jet.puId >= 6 if jet.pt <= 50 else 1))
        """
        
        # Estandar
        self.selecsdict[2016] = lambda jet: jet.jetId > 1
        self.selecsdict[2017] = lambda jet: jet.jetId > 1
        self.selecsdict[2018] = lambda jet: jet.jetId > 1
        

        self.jc     = jetCollection
        self.lc     = lepCollection

        self.jetPts = jetPts
        self.jetPtNoisyFwd = jetPtNoisyFwd

        self.isMC      = isMC
        self.nomjecscaff = "_nom"
        self.jecvars   = jecvars
        self.lepenvars = lepenvars

        if not self.isMC:
            self.nomjecscaff = ""
            self.jecvars     = []
            if "elscale" in self.lepenvars:
                self.lepenvars = ["elscale"]
            else:
                self.lepenvars = []

        self.year      = year_
        self.isSet     = False
        self.selection = lambda j : True

        self.algodict = {"DeepFlav" : "DeepFlav",
                         "DeepCSV"  : "Deep"}

        self.algoscaff = algo + "_{y}_{wp}"
        self.jetBTag   = "btag" + self.algodict[algo] + "B"
        self.btagWP    = btagWP_
        ### WP extracted on 2022-04-27 from:
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL16preVFP  # APV
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL16postVFP
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL17
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL18
        self.btagWPs   = {"DeepCSV_2016apv_L" : 0.2027,
                          "DeepCSV_2016apv_M" : 0.6001,
                          "DeepCSV_2016apv_T" : 0.8819,
                          "DeepCSV_2016_L"    : 0.1918,
                          "DeepCSV_2016_M"    : 0.5847,
                          "DeepCSV_2016_T"    : 0.8767,
                          "DeepCSV_2017_L"    : 0.1355,
                          "DeepCSV_2017_M"    : 0.4506,
                          "DeepCSV_2017_T"    : 0.7738,
                          "DeepCSV_2018_L"    : 0.1208,
                          "DeepCSV_2018_M"    : 0.4168,
                          "DeepCSV_2018_T"    : 0.7665,
                          "DeepFlav_2016apv_L": 0.0508,
                          "DeepFlav_2016apv_M": 0.2598,
                          "DeepFlav_2016apv_T": 0.6502,
                          "DeepFlav_2016_L"   : 0.0480,
                          "DeepFlav_2016_M"   : 0.2489,
                          "DeepFlav_2016_T"   : 0.6377,
                          "DeepFlav_2017_L"   : 0.0532,
                          "DeepFlav_2017_M"   : 0.3040,
                          "DeepFlav_2017_T"   : 0.7476,
                          "DeepFlav_2018_L"   : 0.0490,
                          "DeepFlav_2018_M"   : 0.2783,
                          "DeepFlav_2018_T"   : 0.7100,
                         }
        self.btagWPcut = 0.
        self.deltaRcut = deltaRcut
        self.debug     = debug
        return


    #### Running methods
    def beginJob(self, histFile = None, histDirName = None):
        # Configuring variations (lepton energy corrs. and JEC).
        self.systsLepEn = {}
        self.systsJEC   = {0 : ""}
        for sys in range(len(self.jecvars)):
            self.systsJEC[sys+1]    = '_' + self.jecvars[sys] + 'Up'
            self.systsJEC[-(sys+1)] = '_' + self.jecvars[sys] + 'Down'
        for sys in range(len(self.lepenvars)):
            self.systsLepEn[sys+1]    = '_' + self.lepenvars[sys] + 'Up'
            self.systsLepEn[-(sys+1)] = '_' + self.lepenvars[sys] + 'Down'


        # Configuring output branches
        self._outjetvars    = [x%self.jc for x in ['nB%sSelMedium%%d','n%sSel%%d']]
        self._outfwdjetvars = [x%self.jc for x in ['nFwd%sSel%%d']]
        self.outjetvars  = []
        for jetPt in self.jetPts:
            self.outjetvars.extend([(x%jetPt+y, 'I') for x in self._outjetvars    for y in list(self.systsJEC.values())   if "nB" in x])
            self.outjetvars.extend([(x%jetPt+y, 'I') for x in self._outfwdjetvars for y in list(self.systsJEC.values())   if "nB" in x])

            self.outjetvars.extend([(x%jetPt+y, 'I') for x in self._outjetvars    for y in list(self.systsLepEn.values()) if "nB" in x])
            self.outjetvars.extend([(x%jetPt+y, 'I') for x in self._outfwdjetvars for y in list(self.systsLepEn.values()) if "nB" in x])

        self.branches = [(var + self.label, _type) for var,_type in self.outjetvars]

        # Configuring C++ helpers for JEC variations as well as nominal results
        self.colls         = {}
        for delta,sys in self.systsJEC.items():
            self.colls["jetsSup"    + sys] = CollectionSkimmer("{col}Sel{pt}".format(   col = self.jc, pt = self.jetPts[0]) + sys + self.label, self.jc, floats = [], saveSelectedIndices = True)
            self.colls["jetsInf"    + sys] = CollectionSkimmer("{col}Sel{pt}".format(   col = self.jc, pt = self.jetPts[1]) + sys + self.label, self.jc, floats = [], saveSelectedIndices = True)
            self.colls["fwdjetsSup" + sys] = CollectionSkimmer("Fwd{col}Sel{pt}".format(col = self.jc, pt = self.jetPts[0]) + sys + self.label, self.jc, floats = [], saveSelectedIndices = True)
            self.colls["fwdjetsInf" + sys] = CollectionSkimmer("Fwd{col}Sel{pt}".format(col = self.jc, pt = self.jetPts[1]) + sys + self.label, self.jc, floats = [], saveSelectedIndices = True)

        for delta,sys in self.systsLepEn.items():
            self.colls["jetsSup"    + sys] = CollectionSkimmer("{col}Sel{pt}".format(   col = self.jc, pt = self.jetPts[0]) + sys + self.label, self.jc, floats = [], saveSelectedIndices = True)
            self.colls["jetsInf"    + sys] = CollectionSkimmer("{col}Sel{pt}".format(   col = self.jc, pt = self.jetPts[1]) + sys + self.label, self.jc, floats = [], saveSelectedIndices = True)
            self.colls["fwdjetsSup" + sys] = CollectionSkimmer("Fwd{col}Sel{pt}".format(col = self.jc, pt = self.jetPts[0]) + sys + self.label, self.jc, floats = [], saveSelectedIndices = True)
            self.colls["fwdjetsInf" + sys] = CollectionSkimmer("Fwd{col}Sel{pt}".format(col = self.jc, pt = self.jetPts[1]) + sys + self.label, self.jc, floats = [], saveSelectedIndices = True)
        return



    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        for tag,helper in self.colls.items():
            helper.initInputTree(inputTree)

        declareOutput(self, wrappedOutputTree, self.branches)

        for tag,helper in self.colls.items():
            helper.initOutputTree(wrappedOutputTree.tree(), True)
        return



    def analyze(self, event):
        if self.debug: print("[pythonCleaningTopRun2::analyze]")
        # 1) Configure settings for the first event processed
        if not self.isSet:
            if self.debug: print("[pythonCleaningTopRun2::analyze] Configuring settings")
            self.configureCleaning(event)
            self.isSet = True


        # 2) Process
        if self.debug: print("[pythonCleaningTopRun2::analyze] Beginning processing")

        leps = [l for l in Collection(event, self.lc)]
        jets = [j for j in Collection(event, self.jc)]

        cleanedandgoodjets = {}
        cleanedandgoodjets[""] = self.areMyJetsCleanAndGood(jets, leps)

        finalJetDict     = {}
        finalnBtagsDict0 = {}
        finalnBtagsDict1 = {}

        for key in self.colls:
            finalJetDict[key] = []
            self.colls[key].initEventOpt(event) #### NOTE: using optimised set of methods. DO NOT MIX THEM WITH LEGACY ONES!

        for delta,jvar in self.systsJEC.items():
            finalnBtagsDict0[jvar] = 0
            finalnBtagsDict1[jvar] = 0

        for delta,lvar in self.systsLepEn.items():
            finalnBtagsDict0[lvar] = 0
            finalnBtagsDict1[lvar] = 0
            leps = [l for l in Collection(event, self.lc + lvar[1:])]
            cleanedandgoodjets[lvar] = self.areMyJetsCleanAndGood(jets, leps)


        if self.debug: print("[pythonCleaningTopRun2::analyze] Beginning jet loop")
        for iJ in range(len(jets)):
            if cleanedandgoodjets[""][iJ]:
                abseta = abs(jets[iJ].eta)
                if   abseta < 2.4:
                    tmpisbtag = getattr(jets[iJ], self.jetBTag) > self.btagWPcut
                    for delta,jvar in self.systsJEC.items():
                        if   getattr(jets[iJ], "pt" + (self.nomjecscaff if jvar == "" else jvar)) > self.jetPts[0]:
                            finalJetDict["jetsSup"    + jvar].append(iJ)
                            if tmpisbtag: finalnBtagsDict0[jvar] += 1
                        elif getattr(jets[iJ], "pt" + (self.nomjecscaff if jvar == "" else jvar)) > self.jetPts[1]:
                            finalJetDict["jetsInf"    + jvar].append(iJ)
                            if tmpisbtag: finalnBtagsDict1[jvar] += 1

                elif abseta < 2.7:
                    for delta,jvar in self.systsJEC.items():
                        if   getattr(jets[iJ], "pt" + (self.nomjecscaff if jvar == "" else jvar)) > self.jetPts[0]:
                            finalJetDict["fwdjetsSup" + jvar].append(iJ)
                        elif getattr(jets[iJ], "pt" + (self.nomjecscaff if jvar == "" else jvar)) > self.jetPts[1]:
                            finalJetDict["fwdjetsInf" + jvar].append(iJ)

                elif abseta < 3.0:
                    for delta,jvar in self.systsJEC.items():
                        if   getattr(jets[iJ], "pt" + (self.nomjecscaff if jvar == "" else jvar)) > self.jetPtNoisyFwd:
                            finalJetDict["fwdjetsSup" + jvar].append(iJ)

                elif abseta < 5.0:
                    for delta,jvar in self.systsJEC.items():
                        if   getattr(jets[iJ], "pt" + (self.nomjecscaff if jvar == "" else jvar)) > self.jetPts[0]:
                            finalJetDict["fwdjetsSup" + jvar].append(iJ)
                        elif getattr(jets[iJ], "pt" + (self.nomjecscaff if jvar == "" else jvar)) > self.jetPts[1]:
                            finalJetDict["fwdjetsInf" + jvar].append(iJ)

            for delta,lvar in self.systsLepEn.items():
                if cleanedandgoodjets[lvar][iJ]:
                    abseta = abs(jets[iJ].eta)
                    if   abseta < 2.4:
                        tmpisbtag = getattr(jets[iJ], self.jetBTag) > self.btagWPcut
                        if   getattr(jets[iJ], "pt" + self.nomjecscaff) > self.jetPts[0]:
                            finalJetDict["jetsSup"    + lvar].append(iJ)
                            if tmpisbtag: finalnBtagsDict0[lvar] += 1
                        elif getattr(jets[iJ], "pt" + self.nomjecscaff) > self.jetPts[1]:
                            finalJetDict["jetsInf"    + lvar].append(iJ)
                            if tmpisbtag: finalnBtagsDict1[lvar] += 1

                    elif abseta < 2.7:
                        if   getattr(jets[iJ], "pt" + self.nomjecscaff) > self.jetPts[0]:
                            finalJetDict["fwdjetsSup" + lvar].append(iJ)
                        elif getattr(jets[iJ], "pt" + self.nomjecscaff) > self.jetPts[1]:
                            finalJetDict["fwdjetsInf" + lvar].append(iJ)

                    elif abseta < 3.0:
                        if   getattr(jets[iJ], "pt" + self.nomjecscaff) > self.jetPtNoisyFwd:
                            finalJetDict["fwdjetsSup" + lvar].append(iJ)

                    elif abseta < 5.0:
                        if   getattr(jets[iJ], "pt" + self.nomjecscaff) > self.jetPts[0]:
                            finalJetDict["fwdjetsSup" + lvar].append(iJ)
                        elif getattr(jets[iJ], "pt" + self.nomjecscaff) > self.jetPts[1]:
                            finalJetDict["fwdjetsInf" + lvar].append(iJ)


        if self.debug: print("[pythonCleaningTopRun2::analyze] Finished jet loop")

        # 3) Save results
        # First, in the collectionskimmers
        if self.debug: print("[pythonCleaningTopRun2::analyze] Saving results in collectionskimmers and non-collectionskimmer variables")

        for delta,sys in self.systsJEC.items():
            self.colls["jetsSup"    + sys].push_back_allOpt(finalJetDict["jetsSup"    + sys])
            self.colls["jetsInf"    + sys].push_back_allOpt(finalJetDict["jetsInf"    + sys])
            self.colls["fwdjetsSup" + sys].push_back_allOpt(finalJetDict["fwdjetsSup" + sys])
            self.colls["fwdjetsInf" + sys].push_back_allOpt(finalJetDict["fwdjetsInf" + sys])

            # Afterwards, those variables not in CS or their counts
            self.wrappedOutputTree.fillBranch('nBJetSelMedium' + str(self.jetPts[0]) + sys + self.label, finalnBtagsDict0[sys])
            self.wrappedOutputTree.fillBranch('nBJetSelMedium' + str(self.jetPts[1]) + sys + self.label, finalnBtagsDict1[sys])
            
            self.colls["jetsSup"    + sys].finishEventOpt()
            self.colls["jetsInf"    + sys].finishEventOpt()
            self.colls["fwdjetsSup" + sys].finishEventOpt()
            self.colls["fwdjetsInf" + sys].finishEventOpt()

        for delta,sys in self.systsLepEn.items():
            self.colls["jetsSup"    + sys].push_back_allOpt(finalJetDict["jetsSup"    + sys])
            self.colls["jetsInf"    + sys].push_back_allOpt(finalJetDict["jetsInf"    + sys])
            self.colls["fwdjetsSup" + sys].push_back_allOpt(finalJetDict["fwdjetsSup" + sys])
            self.colls["fwdjetsInf" + sys].push_back_allOpt(finalJetDict["fwdjetsInf" + sys])

            # Afterwards, those variables not in CS or their counts
            self.wrappedOutputTree.fillBranch('nBJetSelMedium' + str(self.jetPts[0]) + sys + self.label, finalnBtagsDict0[sys])
            self.wrappedOutputTree.fillBranch('nBJetSelMedium' + str(self.jetPts[1]) + sys + self.label, finalnBtagsDict1[sys])

            self.colls["jetsSup"    + sys].finishEventOpt()
            self.colls["jetsInf"    + sys].finishEventOpt()
            self.colls["fwdjetsSup" + sys].finishEventOpt()
            self.colls["fwdjetsInf" + sys].finishEventOpt()

        if self.debug: print("[pythonCleaningTopRun2::analyze] Event finished.")

        return True



    #### Other methods
    def configureCleaning(self, ev):
        if self.year != None:
            self.selection = self.selecsdict[int(self.year) if not "apv" in self.year else 2016]
            self.btagWPcut = self.btagWPs[self.algoscaff.format(y = self.year, wp = self.btagWP)]
        else:
            wr.warn("WARNING: as you did not choose year, we will use deepcsv as algorithm with random values for the working points.")
            self.btagWPcut = self.btagWPs["DeepCSVM"]
        return


#### WEIRD CLEANING (ttH-like)
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


#### NORMAL CLEANING (top-like)
    #def areMyJetsCleanAndGood(self, thejets, theleps):
        #clist = []
        #for iJ in range(len(thejets)):
            #clist.append(self.isCleanAndGood(thejets[iJ], theleps))
        #return clist


    def isCleanAndGood(self, jet, leps):
        if not self.selection(jet): return False
        for iL in range(len(leps)):
            if abs(jet.p4().DeltaR(leps[iL].p4())) < self.deltaRcut:
                return False
        return True
