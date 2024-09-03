from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop   import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel   import Collection
from CMGTools.TTHAnalysis.tools.nanoAOD.friendVariableProducerTools import writeOutput
import ROOT as r
from copy import deepcopy


class jetPUid_weighterUL(Module):
    def __init__(self, fich, wp = "T", label = "", year = 2017, lepCollection = "LepGood",
                 jecvars = ["jesTotal", "jer"], lepenvars = ["mu"], debug = False):

        self.wp    = wp
        self.label = label
        self.year  = year
        self.debug = debug

        #### HERE W/O JETPUID AND ALL OTHER CUTS except pt
        self.selecsdict = {}
        self.selecsdict[2016] = lambda jet: (jet.jetId > 1 and (jet.puId >= 3 if jet.pt_nom <= 50 else 1))
        self.selecsdict[2017] = lambda jet: (jet.jetId > 1 and (jet.puId >= 6 if jet.pt_nom <= 50 else 1))
        self.selecsdict[2018] = lambda jet: (jet.jetId > 1 and (jet.puId >= 6 if jet.pt_nom <= 50 else 1))
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

        #### Get effs.
        thef = r.TFile.Open(fich, "read")
        self.h_eff_tag     = deepcopy(thef.Get("hPUJetID_{w}_MCEff".format(w = self.wp)).Clone())
        #self.h_eff_mistag = deepcopy(thef.Get("h2_mistag_mc{y}_{w}".format(w = self.wp)).Clone())

        #### Get SF
        self.h_sf_tag      = deepcopy(thef.Get("hPUJetID_{w}_nom".format(w = self.wp)).Clone())
        #self.h_sf_mistag = deepcopy(thef.Get("h2_mistag_sf{y}_{w}".format(w = self.wp)).Clone())

        # and their systematic uncertainties...
        self.h_sf_tag_systup = deepcopy(thef.Get("hPUJetID_{w}_up".format(w = self.wp)).Clone())
        self.h_sf_tag_systdn = deepcopy(thef.Get("hPUJetID_{w}_down".format(w = self.wp)).Clone())
        #self.h_sf_mistag_syst = deepcopy(thef.Get("h2_mistag_sf{y}_{w}_Systuncty".format(w = self.wp)).Clone())
        thef.Close()
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
        self.jets = []
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
        if   hasattr(ev, "year"):
            self.selection = self.selecsdict[ev.year]
        elif self.year != None:
            self.selection = self.selecsdict[self.year]
        return

    def computeWeights(self, event):
        cleanedandgoodjets = {}
        cleanedandgoodjets[""] = self.areMyJetsCleanAndGood(self.all_jets, self.leps)

        for delta,lvar in self.systsLepEn.items():
            tmpleps = [l for l in Collection(event, self.lc + lvar[1:])]
            cleanedandgoodjets[lvar] = self.areMyJetsCleanAndGood(self.all_jets, tmpleps)

        for delta,jecVar in self.systsJEC.items():
            mcTag       = 1.
            mcNoTag     = 1.
            dataTag     = 1.
            dataNoTag   = 1.
            unctagup    = 1.
            unctagdn    = 1.
            uncmistagup = 1.
            uncmistagdn = 1.

            jetjecsysscaff = (jecVar if jecVar != "" else self.nominaljecscaff)

            for iJ in range(len(self.all_jets)):
                if cleanedandgoodjets[""][iJ]:
                    thept   = getattr(self.all_jets[iJ], "pt" + jetjecsysscaff)
                    if thept <= 20 or thept >= 50: continue
                    isPUjet = getattr(self.all_jets[iJ], "genJetIdx") == -1

                    eff = self.getEffFromHisto(thept, self.all_jets[iJ].eta, isPUjet) # WARNING: do not use later without changing it back!
                    SF  = self.getSF (thept, self.all_jets[iJ].eta, isPUjet)

                    istag = self.all_jets[iJ].puId == 7

                    if self.debug and jecVar == "":
                        print("[jetPUid_weighter::computeWeights] -", thept, self.all_jets[iJ].eta, eff, SF)

                    if eff == 0:
                        print("[jetPUid_weighter::computeWeights] - WARNING: event efficiency is zero. To prevent a ZeroDivision error, we will fix its value to 10^-5.")
                        eff = 1e-5

                    if istag:
                        mcTag   *= eff
                        dataTag *= eff*SF[0]

                        if (isPUjet and jecVar == ""): #### PU jet
                            unctagup    *= eff*SF[1]
                            unctagdn    *= eff*SF[2]
                            uncmistagup *= eff*SF[0]
                            uncmistagdn *= eff*SF[0]

                        elif jecVar == "":              #### non-PU jet
                            unctagup    *= eff*SF[0]
                            unctagdn    *= eff*SF[0]
                            uncmistagup *= eff*SF[1]
                            uncmistagdn *= eff*SF[2]

                    else:
                        mcNoTag   *= (1 - eff      )
                        dataNoTag *= (1 - eff*SF[0])

                        if (isPUjet and jecVar == ""):
                            unctagup    *= (1 - eff*SF[1])
                            unctagdn    *= (1 - eff*SF[2])
                            uncmistagup *= (1 - eff*SF[0])
                            uncmistagdn *= (1 - eff*SF[0])

                        elif jecVar == "":
                            unctagup    *= (1 - eff*SF[0])
                            unctagdn    *= (1 - eff*SF[0])
                            uncmistagup *= (1 - eff*SF[1])
                            uncmistagdn *= (1 - eff*SF[2])


            #central = (dataNoTag * dataTag ) / ( mcNoTag * mcTag )
            central = ( dataTag ) / ( mcTag )

            self.ret["jetPUidWeight" + self.label + jecVar] = central


            # b-tagging & mistagging eff. unc. varied weights
            if jecVar == "":
                #self.ret["jetPUidWeight" + self.label + "_tag_Up"]    = unctagup    / ( mcNoTag * mcTag )
                self.ret["jetPUidWeight" + self.label + "_tag_Up"]    = unctagup    / ( mcTag )
                #self.ret["jetPUidWeight" + self.label + "_tag_Dn"]    = unctagdn    / ( mcNoTag * mcTag )
                self.ret["jetPUidWeight" + self.label + "_tag_Dn"]    = unctagdn    / ( mcTag )
                #self.ret["jetPUidWeight" + self.label + "_mistag_Up"] = uncmistagup / ( mcNoTag * mcTag )
                #self.ret["jetPUidWeight" + self.label + "_mistag_Dn"] = uncmistagdn / ( mcNoTag * mcTag )

                for ldelta,lepVar in self.systsLepEn.items():
                    mcTag     = 1.
                    mcNoTag   = 1.
                    dataTag   = 1.
                    dataNoTag = 1.

                    self.jets    = [self.all_jets[getattr(event, 'iJetSel30{v}_Recl'.format(v = lepVar))[j]]
                                for j in range(min([getattr(event, 'nJetSel30{v}_Recl'.format(v = lepVar)), 5]))]

                    jetjecsysscaff = self.nominaljecscaff

                    for iJ in range(len(self.all_jets)):
                        if cleanedandgoodjets[lepVar][iJ]:
                            thept   = getattr(self.all_jets[iJ],  "pt" + jetjecsysscaff)
                            if thept <= 20 or thept >= 50: continue
                            isPUjet = getattr(self.all_jets[iJ], "genJetIdx") == -1

                            eff = self.getEffFromHisto(thept, self.all_jets[iJ].eta, isPUjet) # WARNING: do not use later without changing it back!
                            SF  = self.getSF (thept, self.all_jets[iJ].eta, isPUjet)

                            istag = self.all_jets[iJ].puId == 7

                            if eff == 0:
                                print("[jetPUid_weighter::computeWeights] - WARNING: jet eff is zero. To prevent a ZeroDivision error, we will fix its value to 10^-5.")
                                eff = 1e-5

                            if istag:
                                mcTag   *= eff
                                dataTag *= eff*SF[0]
                            else:
                                mcNoTag   *= (1 - eff      )
                                dataNoTag *= (1 - eff*SF[0])

                    #central = (dataNoTag * dataTag ) / ( mcNoTag * mcTag )
                    central = (dataTag ) / ( mcTag )

                    self.ret["jetPUidWeight" + self.label + lepVar] = central

        return


    def getEffFromHisto(self, pt, eta, ispujet):
        if ispujet:
            histo = self.h_eff_tag
        else:
            #histo = self.h_eff_mistag
            return 1.


        normal = True # normal := {x: pt, y: eta}
        if (histo.GetXaxis().GetBinUpEdge(histo.GetNbinsX())) < 10: normal = False

        tmpeta = eta

        if    ((normal     and histo.GetYaxis().GetBinLowEdge(1) >= 0)
            or (not normal and histo.GetXaxis().GetBinLowEdge(1) >= 0)):
            tmpeta = abs(tmpeta)

        if normal:
            xbin = max(1, min(histo.GetNbinsX(), histo.GetXaxis().FindBin(pt)))
            ybin = max(1, min(histo.GetNbinsY(), histo.GetYaxis().FindBin(tmpeta)))
        else:
            xbin = max(1, min(histo.GetNbinsX(), histo.GetXaxis().FindBin(tmpeta)))
            ybin = max(1, min(histo.GetNbinsY(), histo.GetYaxis().FindBin(pt)))

        #print("pt:", pt, "eta:", eta, xbin, ybin)
        return histo.GetBinContent(xbin, ybin)


    def getSF(self, pt, eta, ispujet):
        if ispujet:
            histo   = self.h_sf_tag
            histoup = self.h_sf_tag_systup
            histodn = self.h_sf_tag_systdn
        else:
            return [1., 1., 1.]
            #histo   = self.h_sf_mistag
            #histoup = self.h_sf_tag_systup
            #histodn = self.h_sf_tag_systdn

        normal = True # normal := {x: pt, y: eta}
        if (histo.GetXaxis().GetBinUpEdge(histo.GetNbinsX())) < 10: normal = False

        tmpeta = eta

        if    ((normal     and histo.GetYaxis().GetBinLowEdge(1) >= 0)
            or (not normal and histo.GetXaxis().GetBinLowEdge(1) >= 0)):
            tmpeta = abs(tmpeta)

        if normal:
            xbin = max(1, min(histo.GetNbinsX(), histo.GetXaxis().FindBin(pt)))
            ybin = max(1, min(histo.GetNbinsY(), histo.GetYaxis().FindBin(tmpeta)))
        else:
            xbin = max(1, min(histo.GetNbinsX(), histo.GetXaxis().FindBin(tmpeta)))
            ybin = max(1, min(histo.GetNbinsY(), histo.GetYaxis().FindBin(pt)))

        SF   = histo.GetBinContent(xbin, ybin)
        SFup = SF + histoup.GetBinContent(xbin, ybin)
        SFdn = SF - histodn.GetBinContent(xbin, ybin)

        return [SF, SFup, SFdn]

