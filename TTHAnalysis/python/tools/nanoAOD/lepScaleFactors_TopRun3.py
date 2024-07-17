import correctionlib._core as core
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from CMGTools.TTHAnalysis.tools.nanoAOD.TopRun2UL_modules import ch
from copy import deepcopy
import ROOT as r
import os


class lepScaleFactors_TopRun3(Module):
    # =================== ### Main methods
    def __init__(self, year_, lepenvars = ["mu"]):
        self.basepath     = os.environ["CMSSW_BASE"] + "/src/CMGTools/TTHAnalysis/data/TopRun3/"
        self.basepathlep  = self.basepath + "/lepton/"
        self.basepathtrig = self.basepath + "/trigger/"
        self.year         = year_

        self.leptonSF = {}; self.triggerSF = {}
        self.leptonSF["m"] = {}; self.leptonSF["e"] = {}
        self.triggerSF[ch.ElMu] = {}; self.triggerSF[ch.Elec] = {}; self.triggerSF[ch.Muon] = {}
        for y in ["2022", "2022PostEE"]:
            self.leptonSF["m"][y] = {}
            self.leptonSF["e"][y] = {}
            for chan in [ch.ElMu, ch.Elec, ch.Muon]:
                self.triggerSF[chan][y] = {}

        self.systsLepEn = {}
        if len(lepenvars):
            for i, var in enumerate(lepenvars):
                self.systsLepEn[i+1]    = "_%sUp"%var
                self.systsLepEn[-(i+1)] = "_%sDown"%var

        print("[lepScaleFactors_TopRun2::constructor] - Loading histograms")
        
        self.jsonDictNames = {"": "", "_Up": "up", "_Dn": "down"}
        self.jsonDictNamesMu = {"": "nominal", "_statUp": "stat", "_statDn": "stat", "_systUp": "syst", "_systDn": "syst"}
        
        ## Muon ID
        self.leptonSF["m"]["2022"]["evaluatorID"]    = core.CorrectionSet.from_file(self.basepathlep + "ScaleFactors_Muon_Z_ID_ISO_2022_schemaV2.json")["NUM_TightID_DEN_TrackerMuons"]
        self.leptonSF["m"]["2022PostEE"]["evaluatorID"]    = core.CorrectionSet.from_file(self.basepathlep + "ScaleFactors_Muon_Z_ID_ISO_2022_EE_schemaV2.json")["NUM_TightID_DEN_TrackerMuons"]
        ##self.leptonSF["m"]["2022PostEE"]["evaluatorID"]    = core.CorrectionSet.from_file(self.basepathlep + "ScaleFactors_Muon_trackerMuons_Z_2022EE_Prompt_ID_ISO_schemaV2.json")["NUM_TightID_DEN_TrackerMuons"] #old onlyFG
        
        # Muon iso
        self.leptonSF["m"]["2022"]["evaluatorISO"]    = core.CorrectionSet.from_file(self.basepathlep + "ScaleFactors_Muon_Z_ID_ISO_2022_schemaV2.json")["NUM_TightPFIso_DEN_TightID"]
        self.leptonSF["m"]["2022PostEE"]["evaluatorISO"]    = core.CorrectionSet.from_file(self.basepathlep + "ScaleFactors_Muon_Z_ID_ISO_2022_EE_schemaV2.json")["NUM_TightPFIso_DEN_TightID"]
        ##self.leptonSF["m"]["2022PostEE"]["evaluatorISO"]    = core.CorrectionSet.from_file(self.basepathlep + "ScaleFactors_Muon_trackerMuons_Z_2022EE_Prompt_ID_ISO_schemaV2.json")["NUM_TightPFIso_DEN_TightID"] #old onlyFG

        # Elec ID-RECO
        self.leptonSF["e"]["2022"]["evaluatorID-RECO"]    = core.CorrectionSet.from_file(self.basepathlep + "electron_reRecoBCD.json")["Electron-ID-SF"]
        self.leptonSF["e"]["2022PostEE"]["evaluatorID-RECO"]    = core.CorrectionSet.from_file(self.basepathlep + "electron_reRecoE_promptFG.json")["Electron-ID-SF"]
        ##self.leptonSF["e"]["2022PostEE"]["evaluatorID-RECO"]    = core.CorrectionSet.from_file(self.basepathlep + "electronID_FG.json")["2022FG-Electron-ID-SF"] #old onlyFG
        self.leptonSF["e"]["2022"]["key"]       = "2022Re-recoBCD" 
        self.leptonSF["e"]["2022PostEE"]["key"] = "2022Re-recoE+PromptFG"
        ##self.leptonSF["e"]["2022PostEE"]["key"] = "2022FG"  #old onlyFG

        # Trigger elmu
        self.triggerSF[ch.ElMu]["2022"]    = self.loadHisto(self.basepathtrig + "triggerSFs_CD.root",        "h2D_SF_emu_lepABpt_FullError")
        self.triggerSF[ch.ElMu]["2022PostEE"]    = self.loadHisto(self.basepathtrig + "triggerSFs_EFG.root",        "h2D_SF_emu_lepABpt_FullError")
        ##self.triggerSF[ch.ElMu]["2022PostEE"]    = self.loadHisto(self.basepathtrig + "triggerSFs.root",        "h2D_SF_emu_lepABpt_FullError")
        # Trigger elel
        self.triggerSF[ch.Elec]["2022"]    = self.loadHisto(self.basepathtrig + "triggerSFs_CD.root",        "h2D_SF_ee_lepABpt_FullError")
        self.triggerSF[ch.Elec]["2022PostEE"]    = self.loadHisto(self.basepathtrig + "triggerSFs_EFG.root",        "h2D_SF_ee_lepABpt_FullError")
        ##self.triggerSF[ch.Elec]["2022PostEE"]    = self.loadHisto(self.basepathtrig + "triggerSFs.root",        "h2D_SF_ee_lepABpt_FullError")
        # Trigger mumu
        self.triggerSF[ch.Muon]["2022"]    = self.loadHisto(self.basepathtrig + "triggerSFs_CD.root",        "h2D_SF_mumu_lepABpt_FullError")
        self.triggerSF[ch.Muon]["2022PostEE"]    = self.loadHisto(self.basepathtrig + "triggerSFs_EFG.root",        "h2D_SF_mumu_lepABpt_FullError")
        ##self.triggerSF[ch.Muon]["2022PostEE"]    = self.loadHisto(self.basepathtrig + "triggerSFs.root",        "h2D_SF_mumu_lepABpt_FullError")


        print("[lepScaleFactors_TopRun2::constructor] - Finished loading histograms")
        return


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for var in ["", "_Up", "_Dn"]:
            self.out.branch('ElecIDSF'   + var, 'F')
            self.out.branch('ElecRECOSF' + var, 'F')
            self.out.branch('TrigSF'     + var, 'F')
            # Muon SFs with stat and syst variations no separated (ttbar Run3)
            #self.out.branch('MuonIDSF'   + var, 'F')
            #self.out.branch('MuonISOSF'  + var, 'F')
            if var == "":
                for delta,sys in self.systsLepEn.items():
                    self.out.branch('MuonIDSF'   + sys, 'F')
                    self.out.branch('MuonISOSF'  + sys, 'F')
                    self.out.branch('ElecIDSF'   + sys, 'F')
                    self.out.branch('ElecRECOSF' + sys, 'F')
                    self.out.branch('TrigSF'     + sys, 'F')

        for var in ["", "_statUp", "_statDn", "_systUp", "_systDn"]:
            self.out.branch('MuonIDSF'   + var, 'F')
            self.out.branch('MuonISOSF'  + var, 'F')

        return


    def analyze(self, event):
        year = self.year
        chan = event.channel
        run  = event.run
        leps = [l for l in Collection(event, "LepGood")]
        varleps = []

        # leptons
        # muons (separated into stat & syst)
        for var in ["", "_statUp", "_statDn", "_systUp", "_systDn"]:
            muonidsf = 1; muonisosf  = 1

            if event.nLepGood > 0:
                if abs(leps[0].pdgId) == 13: # muon
                    muon0IdSfNom = self.leptonSF["m"][year]["evaluatorID"].evaluate(abs(leps[0].eta), min(leps[0].pt_corrAll,199.), self.jsonDictNamesMu[""])
                    muon0IsoSfNom = self.leptonSF["m"][year]["evaluatorISO"].evaluate(abs(leps[0].eta), min(leps[0].pt_corrAll,199.), self.jsonDictNamesMu[""])
                    if var  == "":
                        muonidsf   *= muon0IdSfNom
                        muonisosf   *= muon0IsoSfNom
                    else: # The json contains the systematics variations that should be applied to the nominal SF
                        # We have to add 0.5% uncertainty to the nominal SF only for the systematic variations
                        if "syst" in var:
                            muon0IsoExtrpUnc = (muon0IsoSfNom * 0.005) * (muon0IsoSfNom * 0.005)
                        else:
                            muon0IsoExtrpUnc = 0
                        if "Up" in var:
                            muonidsf   *= (muon0IdSfNom + 
                                        self.leptonSF["m"][year]["evaluatorID"].evaluate(abs(leps[0].eta), min(leps[0].pt_corrAll,199.), self.jsonDictNamesMu[var]))
                            muonisosf *= (muon0IsoSfNom + 
                                        (self.leptonSF["m"][year]["evaluatorISO"].evaluate(abs(leps[0].eta), min(leps[0].pt_corrAll,199.), self.jsonDictNamesMu[var])**2 + muon0IsoExtrpUnc)**0.5)
                        elif "Dn" in var:
                            muonidsf   *= (muon0IdSfNom - 
                                        self.leptonSF["m"][year]["evaluatorID"].evaluate(abs(leps[0].eta), min(leps[0].pt_corrAll,199.), self.jsonDictNamesMu[var]))
                            muonisosf   *= (muon0IsoSfNom - 
                                        (self.leptonSF["m"][year]["evaluatorISO"].evaluate(abs(leps[0].eta), min(leps[0].pt_corrAll,199.), self.jsonDictNamesMu[var])**2 + muon0IsoExtrpUnc)**0.5)
                        else:
                            print("[lepScaleFactors_TopRun3::analyze] - ERROR: variation not recognized: ", var)                             

                if event.nLepGood > 1:
                    if abs(leps[1].pdgId) == 13: # muon
                        muon1IdSfNom = self.leptonSF["m"][year]["evaluatorID"].evaluate(abs(leps[1].eta), min(leps[1].pt_corrAll,199.), self.jsonDictNamesMu[""])
                        muon1IsoSfNom = self.leptonSF["m"][year]["evaluatorISO"].evaluate(abs(leps[1].eta), min(leps[1].pt_corrAll,199.), self.jsonDictNamesMu[""])
                        if var  == "":
                            muonidsf   *= muon1IdSfNom
                            muonisosf  *= muon1IsoSfNom
                        else: # The json contains the systematics variations that should be applied to the nominal SF
                            # We have to add 0.5% uncertainty to the nominal SF only for the systematic variations
                            if "syst" in var:
                                muon1IsoExtrpUnc = (muon1IsoSfNom * 0.005) * (muon1IsoSfNom * 0.005)
                            else:
                                muon1IsoExtrpUnc = 0
                            if "Up" in var:
                                muonidsf   *= (muon1IdSfNom + 
                                            self.leptonSF["m"][year]["evaluatorID"].evaluate(abs(leps[1].eta), min(leps[1].pt_corrAll,199.), self.jsonDictNamesMu[var]))
                                muonisosf *= (muon1IsoSfNom + 
                                            (self.leptonSF["m"][year]["evaluatorISO"].evaluate(abs(leps[1].eta), min(leps[1].pt_corrAll,199.), self.jsonDictNamesMu[var])**2 + muon1IsoExtrpUnc)**0.5)
                            elif "Dn" in var:
                                muonidsf   *= (muon1IdSfNom - 
                                            self.leptonSF["m"][year]["evaluatorID"].evaluate(abs(leps[1].eta), min(leps[1].pt_corrAll,199.), self.jsonDictNamesMu[var]))
                                muonisosf   *= (muon1IsoSfNom - 
                                            (self.leptonSF["m"][year]["evaluatorISO"].evaluate(abs(leps[1].eta), min(leps[1].pt_corrAll,199.), self.jsonDictNamesMu[var])**2 + muon1IsoExtrpUnc)**0.5)
                            else:
                                print("[lepScaleFactors_TopRun3::analyze] - ERROR: variation not recognized: ", var)  
            self.out.fillBranch('MuonIDSF'   + var, muonidsf)
            self.out.fillBranch('MuonISOSF'  + var, muonisosf)


        for var in ["", "_Up", "_Dn"]:
            # electrons (they are not separated into stat/syst)
            #### NOTE: electron SF are in terms of the eta of the associated supercluster and
            # the pT of the electron. Take also into account that deltaEtaSC = etaSC - eta, as
            # indicated in PhysicsTools/NanoAOD/python/electrons_cff.py.
            elecidsf = 1; elecrecosf = 1
            if event.nLepGood > 0:
                if abs(leps[0].pdgId) == 11: # electron
                    elec0IdSfNom = self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[""], "Tight", leps[0].eta + leps[0].deltaEtaSC, leps[0].pt_corrAll)
                    elec0IdSfExtrpUnc = (elec0IdSfNom * 0.01) * (elec0IdSfNom * 0.01)
                    if var == "":
                        elecidsf   *= elec0IdSfNom
                    elif var == "_Up":
                        elec0IdSfUp = self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[var], "Tight", leps[0].eta + leps[0].deltaEtaSC, leps[0].pt_corrAll)
                        elecidsf   *= (elec0IdSfNom + ((elec0IdSfUp - elec0IdSfNom)**2 + elec0IdSfExtrpUnc)**0.5)
                    elif var == "_Dn":
                        elec0IdSfDn = self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[var], "Tight", leps[0].eta + leps[0].deltaEtaSC, leps[0].pt_corrAll)
                        elecidsf   *= (elec0IdSfNom - ((elec0IdSfDn - elec0IdSfNom)**2 + elec0IdSfExtrpUnc)**0.5)                    
                    
                    if leps[0].pt_corrAll >= 20 and leps[0].pt_corrAll < 75:
                        elecrecosf *= self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[var], "Reco20to75", leps[0].eta + leps[0].deltaEtaSC, leps[0].pt_corrAll)
                    else:
                        elecrecosf *= self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[var], "RecoAbove75", leps[0].eta + leps[0].deltaEtaSC, leps[0].pt_corrAll)
                    #elecidsf   *= self.getLepSF(leps[0].pt_corrAll, leps[0].eta + leps[0].deltaEtaSC, var, "e", year, event, "id")
                    #elecrecosf *= self.getLepSF(leps[0].pt_corrAll, leps[0].eta + leps[0].deltaEtaSC, var, "e", year, event, "reco")
                    # Debug this, print electron pt, eta+deltaEtaSC, SF
                    #print(var)
                    #print("Electron pt: ", leps[0].pt_corrAll, " eta+deltaEtaSC: ", leps[0].eta + leps[0].deltaEtaSC, " SF: ", self.getLepSF(leps[0].pt_corrAll, leps[0].eta + leps[0].deltaEtaSC, var, "e", year, event, "id"))

                if event.nLepGood > 1:
                    if abs(leps[1].pdgId) == 11: # electron
                        #elecidsf   *= self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[var], "Tight", leps[1].eta + leps[1].deltaEtaSC, leps[1].pt_corrAll)
                        elec1IdSfNom = self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[""], "Tight", leps[1].eta + leps[1].deltaEtaSC, leps[1].pt_corrAll)
                        elec1IdSfExtrpUnc = (elec1IdSfNom * 0.01) * (elec1IdSfNom * 0.01)
                        if var == "":
                            elecidsf   *= elec1IdSfNom
                        elif var == "_Up":
                            elec1IdSfUp = self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[var], "Tight", leps[1].eta + leps[1].deltaEtaSC, leps[1].pt_corrAll)
                            elecidsf   *= (elec1IdSfNom + ((elec1IdSfUp - elec1IdSfNom)**2 + elec1IdSfExtrpUnc)**0.5)
                        elif var == "_Dn":
                            elec1IdSfDn = self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[var], "Tight", leps[1].eta + leps[1].deltaEtaSC, leps[1].pt_corrAll)
                            elecidsf   *= (elec1IdSfNom - ((elec1IdSfDn - elec1IdSfNom)**2 + elec1IdSfExtrpUnc)**0.5)     

                        if leps[1].pt_corrAll >= 20 and leps[1].pt_corrAll < 75:
                            elecrecosf *= self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[var], "Reco20to75", leps[1].eta + leps[1].deltaEtaSC, leps[1].pt_corrAll)
                        else:
                            elecrecosf *= self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[var], "RecoAbove75", leps[1].eta + leps[1].deltaEtaSC, leps[1].pt_corrAll)
                        #elecidsf   *= self.getLepSF(leps[1].pt_corrAll, leps[1].eta + leps[1].deltaEtaSC, var, "e", year, event, "id")
                        #elecrecosf *= self.getLepSF(leps[1].pt_corrAll, leps[1].eta + leps[1].deltaEtaSC, var, "e", year, event, "reco")


            self.out.fillBranch('ElecIDSF'   + var, elecidsf)
            self.out.fillBranch('ElecRECOSF' + var, elecrecosf)

            # triggers
            trigsf = 1
            #if len(leps) > 1 and chan != ch.NoChan:
            if event.nLepGood > 1 and chan != ch.NoChan:
                if chan == ch.ElMu:
                    trigsf *= self.getTrigSFdependingOnFlavour((leps[0].pt_corrAll, leps[0].pdgId),
                                                               (leps[1].pt_corrAll, leps[1].pdgId),
                                                               var, chan, year, ev = event)
                else:
                    trigsf *= self.getTrigSF(leps[0].pt_corrAll, leps[1].pt_corrAll, var, chan, year, ev = event)

            self.out.fillBranch('TrigSF' + var, trigsf)

            # Now for the lepton energy variations...
            if var == "":
                for delta,sys in self.systsLepEn.items():
                    # leptons
                    muonidsf = 1; muonisosf  = 1
                    elecidsf = 1; elecrecosf = 1
                    varleps = [l for l in Collection(event, "LepGood" + sys[1:])]

                    if len(varleps) > 0:
                        if   abs(varleps[0].pdgId) == 11: # electron
                            #elecidsf   *= self.getLepSF(getattr(varleps[0], "pt" + sys), varleps[0].eta + varleps[0].deltaEtaSC, var, "e", year, event, "id")
                            elecidsf   *= self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[var], "Tight", varleps[0].eta + varleps[0].deltaEtaSC, getattr(varleps[0], "pt" + sys))
                            #elecrecosf *= self.getLepSF(getattr(varleps[0], "pt" + sys), varleps[0].eta + varleps[0].deltaEtaSC, var, "e", year, event, "reco")
                            if getattr(varleps[0], "pt" + sys) >= 20 and getattr(varleps[0], "pt" + sys) < 75:
                                elecrecosf *= self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[var], "Reco20to75", varleps[0].eta + varleps[0].deltaEtaSC, getattr(varleps[0], "pt" + sys))
                            else:
                                elecrecosf *= self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[var], "RecoAbove75", varleps[0].eta + varleps[0].deltaEtaSC, getattr(varleps[0], "pt" + sys))
                        elif abs(varleps[0].pdgId) == 13: # muon
                            #muonidsf  *= self.getLepSF(getattr(varleps[0], "pt" + sys), varleps[0].eta, var, "m", year, event, "id")
                            #muonisosf *= self.getLepSF(getattr(varleps[0], "pt" + sys), varleps[0].eta, var, "m", year, event, "iso")
                            muonidsf   *= self.leptonSF["m"][year]["evaluatorID"].evaluate(abs(varleps[0].eta), min(getattr(varleps[0], "pt" + sys), 199.), self.jsonDictNamesMu[var])
                            muonisosf   *= self.leptonSF["m"][year]["evaluatorISO"].evaluate(abs(varleps[0].eta), min(getattr(varleps[0], "pt" + sys), 199.), self.jsonDictNamesMu[var])            

                        if len(varleps) > 1:
                            if   abs(varleps[1].pdgId) == 11: # electron
                                elecidsf   *= self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[var], "Tight", varleps[1].eta + varleps[1].deltaEtaSC, getattr(varleps[1], "pt" + sys))
                                #elecidsf   *= self.getLepSF(getattr(varleps[1], "pt" + sys), varleps[1].eta + varleps[1].deltaEtaSC, var, "e", year, event, "id")
                                #elecrecosf *= self.getLepSF(getattr(varleps[1], "pt" + sys), varleps[1].eta + varleps[1].deltaEtaSC, var, "e", year, event, "reco")
                                if getattr(varleps[1], "pt" + sys) >= 20 and getattr(varleps[1], "pt" + sys) < 75:
                                    elecrecosf *= self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[var], "Reco20to75", varleps[1].eta + varleps[1].deltaEtaSC, getattr(varleps[1], "pt" + sys))
                                else:
                                    elecrecosf *= self.leptonSF["e"][year]["evaluatorID-RECO"].evaluate(self.leptonSF["e"][year]["key"], "sf" + self.jsonDictNames[var], "RecoAbove75", varleps[1].eta + varleps[1].deltaEtaSC, getattr(varleps[1], "pt" + sys))
                            elif abs(varleps[1].pdgId) == 13: # muon
                                #muonidsf  *= self.getLepSF(getattr(varleps[1], "pt" + sys), varleps[1].eta, var, "m", year, event, "id")
                                #muonisosf *= self.getLepSF(getattr(varleps[1], "pt" + sys), varleps[1].eta, var, "m", year, event, "iso")
                                muonidsf   *= self.leptonSF["m"][year]["evaluatorID"].evaluate(abs(varleps[1].eta), min(getattr(varleps[1], "pt" + sys), 199.), self.jsonDictNamesMu[var])
                                muonisosf   *= self.leptonSF["m"][year]["evaluatorISO"].evaluate(abs(varleps[1].eta), min(getattr(varleps[1], "pt" + sys), 199.), self.jsonDictNamesMu[var])

                    self.out.fillBranch('MuonIDSF'   + sys, muonidsf)
                    self.out.fillBranch('MuonISOSF'  + sys, muonisosf)
                    self.out.fillBranch('ElecIDSF'   + sys, elecidsf)
                    self.out.fillBranch('ElecRECOSF' + sys, elecrecosf)

                    # triggers
                    trigsf = 1
                    if len(varleps) > 1 and getattr(event, "channel" + sys) != ch.NoChan:
                        if getattr(event, "channel" + sys) == ch.ElMu:
                            trigsf *= self.getTrigSFdependingOnFlavour((getattr(varleps[0], "pt" + sys), varleps[0].pdgId),
                                                                       (getattr(varleps[1], "pt" + sys), varleps[1].pdgId),
                                                                       var, getattr(event, "channel" + sys),
                                                                       year, ev = event)
                        else:
                            trigsf *= self.getTrigSF(getattr(varleps[0], "pt" + sys),
                                                     getattr(varleps[1], "pt" + sys),
                                                     var, getattr(event, "channel" + sys),
                                                     year, ev = event)

                    self.out.fillBranch('TrigSF' + sys, trigsf)

        return True


    # =================== ### Other methods
    def loadHisto(self, fil, hist):
        tf = r.TFile.Open(fil)
        if not tf: raise RuntimeError("[lepScaleFactors_TopRun2::loadHisto] FATAL: no such file %s"%fil)
        hist = tf.Get(hist)
        if not hist: raise RuntimeError("[lepScaleFactors_TopRun2::loadHisto] FATAL: no such object %s in %s"%(hist,fil))
        ret = deepcopy(hist)
        tf.Close()
        return ret


#    def loadHistoWithUncs(self, fil, hist, uncs = []):
#        tf = r.TFile.Open(fil)
#        if not tf: raise RuntimeError("[lepScaleFactors_TopRun2::loadHistoWithVars] FATAL: no such file %s"%fil)
#        histnom = tf.Get(hist)
#        if not hist: raise RuntimeError("[lepScaleFactors_TopRun2::loadHistoWithVars] FATAL: no such object %s in %s"%(hist,fil))
#        ret = [deepcopy(histnom)]
#
#        for var in uncs:
#            tmphist = tf.Get(hist + var)
#            ret.append(deepcopy(tmphist))
#        tf.Close()
#        return ret


    def getTrigSF(self, pt1, pt2, var, fl, yr, ev = None):
        # NOTE: it will always use as the first lepton the one with largest pt
        hist = self.triggerSF[fl][yr]
        tmppt1 = min(max([pt1, pt2]), 349.); tmppt2 = min(min([pt1, pt2]), 349.)
        pt1bin = max(1, min(hist.GetNbinsX(), hist.GetXaxis().FindBin(tmppt1)))
        pt2bin = max(1, min(hist.GetNbinsY(), hist.GetYaxis().FindBin(tmppt2)))

        #if ev.event == 467503 or ev.event == 486099:
            #print "\n======================================"

            #print "pt1:", pt1,       "pt2:", pt2
            #print "tmppt1:", tmppt1, "tmppt2:", tmppt2
            #print "pt1bin:", pt1bin, "pt2bin:", pt2bin
            #print "histname:", hist.GetName()
            #print "yr:", yr
            #print "dict:", self.triggerSF

            #print "======================================\n"

        out = hist.GetBinContent(pt1bin, pt2bin)
        if   "up" in var.lower():
            out += hist.GetBinError(pt1bin, pt2bin)
        elif "dn" in var.lower():
            out -= hist.GetBinError(pt1bin, pt2bin)

        return out


    def getTrigSFdependingOnFlavour(self, l1tp, l2tp, var, fl, yr, ev = None):
        # NOTE: X -> electron, Y -> muon
        hist = self.triggerSF[fl][yr]

        if abs(l1tp[1]) == 11:
            elpt = min(l1tp[0], 349.)
            mupt = min(l2tp[0], 349.)
        else:
            elpt = min(l2tp[0], 349.)
            mupt = min(l1tp[0], 349.)
        pt1bin = max(1, min(hist.GetNbinsX(), hist.GetXaxis().FindBin(elpt)))
        pt2bin = max(1, min(hist.GetNbinsY(), hist.GetYaxis().FindBin(mupt)))

        out = hist.GetBinContent(pt1bin, pt2bin)
        if   "up" in var.lower():
            out += hist.GetBinError(pt1bin, pt2bin)
        elif "dn" in var.lower():
            out -= hist.GetBinError(pt1bin, pt2bin)

        return out

#    def getEfficiency(self, histo, pt, eta):
#        '''
#        This method returns the efficiency for a given pt and eta, ttbar run3
#        '''
#        eta = abs(eta)
#        ptbin = max(1, min(histo.GetNbinsX(), histo.GetXaxis().FindBin(pt)))
#        etabin = max(1, min(histo.GetNbinsY(), histo.GetYaxis().FindBin(eta)))
#        eff = histo.GetBinContent(ptbin, etabin)
#        err = histo.GetBinError(ptbin, etabin)
#        return eff, err

    
#    def getTrigSFfromSingleLepton(self, mupt, mueta, ept, eeta, var, fl, yr, ev = None):
#        '''
#        It computes the SFs from the single lepton trigger efficiencies
#        eff_data = mu_eff_data + e_eff_data - mu_eff_data*e_eff_data
#        eff_MC = mu_eff_MC + e_eff_MC - mu_eff_MC*e_eff_MC
#        SF = eff_data/eff_MC
#
#        We separate the errors in electron and muon triggers
#        SF_eUp=(effm_data + effe_data+er_effe_data - effm_data*(er_effe_data+effe_data)) /(effm_MC + effe_MC+er_effe_MC - effm_MC*(er_effe_MC+effe_MC)) 
#        SF_e_Down = (effm_data + effe_data-er_effe_data - effm_data*(er_effe_data-effe_data)) /(effm_MC + effe_MC-er_effe_MC - effm_MC*(er_effe_MC-effe_MC)) 
#        SF_mUp= (effm_data+er_effm_data + effe_data - (effm_data+er_effm_data)*effe_data) /(effm_MC+er_effm_MC + effe_MC - (effm_MC+er_effm_MC)*effe_MC) 
#        SF_mDo=(effm_data-er_effm_data + effe_data - (effm_data-er_effm_data)*effe_data) /(effm_MC-er_effm_MC + effe_MC - (effm_MC-er_effm_MC)*effe_MC) 
#        '''
#        hist_mu_Data = self.triggerSF["m"][yr]["effData"]
#        hist_mu_MC   = self.triggerSF["m"][yr]["effMC"]
#        hist_e_Data = self.triggerSF["e"][yr]["effData"]
#        hist_e_MC   = self.triggerSF["e"][yr]["effMC"]
#
#        if fl == ch.ElMu:
#            eff_Data = self.getEfficiency(hist_mu_Data, mupt, mueta)[0] + self.getEfficiency(hist_e_Data, ept, eeta)[0] - self.getEfficiency(hist_mu_Data, mupt, mueta)[0]*self.getEfficiency(hist_e_Data, ept, eeta)[0]
#            eff_MC = self.getEfficiency(hist_mu_MC, mupt, mueta)[0] + self.getEfficiency(hist_e_MC, ept, eeta)[0] - self.getEfficiency(hist_mu_MC, mupt, mueta)[0]*self.getEfficiency(hist_e_MC, ept, eeta)[0]
#            out = eff_Data/eff_MC
#            if   "eup" in var.lower():
#                eff_Data = self.getEfficiency(hist_mu_Data, mupt, mueta)[0] + self.getEfficiency(hist_e_Data, ept, eeta)[0] + self.getEfficiency(hist_e_Data, ept, eeta)[1] - self.getEfficiency(hist_mu_Data, mupt, mueta)[0]*(self.getEfficiency(hist_e_Data, ept, eeta)[1]+self.getEfficiency(hist_e_Data, ept, eeta)[0])
#                eff_MC = self.getEfficiency(hist_mu_MC, mupt, mueta)[0] + self.getEfficiency(hist_e_MC, ept, eeta)[0] + self.getEfficiency(hist_e_MC, ept, eeta)[1] - self.getEfficiency(hist_mu_MC, mupt, mueta)[0]*(self.getEfficiency(hist_e_MC, ept, eeta)[1]+self.getEfficiency(hist_e_MC, ept, eeta)[0])
#                out = eff_Data/eff_MC
#            elif "edn" in var.lower():
#                eff_Data = self.getEfficiency(hist_mu_Data, mupt, mueta)[0] + self.getEfficiency(hist_e_Data, ept, eeta)[0] - self.getEfficiency(hist_e_Data, ept, eeta)[1] - self.getEfficiency(hist_mu_Data, mupt, mueta)[0]*(self.getEfficiency(hist_e_Data, ept, eeta)[1]-self.getEfficiency(hist_e_Data, ept, eeta)[0])
#                eff_MC = self.getEfficiency(hist_mu_MC, mupt, mueta)[0] + self.getEfficiency(hist_e_MC, ept, eeta)[0] - self.getEfficiency(hist_e_MC, ept, eeta)[1] - self.getEfficiency(hist_mu_MC, mupt, mueta)[0]*(self.getEfficiency(hist_e_MC, ept, eeta)[1]-self.getEfficiency(hist_e_MC, ept, eeta)[0])
#                out = eff_Data/eff_MC
#            elif "mup" in var.lower():
#                eff_Data = self.getEfficiency(hist_mu_Data, mupt, mueta)[0] + self.getEfficiency(hist_mu_Data, mupt, mueta)[1] + self.getEfficiency(hist_e_Data, ept, eeta)[0] - (self.getEfficiency(hist_mu_Data, mupt, mueta)[0]+self.getEfficiency(hist_mu_Data, mupt, mueta)[1])*self.getEfficiency(hist_e_Data, ept, eeta)[0]
#                eff_MC = self.getEfficiency(hist_mu_MC, mupt, mueta)[0] + self.getEfficiency(hist_mu_MC, mupt, mueta)[1] + self.getEfficiency(hist_e_MC, ept, eeta)[0] - (self.getEfficiency(hist_mu_MC, mupt, mueta)[0]+self.getEfficiency(hist_mu_MC, mupt, mueta)[1])*self.getEfficiency(hist_e_MC, ept, eeta)[0]
#                out = eff_Data/eff_MC
#            elif "mdn" in var.lower():
#                eff_Data = self.getEfficiency(hist_mu_Data, mupt, mueta)[0] - self.getEfficiency(hist_mu_Data, mupt, mueta)[1] + self.getEfficiency(hist_e_Data, ept, eeta)[0] - (self.getEfficiency(hist_mu_Data, mupt, mueta)[0]-self.getEfficiency(hist_mu_Data, mupt, mueta)[1])*self.getEfficiency(hist_e_Data, ept, eeta)[0]
#                eff_MC = self.getEfficiency(hist_mu_MC, mupt, mueta)[0] - self.getEfficiency(hist_mu_MC, mupt, mueta)[1] + self.getEfficiency(hist_e_MC, ept, eeta)[0] - (self.getEfficiency(hist_mu_MC, mupt, mueta)[0]-self.getEfficiency(hist_mu_MC, mupt, mueta)[1])*self.getEfficiency(hist_e_MC, ept, eeta)[0]
#                out = eff_Data/eff_MC
#            return out
    
                


#    def getFromHisto(self, pt, eta, histo, err = False):
#        #### NOTE: this method obtains the information for a given pt and eta. In the
#        # case, like e.g. electron SF, where the histogram is in terms of the supercluster
#        # pseudorapidity rather than the final pseudorapidity, this script will not
#        # automatically change it.
#        normal = True # normal := {x: pt, y: eta}
#        if (histo.GetXaxis().GetBinUpEdge(histo.GetNbinsX())) < 10: normal = False
#
#        tmpeta = eta
#        #tmpeta = abs(eta)
#
#        if    ((normal     and histo.GetYaxis().GetBinLowEdge(1) >= 0)
#            or (not normal and histo.GetXaxis().GetBinLowEdge(1) >= 0)):
#            tmpeta = abs(tmpeta)
#
#        if normal:
#            xbin = max(1, min(histo.GetNbinsX(), histo.GetXaxis().FindBin(pt)))
#            ybin = max(1, min(histo.GetNbinsY(), histo.GetYaxis().FindBin(tmpeta)))
#        else:
#            xbin = max(1, min(histo.GetNbinsX(), histo.GetXaxis().FindBin(tmpeta)))
#            ybin = max(1, min(histo.GetNbinsY(), histo.GetYaxis().FindBin(pt)))
#
#        if err: return histo.GetBinError(xbin, ybin)
#        else:   return histo.GetBinContent(xbin, ybin)
#        return



#    def getLepSF(self, pt, eta, var, fl, yr, ev = None, typ = None):
#        out = 1
#        doneSFtype = []
#
#        for sf in self.leptonSF[fl][yr]:
#            if typ:
#                if typ not in sf.lower(): continue
#            
#            if sf not in doneSFtype:
#                out  *= self.getFromHisto(pt, eta, self.leptonSF[fl][yr][sf])
#                doneSFtype.append(sf)
#        
#        if var != "":
#            doneSFtype = []
#            tmpsum = 0
#            for sf in self.leptonSF[fl][yr]:
#                if typ:
#                    if typ not in sf.lower(): continue
#
#                if sf not in doneSFtype:
#                    tmpsum += (self.getFromHisto(pt, eta, self.leptonSFuncs[fl][yr][sf + var[:-2]] if (sf + var[:-2]) in self.leptonSFuncs[fl][yr] else self.leptonSF[fl][yr][sf], err = True))**2
#                    doneSFtype.append(sf)
#
##            if fl == "m" and "syst" in var:
##                if typ:
##                    if "iso" in typ:
##                        tmpsum += (out * 0.005)**2 # Extrapolation to top-like from DY-like phase space (CMS-AN-"2018"-210). Only for muons.
##                else:
##                    tmpsum += (out * 0.005)**2 # Extrapolation to top-like from DY-like phase space (CMS-AN-"2018"-210). Only for muons.
#            
#            if fl == "m":
#                if typ:
#                    tmpsum += (out * 0.005)**2 # Extrapolation to top-like from DY-like phase space (CMS-AN-"2018"-210). Only for muons.
#                else:
#                    tmpsum += (out * 0.005)**2 # Extrapolation to top-like from DY-like phase space (CMS-AN-"2018"-210). Only for muons.
#            if fl == "e":
#                if typ:
#                    tmpsum += (out * 0.01)**2 # Extrapolation to top-like from DY-like phase space (CMS-AN-"2018"-210). Only for electrons.
#                else:
#                    tmpsum += (out * 0.01)**2 # Extrapolation to top-like from DY-like phase space (CMS-AN-"2018"-210). Only for electrons.
#            tmpsum = r.TMath.Sqrt(tmpsum)
#
#            if "up" in var.lower(): out += tmpsum
#            else:                   out -= tmpsum
#        return out
