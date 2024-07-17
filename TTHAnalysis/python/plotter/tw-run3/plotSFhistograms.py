import os, sys, argparse, enum
from multiprocessing import Pool
from copy import deepcopy
import warnings as wr
import ROOT as r
import correctionlib._core as core
import numpy as np

sys.path.append('{cmsswpath}/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/differential/'.format(cmsswpath = os.environ['CMSSW_BASE']))
import tdrstyle
from array import array

r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(True)

class ch(enum.IntEnum):
    NoChan = 0
    ElMu   = 1
    Muon   = 2
    Elec   = 3
    ElMuFromTaus = 4
    MuonFromTaus = 5
    ElecFromTaus = 6
    ElMuMixedFromTaus = 7
    MuonMixedFromTaus = 8
    ElecMixedFromTaus = 9


def loadHisto(fil, hist):
    tf = r.TFile.Open(fil)
    if not tf: raise RuntimeError("[lepScaleFactors_TopRun2::loadHisto] FATAL: no such file %s"%fil)
    hist = tf.Get(hist)
    if not hist: raise RuntimeError("[lepScaleFactors_TopRun2::loadHisto] FATAL: no such object %s in %s"%(hist,fil))
    ret = deepcopy(hist)
    tf.Close()
    return ret

def getHistoFromJson(jsonDict, binEta, binPt, dictNames, isMuonSF = True, whichSF="", electronCorrName = ""):
    '''
    This Function takes a json dictionary and returns a TH2F histogram
    Cada Json es de su padre y su madre asi que hay que distinguir a la hora de llamar a la funcion
    '''
    thehisto = r.TH2F("thehisto", "", len(binEta) - 1, array("f", binEta), len(binPt) - 1, array("f", binPt))
    for iB in range(1, thehisto.GetNbinsX() + 1):
        for jB in range(1, thehisto.GetNbinsY() + 1):
            if isMuonSF:
                thehisto.SetBinContent(iB, jB, jsonDict.evaluate(float(thehisto.GetXaxis().GetBinCenter(iB)), float(thehisto.GetYaxis().GetBinCenter(jB)), dictNames[""]))
                thehisto.SetBinError(iB, jB, 
                    (((jsonDict.evaluate(float(thehisto.GetXaxis().GetBinCenter(iB)), float(thehisto.GetYaxis().GetBinCenter(jB)), dictNames["_statUp"]))**2 +
                    (jsonDict.evaluate(float(thehisto.GetXaxis().GetBinCenter(iB)), float(thehisto.GetYaxis().GetBinCenter(jB)), dictNames["_systUp"]))**2)**0.5 + 
                    ((jsonDict.evaluate(float(thehisto.GetXaxis().GetBinCenter(iB)), float(thehisto.GetYaxis().GetBinCenter(jB)), dictNames["_statDn"]))**2 +
                    (jsonDict.evaluate(float(thehisto.GetXaxis().GetBinCenter(iB)), float(thehisto.GetYaxis().GetBinCenter(jB)), dictNames["_systDn"]))**2)**0.5)/2.)
            else:
                thehisto.SetBinContent(iB, jB, jsonDict.evaluate(electronCorrName,"sf"+dictNames[""],whichSF,float(thehisto.GetXaxis().GetBinCenter(iB)), float(thehisto.GetYaxis().GetBinCenter(jB))))
                # Here we only have the up and down variations
                thehisto.SetBinError(iB, jB,
                    (abs(jsonDict.evaluate(electronCorrName,"sf"+dictNames[""],whichSF,float(thehisto.GetXaxis().GetBinCenter(iB)), float(thehisto.GetYaxis().GetBinCenter(jB))) - 
                    jsonDict.evaluate(electronCorrName,"sf"+dictNames["_Up"],whichSF,float(thehisto.GetXaxis().GetBinCenter(iB)), float(thehisto.GetYaxis().GetBinCenter(jB)))) +
                    abs(jsonDict.evaluate(electronCorrName,"sf"+dictNames[""],whichSF,float(thehisto.GetXaxis().GetBinCenter(iB)), float(thehisto.GetYaxis().GetBinCenter(jB))) -
                    jsonDict.evaluate(electronCorrName,"sf"+dictNames["_Dn"],whichSF,float(thehisto.GetXaxis().GetBinCenter(iB)), float(thehisto.GetYaxis().GetBinCenter(jB)))))/2.
                )         
    return thehisto

#def loadHistoWithUncs(fil, hist, uncs = []):
#    tf = r.TFile.Open(fil)
#    if not tf: raise RuntimeError("[lepScaleFactors_TopRun2::loadHistoWithVars] FATAL: no such file %s"%fil)
#    histnom = tf.Get(hist)
#    if not hist: raise RuntimeError("[lepScaleFactors_TopRun2::loadHistoWithVars] FATAL: no such object %s in %s"%(hist,fil))
#    ret = [deepcopy(histnom)]
#
#    for var in uncs:
#        tmphist = tf.Get(hist + var)
#        ret.append(deepcopy(tmphist))
#    tf.Close()
#    return ret

sfpath       = os.environ['CMSSW_BASE'] + "/src/CMGTools/TTHAnalysis/data/TopRun3/"
basepathlep  = sfpath + "/lepton/"
basepathtrig = sfpath + "/trigger/"
SFdict = {}
SFdict["lepton"] = {}; SFdict["lepton"] = {}; SFdict["trigger"] = {}
SFdict["lepton"]["m"] = {}; SFdict["lepton"]["e"] = {}
SFdict["trigger"][ch.ElMu] = {}; SFdict["trigger"][ch.Elec] = {}; SFdict["trigger"][ch.Muon] = {}
for y in ["2022", "2022PostEE"]:
    SFdict["lepton"]["m"][y] = {}
    SFdict["lepton"]["e"][y] = {}
    for chan in [ch.ElMu, ch.Elec, ch.Muon]:
        SFdict["trigger"][chan][y] = {}

jsonDictNames = {"": "", "_Up": "up", "_Dn": "down"}
jsonDictNamesMu = {"": "nominal", "_statUp": "stat", "_statDn": "stat", "_systUp": "syst", "_systDn": "syst"}

## Muon ID
#SFdict["lepton"]["m"]["2022"]["idtight"]    = getHistoFromJson(core.CorrectionSet.from_file(basepathlep + "ScaleFactors_Muon_trackerMuons_Z_2022EE_Prompt_ID_ISO_schemaV2.json")["NUM_TightID_DEN_TrackerMuons"], binEta = [0.0,0.9,1.2,2.1,2.4], binPt = [15.0,20.0,25.0,30.0,40.0,50.0,60.0,120.0,200.0],dictNames=jsonDictNamesMu, isMuonSF = True)
SFdict["lepton"]["m"]["2022"]["idtight"]    = getHistoFromJson(core.CorrectionSet.from_file(basepathlep + "ScaleFactors_Muon_Z_ID_ISO_2022_schemaV2.json")["NUM_TightID_DEN_TrackerMuons"], binEta = [0.0,0.9,1.2,2.1,2.4], binPt = [15.0,20.0,25.0,30.0,40.0,50.0,60.0,120.0,200.0],dictNames=jsonDictNamesMu, isMuonSF = True)
SFdict["lepton"]["m"]["2022PostEE"]["idtight"]    = getHistoFromJson(core.CorrectionSet.from_file(basepathlep + "ScaleFactors_Muon_Z_ID_ISO_2022_EE_schemaV2.json")["NUM_TightID_DEN_TrackerMuons"], binEta = [0.0,0.9,1.2,2.1,2.4], binPt = [15.0,20.0,25.0,30.0,40.0,50.0,60.0,120.0,200.0],dictNames=jsonDictNamesMu, isMuonSF = True)

# Muon iso
#SFdict["lepton"]["m"]["2022"]["iso"]    = getHistoFromJson(core.CorrectionSet.from_file(basepathlep + "ScaleFactors_Muon_trackerMuons_Z_2022EE_Prompt_ID_ISO_schemaV2.json")["NUM_TightPFIso_DEN_TightID"], binEta = [0.0,0.9,1.2,2.1,2.4], binPt = [15.0,20.0,25.0,30.0,40.0,50.0,60.0,120.0,200.0],dictNames=jsonDictNamesMu, isMuonSF = True)
SFdict["lepton"]["m"]["2022"]["iso"]    = getHistoFromJson(core.CorrectionSet.from_file(basepathlep + "ScaleFactors_Muon_Z_ID_ISO_2022_schemaV2.json")["NUM_TightPFIso_DEN_TightID"], binEta = [0.0,0.9,1.2,2.1,2.4], binPt = [15.0,20.0,25.0,30.0,40.0,50.0,60.0,120.0,200.0],dictNames=jsonDictNamesMu, isMuonSF = True)
SFdict["lepton"]["m"]["2022PostEE"]["iso"]    = getHistoFromJson(core.CorrectionSet.from_file(basepathlep + "ScaleFactors_Muon_Z_ID_ISO_2022_EE_schemaV2.json")["NUM_TightPFIso_DEN_TightID"], binEta = [0.0,0.9,1.2,2.1,2.4], binPt = [15.0,20.0,25.0,30.0,40.0,50.0,60.0,120.0,200.0],dictNames=jsonDictNamesMu, isMuonSF = True)

# Elec ID
#SFdict["lepton"]["e"]["2022"]["idtight"]    = getHistoFromJson(core.CorrectionSet.from_file(basepathlep + "electronID_FG.json")["2022FG-Electron-ID-SF"], binEta = [-2.4,-2.0,-1.566,-1.444,-0.8,0.0,0.8,1.444,1.566,2.0,2.4], binPt = [10.0,20.0,35.0,50.0,100.0,200.0],dictNames=jsonDictNames, isMuonSF = False, whichSF = "Tight")
SFdict["lepton"]["e"]["2022"]["idtight"]    = getHistoFromJson(core.CorrectionSet.from_file(basepathlep + "electron_reRecoBCD.json")["Electron-ID-SF"], binEta = [-2.4,-2.0,-1.566,-1.444,-0.8,0.0,0.8,1.444,1.566,2.0,2.4], binPt = [10.0,20.0,35.0,50.0,100.0,200.0],dictNames=jsonDictNames, isMuonSF = False, whichSF = "Tight", electronCorrName = "2022Re-recoBCD")
SFdict["lepton"]["e"]["2022PostEE"]["idtight"]    = getHistoFromJson(core.CorrectionSet.from_file(basepathlep + "electron_reRecoE_promptFG.json")["Electron-ID-SF"], binEta = [-2.4,-2.0,-1.566,-1.444,-0.8,0.0,0.8,1.444,1.566,2.0,2.4], binPt = [10.0,20.0,35.0,50.0,100.0,200.0],dictNames=jsonDictNames, isMuonSF = False, whichSF = "Tight", electronCorrName = "2022Re-recoE+PromptFG")

# Elec reco
#SFdict["lepton"]["e"]["2022"]["Reco20to75"]    = getHistoFromJson(core.CorrectionSet.from_file(basepathlep + "electronID_FG.json")["2022FG-Electron-ID-SF"], binEta = [-2.4,-2.0,-1.566,-1.444,-0.8,0.0,0.8,1.444,1.566,2.0,2.4], binPt = [20.0,45.0,75.0],dictNames=jsonDictNames, isMuonSF = False, whichSF = "Reco20to75")
#SFdict["lepton"]["e"]["2022"]["RecoAbove75"]    = getHistoFromJson(core.CorrectionSet.from_file(basepathlep + "electronID_FG.json")["2022FG-Electron-ID-SF"], binEta = [-2.4,-2.0,-1.566,-1.444,-0.8,0.0,0.8,1.444,1.566,2.0,2.4], binPt = [75.0,100.0,200.0],dictNames=jsonDictNames, isMuonSF = False, whichSF = "RecoAbove75")
SFdict["lepton"]["e"]["2022"]["Reco20to75"]    = getHistoFromJson(core.CorrectionSet.from_file(basepathlep + "electron_reRecoBCD.json")["Electron-ID-SF"], binEta = [-2.4,-2.0,-1.566,-1.444,-0.8,0.0,0.8,1.444,1.566,2.0,2.4], binPt = [20.0,45.0,75.0],dictNames=jsonDictNames, isMuonSF = False, whichSF = "Reco20to75", electronCorrName = "2022Re-recoBCD")
SFdict["lepton"]["e"]["2022"]["RecoAbove75"]    = getHistoFromJson(core.CorrectionSet.from_file(basepathlep + "electron_reRecoBCD.json")["Electron-ID-SF"], binEta = [-2.4,-2.0,-1.566,-1.444,-0.8,0.0,0.8,1.444,1.566,2.0,2.4], binPt = [75.0,100.0,200.0],dictNames=jsonDictNames, isMuonSF = False, whichSF = "RecoAbove75", electronCorrName = "2022Re-recoBCD")
SFdict["lepton"]["e"]["2022PostEE"]["Reco20to75"]    = getHistoFromJson(core.CorrectionSet.from_file(basepathlep + "electron_reRecoE_promptFG.json")["Electron-ID-SF"], binEta = [-2.4,-2.0,-1.566,-1.444,-0.8,0.0,0.8,1.444,1.566,2.0,2.4], binPt = [20.0,45.0,75.0],dictNames=jsonDictNames, isMuonSF = False, whichSF = "Reco20to75", electronCorrName = "2022Re-recoE+PromptFG")
SFdict["lepton"]["e"]["2022PostEE"]["RecoAbove75"]    = getHistoFromJson(core.CorrectionSet.from_file(basepathlep + "electron_reRecoE_promptFG.json")["Electron-ID-SF"], binEta = [-2.4,-2.0,-1.566,-1.444,-0.8,0.0,0.8,1.444,1.566,2.0,2.4], binPt = [75.0,100.0,200.0],dictNames=jsonDictNames, isMuonSF = False, whichSF = "RecoAbove75", electronCorrName = "2022Re-recoE+PromptFG")

# Trigger elmu
#SFdict["trigger"][ch.ElMu]["2022"] = loadHisto(basepathtrig + "triggerSFs.root", "h2D_SF_emu_lepABpt_FullError")
SFdict["trigger"][ch.ElMu]["2022"] = loadHisto(basepathtrig + "triggerSFs_CD.root", "h2D_SF_emu_lepABpt_FullError")
SFdict["trigger"][ch.ElMu]["2022PostEE"] = loadHisto(basepathtrig + "triggerSFs_EFG.root", "h2D_SF_emu_lepABpt_FullError")

# Trigger elel
#SFdict["trigger"][ch.Elec]["2022"] = loadHisto(basepathtrig + "triggerSFs.root", "h2D_SF_ee_lepABpt_FullError")
SFdict["trigger"][ch.Elec]["2022"] = loadHisto(basepathtrig + "triggerSFs_CD.root", "h2D_SF_ee_lepABpt_FullError")
SFdict["trigger"][ch.Elec]["2022PostEE"] = loadHisto(basepathtrig + "triggerSFs_EFG.root", "h2D_SF_ee_lepABpt_FullError")

# Trigger mumu
#SFdict["trigger"][ch.Muon]["2022"] = loadHisto(basepathtrig + "triggerSFs.root", "h2D_SF_mumu_lepABpt_FullError")
SFdict["trigger"][ch.Muon]["2022"] = loadHisto(basepathtrig + "triggerSFs_CD.root", "h2D_SF_mumu_lepABpt_FullError")
SFdict["trigger"][ch.Muon]["2022PostEE"] = loadHisto(basepathtrig + "triggerSFs_EFG.root", "h2D_SF_mumu_lepABpt_FullError")

SFdict["btagging_deepjet"] = {}
SFdict["btaggingSF_deepjet"] = {}
SFdict["btagging_particlenet"] = {}
SFdict["btaggingSF_particlenet"] = {}
SFdict["btagging_particletransformer"] = {}
SFdict["btaggingSF_particletransformer"] = {}

f_eff_2022PostEE   = r.TFile.Open("temp_Run3_plots/2024_02_28_btagEff_tW_comb/2022PostEE/eff/output.root", "read")
f_eff_2022  = r.TFile.Open("temp_Run3_plots/2024_02_28_btagEff_tW_comb/2022/eff/output.root", "read")
f_eff_dict = {"2022PostEE": f_eff_2022PostEE, "2022": f_eff_2022}
for el in ["B", "C", "L"]:
    SFdict["btagging_deepjet"][el] = {}
    SFdict["btaggingSF_deepjet"][el] = {}
    SFdict["btagging_particlenet"][el] = {}
    SFdict["btaggingSF_particlenet"][el] = {}
    SFdict["btagging_particletransformer"][el] = {}
    SFdict["btaggingSF_particletransformer"][el] = {}
    for year in ["2022", "2022PostEE"]:
        f_eff = f_eff_dict[year]
        SFdict["btagging_deepjet"][el][year]   = deepcopy(f_eff.Get("btageff_deepjet_{t}_btag_{t}_tw".format(t = el)).Clone())
        SFdict["btaggingSF_deepjet"][el][year] = {}
        SFdict["btaggingSF_deepjet"][el][year][""]   = deepcopy(f_eff.Get("btagsf_deepjet_{t}_btag_pt_{t}_tw" .format(t = el)).Clone())
        SFdict["btaggingSF_deepjet"][el][year]["up"] = deepcopy(f_eff.Get("btagsfup_deepjet_{t}_btag_pt_{t}_tw" .format(t = el)).Clone())
        SFdict["btaggingSF_deepjet"][el][year]["dn"] = deepcopy(f_eff.Get("btagsfdn_deepjet_{t}_btag_pt_{t}_tw" .format(t = el)).Clone())
        SFdict["btagging_particlenet"][el][year]   = deepcopy(f_eff.Get("btageff_particlenet_{t}_btag_{t}_tw".format(t = el)).Clone())
        SFdict["btaggingSF_particlenet"][el][year] = {}
        SFdict["btaggingSF_particlenet"][el][year][""]   = deepcopy(f_eff.Get("btagsf_particlenet_{t}_btag_pt_{t}_tw" .format(t = el)).Clone())
        SFdict["btaggingSF_particlenet"][el][year]["up"] = deepcopy(f_eff.Get("btagsfup_particlenet_{t}_btag_pt_{t}_tw" .format(t = el)).Clone())
        SFdict["btaggingSF_particlenet"][el][year]["dn"] = deepcopy(f_eff.Get("btagsfdn_particlenet_{t}_btag_pt_{t}_tw" .format(t = el)).Clone())
        SFdict["btagging_particletransformer"][el][year]   = deepcopy(f_eff.Get("btageff_particletransformer_{t}_btag_{t}_tw".format(t = el)).Clone())
        SFdict["btaggingSF_particletransformer"][el][year] = {}
        SFdict["btaggingSF_particletransformer"][el][year][""]   = deepcopy(f_eff.Get("btagsf_particletransformer_{t}_btag_pt_{t}_tw" .format(t = el)).Clone())
        SFdict["btaggingSF_particletransformer"][el][year]["up"] = deepcopy(f_eff.Get("btagsfup_particletransformer_{t}_btag_pt_{t}_tw" .format(t = el)).Clone())
        SFdict["btaggingSF_particletransformer"][el][year]["dn"] = deepcopy(f_eff.Get("btagsfdn_particletransformer_{t}_btag_pt_{t}_tw" .format(t = el)).Clone())
f_eff.Close()


def plotSFhisto(tsk):
    thehisto, outfolder, outname, xtitle, ytitle = tsk
    print("\n> Plotting " + thehisto.GetName())

    tdrstyle.setTDRStyle()
    thehisto.SetStats(False)
    thehisto.GetYaxis().SetTitleOffset(1.4)
    thehisto.GetXaxis().SetTitle(xtitle)
    thehisto.GetXaxis().SetTitleOffset(1.2)
    thehisto.GetXaxis().SetTitleFont(43)
    thehisto.GetXaxis().SetTitleSize(22)
    thehisto.GetXaxis().SetLabelFont(43)
    thehisto.GetXaxis().SetLabelSize(22)
    thehisto.GetYaxis().SetTitle(ytitle)
    thehisto.GetYaxis().SetTitleFont(43)
    thehisto.GetYaxis().SetTitleSize(22)
    thehisto.GetYaxis().SetLabelFont(43)
    thehisto.GetYaxis().SetLabelSize(22)

    if "trigger" in outname:
        thehisto.GetXaxis().SetRangeUser(25,200)
        thehisto.GetYaxis().SetRangeUser(25,200)
    
    if "btaggingEff" in outname:
        # We want two separate plots, one for the bin content and one for the error
        c = r.TCanvas('c', "", 600, 600)
        plot = c.GetPad(0)
        plot.SetTopMargin(0.0475); plot.SetRightMargin(0.18); plot.SetLeftMargin(0.12); plot.SetBottomMargin(0.1)
        thehisto.SetMarkerSize(1)
        thehisto.SetMarkerColor(r.kRed)
        thehisto.Draw("colz text45")
        r.gStyle.SetPaintTextFormat("4.3f")
        r.gStyle.SetLabelFont(43, "XYZ")
        r.gStyle.SetLabelSize(22, "XYZ")
        r.gPad.Update()
        c.SaveAs(outfolder + "/" + outname + ".png")
        c.SaveAs(outfolder + "/" + outname + ".pdf")
        c.Close(); del c, plot   
        # Repeat for the error
        c = r.TCanvas('c', "", 600, 600)
        plot = c.GetPad(0)
        plot.SetTopMargin(0.0475); plot.SetRightMargin(0.18); plot.SetLeftMargin(0.12); plot.SetBottomMargin(0.1)
        thehistoerr = thehisto.ProjectionXY()
        thehistoerr.SetMarkerSize(1)
        thehistoerr.SetMarkerColor(r.kRed)
        for iB in range(1, thehisto.GetNbinsX() + 1):
            for jB in range(1, thehisto.GetNbinsY() + 1):
                thehistoerr.SetBinContent(iB, jB, thehisto.GetBinError(iB, jB))
        thehistoerr.Draw("colz text45")
        r.gStyle.SetPaintTextFormat("4.3f")
        r.gStyle.SetLabelFont(43, "XYZ")
        r.gStyle.SetLabelSize(22, "XYZ")
        r.gPad.Update()
        c.SaveAs(outfolder + "/" + outname + "_err.png")
        c.SaveAs(outfolder + "/" + outname + "_err.pdf")
        c.Close(); del c, plot     
    else:
        c = r.TCanvas('c', "", 600, 600)
        plot = c.GetPad(0)
        plot.SetTopMargin(0.0475); plot.SetRightMargin(0.18); plot.SetLeftMargin(0.12); plot.SetBottomMargin(0.1)
        thehisto.SetMarkerSize(1)
        thehisto.SetMarkerColor(r.kRed)
        thehisto.Draw("colz text45 e")
        # Create the rootfile with the histos
        #rootFile = r.TFile(outfolder + "/" + outname + ".root", "RECREATE")
        #thehisto.Write()
        #rootFile.Close()
        r.gStyle.SetPaintTextFormat("4.3f")
        r.gStyle.SetLabelFont(43, "XYZ")
        r.gStyle.SetLabelSize(22, "XYZ")
        r.gPad.Update()
        c.SaveAs(outfolder + "/" + outname + ".png")
        c.SaveAs(outfolder + "/" + outname + ".pdf")
        c.Close(); del c, plot
    return


def plotBtagSFhisto(tsk):
    thedict, outfolder, outname, xtitle, ytitle, year = tsk
    print("\n> Plotting btagging SF")

    theprofile   = thedict["B"][year][""]
    theprofileup = thedict["B"][year]["up"]
    theprofiledn = thedict["B"][year]["dn"]
    # project the TH2F into a TH1F
    thehisto = theprofile.ProjectionX()
    thehistoup = theprofileup.ProjectionX()
    thehistodn = theprofiledn.ProjectionX()
    thearr     = []
    
    for iB in range(1, thehisto.GetNbinsX() + 2):
        thearr.append(thehisto.GetXaxis().GetBinLowEdge(iB))
    thearr = array("f", thearr)
    
    horlin     = r.TH1F("horlin", "", thehisto.GetNbinsX(), thearr)
    tdrstyle.setTDRStyle()
#    thehisto.SetStats(False)
    thehisto.GetXaxis().SetTitle(xtitle)
    thehisto.GetXaxis().SetTitleOffset(1.2)
    thehisto.GetXaxis().SetTitleFont(43)
    thehisto.GetXaxis().SetTitleSize(22)
    thehisto.GetXaxis().SetLabelFont(43)
    thehisto.GetXaxis().SetLabelSize(22)
    thehisto.GetYaxis().SetTitleOffset(1.5)
    thehisto.GetYaxis().SetTitle(ytitle)
    thehisto.GetYaxis().SetTitleFont(43)
    thehisto.GetYaxis().SetTitleSize(22)
    thehisto.GetYaxis().SetLabelFont(43)
    thehisto.GetYaxis().SetLabelSize(22)
    thehisto.SetLineColor(r.kRed)
    thehisto.SetLineWidth(3)
    thehisto.SetMarkerSize(0.5)
    horlin.SetLineColor(r.kBlack)
    """
    xtemp   = r.Double(0.)
    ytemp   = r.Double(0.)
    xtempup = r.Double(0.)
    ytempup = r.Double(0.)
    xtempdn = r.Double(0.)
    ytempdn = r.Double(0.)
    for iB in range(thehisto.GetN()):
        thehisto.GetPoint(iB,   xtemp,   ytemp)
        thehistoup.GetPoint(iB, xtempup, ytempup)
        thehistodn.GetPoint(iB, xtempdn, ytempdn)
        print xtemp, ytemp, xtempup, ytempup, xtempdn, ytempdn
        thehisto.SetPointEYhigh(iB, ytempup if ytempup > ytemp else ytempdn)
        thehisto.SetPointEYlow( iB, ytempdn if ytempup > ytemp else ytempup)
    """
    
    for iB in range(1, thehisto.GetNbinsX() + 1):
        thehisto.SetBinError(iB, (abs(thehisto.GetBinContent(iB) - thehistoup.GetBinContent(iB)) + abs(thehisto.GetBinContent(iB) - thehistodn.GetBinContent(iB)))/2.)
        horlin.SetBinContent(iB, 1.); horlin.SetBinError(iB, 0.);
    c = r.TCanvas('c', "", 600, 600)
    plot = c.GetPad(0);
    plot.SetTopMargin(0.0475); plot.SetRightMargin(0.15); plot.SetLeftMargin(0.12); plot.SetBottomMargin(0.1)
    thehisto.GetYaxis().SetRangeUser(0.8, 1.8)
#    thehisto.Draw("A,P0")
    thehisto.Draw("P0")
    horlin.Draw("E2,same")
    
    tmpCprofile   = thedict["C"][year][""]
    tmpCprofileup = thedict["C"][year]["up"]
    tmpCprofiledn = thedict["C"][year]["dn"]
    tmpC     = tmpCprofile.ProjectionX()
    tmpCup   = tmpCprofileup.ProjectionX()
    tmpCdn   = tmpCprofiledn.ProjectionX()
    tmpC.SetLineColor(r.kGreen + 3)
    tmpC.SetMarkerSize(0.5)
    tmpC.SetLineWidth(2)
    """
    for iB in range(tmpC.GetN()):
        tmpC.GetPoint(iB,   xtemp,   ytemp)
        tmpCup.GetPoint(iB, xtempup, ytempup)
        tmpCdn.GetPoint(iB, xtempdn, ytempdn)
        tmpC.SetPointEYhigh(iB, ytempup if ytempup > ytemp else ytempdn)
        tmpC.SetPointEYlow( iB, ytempdn if ytempup > ytemp else ytempup)"""
    for iB in range(1, tmpC.GetNbinsX() + 1):
        tmpC.SetBinError(iB, (abs(tmpC.GetBinContent(iB) - tmpCup.GetBinContent(iB)) + abs(tmpC.GetBinContent(iB) - tmpCdn.GetBinContent(iB)))/2.)
    tmpC.Draw("P0,same")
    
    tmpLprofile   = thedict["L"][year][""]
    tmpLprofileup = thedict["L"][year]["up"]
    tmpLprofiledn = thedict["L"][year]["dn"]
    tmpL     = tmpLprofile.ProjectionX()
    tmpLup   = tmpLprofileup.ProjectionX()
    tmpLdn   = tmpLprofiledn.ProjectionX()
    tmpL.SetLineColor(r.kBlue)
    tmpL.SetMarkerSize(0.5)
    """
    for iB in range(tmpL.GetN()):
        tmpL.GetPoint(iB,   xtemp,   ytemp)
        tmpLup.GetPoint(iB, xtempup, ytempup)
        tmpLdn.GetPoint(iB, xtempdn, ytempdn)
        tmpL.SetPointEYhigh(iB, ytempup if ytempup > ytemp else ytempdn)
        tmpL.SetPointEYlow( iB, ytempdn if ytempup > ytemp else ytempup)"""
    for iB in range(1, tmpL.GetNbinsX() + 1):
        tmpL.SetBinError(iB, (abs(tmpL.GetBinContent(iB) - tmpLup.GetBinContent(iB)) + abs(tmpL.GetBinContent(iB) - tmpLdn.GetBinContent(iB)))/2.)
    tmpL.Draw("P0,same")

    r.gStyle.SetLabelFont(43, "XYZ")
    r.gStyle.SetLabelSize(22, "XYZ")
    r.gPad.Update()
    
    leg = r.TLegend(0.5, 0.2, 0.7, 0.3)
    leg.SetTextFont(43)
    leg.SetTextSize(12)
    leg.SetBorderSize(0)
    leg.SetFillColor(10)
    leg.SetFillStyle(0); # transparent legend!
    leg.AddEntry(thehisto, "b jets", "L")
    leg.AddEntry(tmpC,     "c jets", "L")
    leg.AddEntry(tmpL,     "light & gluon jets", "L")
    leg.Draw("same")
    c.SaveAs(outfolder + "/" + outname + ".png")
    c.SaveAs(outfolder + "/" + outname + ".pdf")
    c.Close(); del c, plot
    return


axisdict = {}
axisdict["lepton"]   = {}
axisdict["trigger"]  = {}
axisdict["btagging"] = {}


axisdict["lepton"]["m"] = {"x" : "Muon |#eta| (adim.)",
                           "y" : "Muon p_{T} (GeV)"}
axisdict["lepton"]["e"] = {"x" : "Electron supercluster #eta (adim.)",
                           "y" : "Electron p_{T} (GeV)"}
axisdict["trigger"]     = {"x" : "Electron p_{T} (GeV)",
                           "y" : "Muon p_{T} (GeV)"}
axisdict["btagging"]    = {"x" : "Jet p_{T} (GeV)",
                           "y" : "Jet |#eta| (adim.)",
                           "ysf":"Scale factor (adim.)"}

def getTasks(outdir):
    if not os.path.isdir(outdir):
        os.system("mkdir -p " + outdir)

    alltasks = []
    for y in ["2022", "2022PostEE"]:
        if not os.path.isdir(outdir + "/" + str(y)):
            os.system("mkdir -p " + outdir + "/" + str(y))
        # Lepton
        for ty in ["e", "m"]:
            for iH,theH in SFdict["lepton"][ty][y].items():
                if any([subel in iH for subel in ["stat", "syst"]]): continue
                alltasks.append( (theH,
                                  outdir + "/" + str(y),
                                  "leptonSF_" + ty + "_" + iH.replace(",", "_") + "_" + str(y),
                                   axisdict["lepton"][ty]["x"],
                                   axisdict["lepton"][ty]["y"]) )

        # Trigger
        for iC in [ch.ElMu, ch.Muon, ch.Elec]:
            alltasks.append( (SFdict["trigger"][iC][y],
                              outdir + "/" + str(y),
                              "triggerSF_" + str(iC).replace("ch.", "") + "_" + str(y),
                              axisdict["trigger"]["x"],
                              axisdict["trigger"]["y"]
                              ) )

        # Btagging
        for el in ["B", "C", "L"]:
            alltasks.append( (SFdict["btagging_deepjet"][el][y],
                              outdir + "/" + str(y),
                              "btaggingEff_deepjet_" + el + "_" + str(y),
                              axisdict["btagging"]["x"],
                              axisdict["btagging"]["y"]) )
            alltasks.append( (SFdict["btagging_particlenet"][el][y],
                              outdir + "/" + str(y),
                              "btaggingEff_particlenet_" + el + "_" + str(y),
                              axisdict["btagging"]["x"],
                              axisdict["btagging"]["y"]) )
            alltasks.append( (SFdict["btagging_particletransformer"][el][y],
                              outdir + "/" + str(y),
                              "btaggingEff_particletransformer_" + el + "_" + str(y),
                              axisdict["btagging"]["x"],
                              axisdict["btagging"]["y"]) )
    return alltasks


if __name__=="__main__":
    parser = argparse.ArgumentParser(usage = "python plotUncsVariations.py folder", description = "Helper for plotting.", formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outfolder')
    parser.add_argument("-j", metavar = 'ncores',  dest = "ncores", required = False, default = 1, type = int)

    args     = parser.parse_args()
    outfolder= args.outfolder
    ncores   = args.ncores

    tasks    = getTasks(outfolder)

    #print(tasks)
    #sys.exit()

    if ncores > 1:
        pool = Pool(ncores)
        pool.map(plotSFhisto, tasks)
        pool.close()
        pool.join()
    else:
        for tsk in tasks:
            plotSFhisto(tsk)
    
    btagsftasks = []
    for iY in ["2022", "2022PostEE"]:
        btagsftasks.append( (SFdict["btaggingSF_deepjet"],
                             outfolder + "/" + str(iY),
                             "btaggingSF_deepjet_" + str(iY),
                             axisdict["btagging"]["x"],
                             axisdict["btagging"]["ysf"],
                             iY) )
        btagsftasks.append( (SFdict["btaggingSF_particlenet"],
                             outfolder + "/" + str(iY),
                             "btaggingSF_particlenet_" + str(iY),
                             axisdict["btagging"]["x"],
                             axisdict["btagging"]["ysf"],
                             iY) )
        btagsftasks.append( (SFdict["btaggingSF_particletransformer"],
                             outfolder + "/" + str(iY),
                             "btaggingSF_particletransformer_" + str(iY),
                             axisdict["btagging"]["x"],
                             axisdict["btagging"]["ysf"],
                             iY) )
    
    if ncores > 1:
        pool = Pool(ncores)
        pool.map(plotBtagSFhisto, btagsftasks)
        pool.close()
        pool.join()
    else:
        for tsk in btagsftasks:
            plotBtagSFhisto(tsk)
