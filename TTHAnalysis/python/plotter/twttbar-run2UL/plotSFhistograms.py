import os, sys, argparse, enum
from multiprocessing import Pool
from copy import deepcopy
import warnings as wr
import ROOT as r

sys.path.append('{cmsswpath}/src/CMGTools/TTHAnalysis/python/plotter/tw-run2/differential/'.format(cmsswpath = os.environ['CMSSW_BASE']))
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

def loadHistoWithUncs(fil, hist, uncs = []):
    tf = r.TFile.Open(fil)
    if not tf: raise RuntimeError("[lepScaleFactors_TopRun2::loadHistoWithVars] FATAL: no such file %s"%fil)
    histnom = tf.Get(hist)
    if not hist: raise RuntimeError("[lepScaleFactors_TopRun2::loadHistoWithVars] FATAL: no such object %s in %s"%(hist,fil))
    ret = [deepcopy(histnom)]

    for var in uncs:
        tmphist = tf.Get(hist + var)
        ret.append(deepcopy(tmphist))
    tf.Close()
    return ret

sfpath       = os.environ['CMSSW_BASE'] + "/src/CMGTools/TTHAnalysis/data/TopRun2UL/"
basepathlep  = sfpath + "/lepton/"
basepathtrig = sfpath + "/trigger/"
basepathbtag = sfpath + "/btagging/"
basepathjetpuid = sfpath + "/jetPUid/"
SFdict = {}
SFdict["lepton"] = {}; SFdict["lepton"] = {}; SFdict["trigger"] = {}
SFdict["lepton"]["m"] = {}; SFdict["lepton"]["e"] = {}
SFdict["trigger"][ch.ElMu] = {}; SFdict["trigger"][ch.Elec] = {}; SFdict["trigger"][ch.Muon] = {};
SFdict["jetpuid"] = {}
for y in ["2016", "2016apv", "2017", "2018"]:
    SFdict["lepton"]["m"][y] = {}
    SFdict["lepton"]["e"][y] = {}
    SFdict["jetpuid"][y] = {}
    for chan in [ch.ElMu, ch.Elec, ch.Muon]:
        SFdict["trigger"][chan][y] = {}

## Muon ID
SFdict["lepton"]["m"]["2016apv"]["idtight,BCDEF"], SFdict["lepton"]["m"]["2016apv"]["idtight_stat,BCDEF"], SFdict["lepton"]["m"]["2016apv"]["idtight_syst,BCDEF"] = loadHistoWithUncs(basepathlep + "Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root", "NUM_TightID_DEN_TrackerMuons_abseta_pt"  , ["_stat", "_syst"])
SFdict["lepton"]["m"]["2016"]["idtight,GH"], SFdict["lepton"]["m"]["2016"]["idtight_stat,GH"], SFdict["lepton"]["m"]["2016"]["idtight_syst,GH"] = loadHistoWithUncs(basepathlep + "Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root",    "NUM_TightID_DEN_TrackerMuons_abseta_pt"          , ["_stat", "_syst"])
SFdict["lepton"]["m"]["2017"]["idtight"], SFdict["lepton"]["m"]["2017"]["idtight_stat"], SFdict["lepton"]["m"]["2017"]["idtight_syst"] = loadHistoWithUncs(basepathlep + "Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root", "NUM_TightID_DEN_TrackerMuons_abseta_pt"       , ["_stat", "_syst"])
SFdict["lepton"]["m"]["2018"]["idtight"], SFdict["lepton"]["m"]["2018"]["idtight_stat"], SFdict["lepton"]["m"]["2018"]["idtight_syst"] = loadHistoWithUncs(basepathlep + "Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root",  "NUM_TightID_DEN_TrackerMuons_abseta_pt"    , ["_stat", "_syst"])

# Muon iso
SFdict["lepton"]["m"]["2016apv"]["iso,BCDEF"], SFdict["lepton"]["m"]["2016apv"]["iso_stat,BCDEF"], SFdict["lepton"]["m"]["2016apv"]["iso_syst,BCDEF"] = loadHistoWithUncs(basepathlep + "Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt"   , ["_stat", "_syst"])
SFdict["lepton"]["m"]["2016"]["iso,GH"], SFdict["lepton"]["m"]["2016"]["iso_stat,GH"], SFdict["lepton"]["m"]["2016"]["iso_syst,GH"] = loadHistoWithUncs(basepathlep + "Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root",    "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt"   , ["_stat", "_syst"])
SFdict["lepton"]["m"]["2017"]["iso"], SFdict["lepton"]["m"]["2017"]["iso_stat"], SFdict["lepton"]["m"]["2017"]["iso_syst"] = loadHistoWithUncs(basepathlep + "Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt", ["_stat", "_syst"])
SFdict["lepton"]["m"]["2018"]["iso"], SFdict["lepton"]["m"]["2018"]["iso_stat"], SFdict["lepton"]["m"]["2018"]["iso_syst"] = loadHistoWithUncs(basepathlep + "Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root",  "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt", ["_stat", "_syst"])

# Elec ID
SFdict["lepton"]["e"]["2016apv"]["idtight"] = loadHisto(basepathlep + "Electron_2016apvUL_IDTight.root", "EGamma_SF2D")
SFdict["lepton"]["e"]["2016"]["idtight"] = loadHisto(basepathlep + "Electron_2016UL_IDTight.root", "EGamma_SF2D")
SFdict["lepton"]["e"]["2017"]["idtight"] = loadHisto(basepathlep + "Electron_2017UL_IDTight.root", "EGamma_SF2D")
SFdict["lepton"]["e"]["2018"]["idtight"] = loadHisto(basepathlep + "Electron_2018UL_IDTight.root", "EGamma_SF2D")

# Elec reco
SFdict["lepton"]["e"]["2016apv"]["recotight"] = loadHisto(basepathlep + "Electron_2016apvUL_RECO.root", "EGamma_SF2D")
SFdict["lepton"]["e"]["2016"]["recotight"] = loadHisto(basepathlep + "Electron_2016UL_RECO.root", "EGamma_SF2D")
SFdict["lepton"]["e"]["2017"]["recotight"] = loadHisto(basepathlep + "Electron_2017UL_RECO.root", "EGamma_SF2D")
SFdict["lepton"]["e"]["2018"]["recotight"] = loadHisto(basepathlep + "Electron_2018UL_RECO.root", "EGamma_SF2D")

# Trigger elmu
SFdict["trigger"][ch.ElMu]["2016apv"] = loadHisto(basepathtrig + "TriggerSF_2016preVFP_ULv2.root", "h2D_SF_emu_lepABpt_FullError")
SFdict["trigger"][ch.ElMu]["2016"] = loadHisto(basepathtrig + "TriggerSF_2016postVFP_ULv2.root", "h2D_SF_emu_lepABpt_FullError")
SFdict["trigger"][ch.ElMu]["2017"] = loadHisto(basepathtrig + "TriggerSF_2017_ULv2.root", "h2D_SF_emu_lepABpt_FullError")
SFdict["trigger"][ch.ElMu]["2018"] = loadHisto(basepathtrig + "TriggerSF_2018_ULv2.root", "h2D_SF_emu_lepABpt_FullError")

# Trigger elel
SFdict["trigger"][ch.Elec]["2016apv"] = loadHisto(basepathtrig + "TriggerSF_2016preVFP_ULv2.root", "h2D_SF_ee_lepABpt_FullError")
SFdict["trigger"][ch.Elec]["2016"] = loadHisto(basepathtrig + "TriggerSF_2016postVFP_ULv2.root", "h2D_SF_ee_lepABpt_FullError")
SFdict["trigger"][ch.Elec]["2017"] = loadHisto(basepathtrig + "TriggerSF_2017_ULv2.root", "h2D_SF_ee_lepABpt_FullError")
SFdict["trigger"][ch.Elec]["2018"] = loadHisto(basepathtrig + "TriggerSF_2018_ULv2.root", "h2D_SF_ee_lepABpt_FullError")

# Trigger mumu
SFdict["trigger"][ch.Muon]["2016apv"] = loadHisto(basepathtrig + "TriggerSF_2016preVFP_ULv2.root", "h2D_SF_mumu_lepABpt_FullError")
SFdict["trigger"][ch.Muon]["2016"] = loadHisto(basepathtrig + "TriggerSF_2016postVFP_ULv2.root", "h2D_SF_mumu_lepABpt_FullError")
SFdict["trigger"][ch.Muon]["2017"] = loadHisto(basepathtrig + "TriggerSF_2017_ULv2.root", "h2D_SF_mumu_lepABpt_FullError")
SFdict["trigger"][ch.Muon]["2018"] = loadHisto(basepathtrig + "TriggerSF_2018_ULv2.root", "h2D_SF_mumu_lepABpt_FullError")

# JetPUId
SFdict["jetpuid"]["2016apv"]["idloose"], SFdict["jetpuid"]["2016apv"]["idloose_syst"] = loadHistoWithUncs(basepathjetpuid + "PUID_106XTraining_ULRun2_EffSFandUncties_v1.root", "h2_eff_sfUL2016APV_L", ["_Systuncty"])
SFdict["jetpuid"]["2016"]["idloose"], SFdict["jetpuid"]["2016"]["idloose_syst"]= loadHistoWithUncs(basepathjetpuid + "PUID_106XTraining_ULRun2_EffSFandUncties_v1.root", "h2_eff_sfUL2016_L", ["_Systuncty"])
SFdict["jetpuid"]["2017"]["idloose"], SFdict["jetpuid"]["2017"]["idloose_syst"]= loadHistoWithUncs(basepathjetpuid + "PUID_106XTraining_ULRun2_EffSFandUncties_v1.root", "h2_eff_sfUL2017_L", ["_Systuncty"])
SFdict["jetpuid"]["2018"]["idloose"], SFdict["jetpuid"]["2018"]["idloose_syst"]= loadHistoWithUncs(basepathjetpuid + "PUID_106XTraining_ULRun2_EffSFandUncties_v1.root", "h2_eff_sfUL2018_L", ["_Systuncty"])

SFdict["btagging_deepjet"] = {}
SFdict["btagging_deepcsv"] = {}
SFdict["btaggingSF_deepjet"] = {}
SFdict["btaggingSF_deepcsv"] = {}
"""f_eff   = r.TFile.Open(basepathbtag + "/BtagMCSF.root", "read")
for el in ["B", "C", "L"]:
    SFdict["btagging"][el] = {}
    for year in [2016, 2017, 2018]:
        SFdict["btagging"][el][year] = deepcopy(f_eff.Get("BtagSF{t}_{alg}{wp}_{y}".format(t   = el,
                                                                                           alg = "DFlav",
                                                                                           wp  = "M",
                                                                                           y   =  year)).Clone())
f_eff.Close()"""

effPath = "/nfs/fanae/user/asoto/Proyectos/tW-Run3/almaLinux/CMSSW_14_0_11/src/CMGTools/TTHAnalysis/python/plotter/temp_Run2_plots/2024_09_08_btagEff_ttbar/"
#f_eff   = r.TFile.Open("temp_2022_01_21_pruebasbtageff/varplots/2016/outeffs.root", "read")
f_eff_2016apv   = r.TFile.Open(effPath + "/2016apv/eff/output.root", "read")
f_eff_2016   = r.TFile.Open(effPath + "/2016/eff/output.root", "read")
f_eff_2017   = r.TFile.Open(effPath + "/2017/eff/output.root", "read")
f_eff_2018   = r.TFile.Open(effPath + "/2018/eff/output.root", "read")
f_eff_dict = {"2016apv": f_eff_2016apv, "2016": f_eff_2016, "2017": f_eff_2017, "2018": f_eff_2018}
for el in ["B", "C", "L"]:
    SFdict["btagging_deepjet"][el] = {}
    SFdict["btagging_deepcsv"][el] = {}
    SFdict["btaggingSF_deepjet"][el] = {}
    SFdict["btaggingSF_deepcsv"][el] = {}
    for year in ["2016apv", "2016", "2017", "2018"]:
        f_eff = f_eff_dict[year]
        SFdict["btagging_deepjet"][el][year]   = deepcopy(f_eff.Get("btageff_deepjet_{t}_btag_{t}_ttbar".format(t = el)).Clone())
        #SFdict["btagging_deepcsv"][el][year]   = deepcopy(f_eff.Get("btageff_deepcsv_{t}_btag_{t}_ttbar".format(t = el)).Clone())
        SFdict["btaggingSF_deepjet"][el][year] = {}
        #SFdict["btaggingSF_deepcsv"][el][year] = {}
        SFdict["btaggingSF_deepjet"][el][year][""]   = deepcopy(f_eff.Get("btagsf_deepjet_{t}_btag_pt_{t}_ttbar" .format(t = el)).Clone())
        #SFdict["btaggingSF_deepcsv"][el][year][""]   = deepcopy(f_eff.Get("btagsf_deepcsv_{t}_btag_pt_{t}_ttbar" .format(t = el)).Clone())
        SFdict["btaggingSF_deepjet"][el][year]["up"] = deepcopy(f_eff.Get("btagsfup_deepjet_{t}_btag_pt_{t}_ttbar" .format(t = el)).Clone())
        #SFdict["btaggingSF_deepcsv"][el][year]["up"] = deepcopy(f_eff.Get("btagsfup_deepcsv_{t}_btag_pt_{t}_ttbar" .format(t = el)).Clone())
        SFdict["btaggingSF_deepjet"][el][year]["dn"] = deepcopy(f_eff.Get("btagsfdn_deepjet_{t}_btag_pt_{t}_ttbar" .format(t = el)).Clone())
        #SFdict["btaggingSF_deepcsv"][el][year]["dn"] = deepcopy(f_eff.Get("btagsfdn_deepcsv_{t}_btag_pt_{t}_ttbar" .format(t = el)).Clone())
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
        thehisto.GetZaxis().SetTitle("")
        if "jetpuid" in outname:
            thehisto.GetZaxis().SetRangeUser(0.95,1.05)
            thehisto.Draw("colz text45")
            palette = thehisto.GetListOfFunctions().FindObject("palette")
            # If the palette exists, move and customize it
            if palette:
                # Move the palette to a new position
                palette.SetX1NDC(0.85)  # New X1 in normalized device coordinates (NDC)
                palette.SetX2NDC(0.9)   # New X2 in NDC
                palette.SetY1NDC(0.1)   # New Y1 in NDC
                palette.SetY2NDC(0.95)   # New Y2 in NDC
        else:
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
    thehisto.GetYaxis().SetRangeUser(0.5, 1.6)
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
axisdict["jetpuid"]  = {}
axisdict["btagging"] = {}


axisdict["lepton"]["m"] = {"x" : "Muon |#eta| (adim.)",
                           "y" : "Muon p_{T} (GeV)"}
axisdict["lepton"]["e"] = {"x" : "Electron supercluster #eta (adim.)",
                           "y" : "Electron p_{T} (GeV)"}
axisdict["trigger"]     = {"x" : "Electron p_{T} (GeV)",
                           "y" : "Muon p_{T} (GeV)"}
axisdict["jetpuid"]     = {"x" : "Jet p_{T} (GeV)",
                           "y" : "Jet #eta"}
axisdict["btagging"]    = {"x" : "Jet p_{T} (GeV)",
                           "y" : "Jet |#eta| (adim.)",
                           "ysf":"Scale factor (adim.)"}

def getTasks(outdir):
    if not os.path.isdir(outdir):
        os.system("mkdir -p " + outdir)

    alltasks = []
    for y in ["2016apv", "2016", "2017", "2018"]:
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
        # Jet PU id
        for iH,theH in SFdict["jetpuid"][y].items():
            if any([subel in iH for subel in ["syst"]]): continue
            alltasks.append( (theH,
                              outdir + "/" + str(y),
                              "jetpuidSF_" + "_" + str(y),
                               axisdict["jetpuid"]["x"],
                               axisdict["jetpuid"]["y"]) )


        # Btagging
        for el in ["B", "C", "L"]:
            alltasks.append( (SFdict["btagging_deepjet"][el][y],
                              outdir + "/" + str(y),
                              "btaggingEff_deepjet_" + el + "_" + str(y),
                              axisdict["btagging"]["x"],
                              axisdict["btagging"]["y"]) )
#            alltasks.append( (SFdict["btagging_deepcsv"][el][y],
#                              outdir + "/" + str(y),
#                              "btaggingEff_deepcsv_" + el + "_" + str(y),
#                              axisdict["btagging"]["x"],
#                              axisdict["btagging"]["y"]) )
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
    for iY in ["2016apv", "2016", "2017", "2018"]:
        btagsftasks.append( (SFdict["btaggingSF_deepjet"],
                             outfolder + "/" + str(iY),
                             "btaggingSF_deepjet_" + str(iY),
                             axisdict["btagging"]["x"],
                             axisdict["btagging"]["ysf"],
                             iY) )
        #btagsftasks.append( (SFdict["btaggingSF_deepcsv"],
        #                     outfolder + "/" + str(iY),
        #                     "btaggingSF_deepcsv_" + str(iY),
        #                     axisdict["btagging"]["x"],
        #                     axisdict["btagging"]["ysf"],
        #                     iY) )
    
    if ncores > 1:
        pool = Pool(ncores)
        pool.map(plotBtagSFhisto, btagsftasks)
        pool.close()
        pool.join()
    else:
        for tsk in btagsftasks:
            plotBtagSFhisto(tsk)
