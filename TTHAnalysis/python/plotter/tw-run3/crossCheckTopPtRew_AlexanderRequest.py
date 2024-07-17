import os, sys
import ROOT as r
import cmsstyle as CMS
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

def doSpam(text,x1,y1,x2,y2,align=12,fill=False,textSize=0.033,_noDelete={}):
  cmsprel = r.TPaveText(x1,y1,x2,y2,"NDC");
  cmsprel.SetTextSize(textSize);
  cmsprel.SetFillColor(0);
  cmsprel.SetFillStyle(1001 if fill else 0);
  cmsprel.SetLineStyle(2);
  cmsprel.SetLineColor(0);
  cmsprel.SetLineWidth(0);
  cmsprel.SetTextAlign(align);
  cmsprel.SetTextFont(43);
  cmsprel.AddText(text);
  cmsprel.Draw("same");
  _noDelete[text] = cmsprel; ## so it doesn't get deleted by PyROOT                                                                                                                    
  return cmsprel


# Remove stat box
r.gStyle.SetOptStat(0)

maxEvents = 1000000

# This script tries to reproduce figure from https://twiki.cern.ch/twiki/pub/CMS/TopPtReweighting/Fit.png

year = "2022PostEE"
numberFile = [0]
path = "/lustrefs/hdd_pool_dir/nanoAODv12/tw-run3/checkTopPtRew/TTto2L_13p6/" # https://cmsweb.cern.ch/das/request?instance=prod/global&input=file+dataset%3D%2FTTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8%2FRun3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2%2FNANOAODSIM
rootfile = "296b8455-cd4a-4030-a059-06c11adb857e.root"

year13tev = "2018"
numberFile13tev = [0]
path13tev = "/lustrefs/hdd_pool_dir/nanoAODv12/tw-run3/checkTopPtRew/TTto2L_13/" # https://cmsweb.cern.ch/das/request?instance=prod/global&input=file+dataset%3D%2FTTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8%2FRunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1%2FNANOAODSIM
rootfile13tev = "0804DEBA-97D5-BE46-BB9D-B1125570966E.root"


# Let's define a histogram to study the comparison

hTopPt = r.TH1F("hTopPt", ";p_{T} [GeV];Events", 50, 0, 1000)
hTopPt.SetLineColor(r.kBlue)
hTopPt.SetLineWidth(2)
hTopPt.GetYaxis().SetLabelFont(43)
hTopPt.GetYaxis().SetLabelSize(20)
hTopPt.GetYaxis().SetTitleFont(43)
hTopPt.GetYaxis().SetTitleSize(20)
hTopPt.GetXaxis().SetLabelFont(43)
hTopPt.GetXaxis().SetLabelSize(20)
hTopPt.GetXaxis().SetTitleFont(43)
hTopPt.GetXaxis().SetTitleSize(20)
# offsets y and x axis labels
hTopPt.GetYaxis().SetTitleOffset(1.5)
hTopPt.GetXaxis().SetTitleOffset(3.5)


hTopPt13tev = r.TH1F("hTopPt13tev", ";p_{T} [GeV];Events", 50, 0, 1000)
hTopPt13tev.SetLineColor(r.kRed)
hTopPt13tev.SetLineWidth(2)
hTopPt13tev.GetYaxis().SetLabelFont(43)
hTopPt13tev.GetYaxis().SetLabelSize(20)
hTopPt13tev.GetYaxis().SetTitleFont(43)
hTopPt13tev.GetYaxis().SetTitleSize(20)
hTopPt13tev.GetXaxis().SetLabelFont(43)
hTopPt13tev.GetXaxis().SetLabelSize(20)
hTopPt13tev.GetXaxis().SetTitleFont(43)
hTopPt13tev.GetXaxis().SetTitleSize(20)
hTopPt13tev.GetYaxis().SetTitleOffset(1.5)
hTopPt13tev.GetXaxis().SetTitleOffset(3.5)


# Subleading top
hTopPt2 = r.TH1F("hTopPt2", ";p_{T} [GeV];Events", 50, 0, 1000)
hTopPt2.SetLineColor(r.kBlue)
hTopPt2.SetLineWidth(2)

hTopPt2_13tev = r.TH1F("hTopPt2_13tev", ";p_{T} [GeV];Events", 50, 0, 1000)
hTopPt2_13tev.SetLineColor(r.kRed)
hTopPt2_13tev.SetLineWidth(2)

for iF in numberFile:
    rooFile = r.TFile((path + rootfile).format(year = year, num = iF))
    events = rooFile.Get("Events")
    numberofevents = events.GetEntries()
    for i,event in enumerate(events):
        if i%1000 == 0:
            print("Processing event {0} of {1}".format(i, numberofevents))
        if i > maxEvents:
            break
        partobjs = [p for p in Collection(event, "GenPart")]
        candtops = []
        for i, part in enumerate(partobjs):
            if abs(part.pdgId) == 6 and part.status == 62:
                candtops.append(i)
        if len(candtops) != 2:
            continue
        if partobjs[candtops[0]].pt < partobjs[candtops[1]].pt:
            highptind = 1
            lowptind  = 0
        else:
            highptind = 0
            lowptind  = 1
        hTopPt.Fill(partobjs[candtops[highptind]].pt)
        hTopPt.Fill(partobjs[candtops[lowptind]].pt)
        hTopPt2.Fill(partobjs[candtops[lowptind]].pt)
    rooFile.Close()

############################
############################
############################

for iF in numberFile13tev:
    rooFile = r.TFile((path13tev + rootfile13tev).format(year = year13tev, num = iF))
    events = rooFile.Get("Events")
    numberofevents = events.GetEntries()
    for i,event in enumerate(events):
        if i%1000 == 0:
            print("Processing event {0} of {1}".format(i, numberofevents))
        if i > maxEvents:
            break
        partobjs = [p for p in Collection(event, "GenPart")]
        candtops = []
        for i, part in enumerate(partobjs):
            if abs(part.pdgId) == 6 and part.status == 62:
                candtops.append(i)
        if len(candtops) != 2:
            continue
        if partobjs[candtops[0]].pt < partobjs[candtops[1]].pt:
            highptind = 1
            lowptind  = 0
        else:
            highptind = 0
            lowptind  = 1
        hTopPt13tev.Fill(partobjs[candtops[highptind]].pt)
        hTopPt13tev.Fill(partobjs[candtops[lowptind]].pt)
        hTopPt2_13tev.Fill(partobjs[candtops[lowptind]].pt)
    rooFile.Close()




c = r.TCanvas("c", "c", 800, 600)

# Add legends
leg = r.TLegend(0.6, 0.6, 0.9, 0.9)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.AddEntry(hTopPt, "Top Pt (13.6 TeV)", "l")
leg.AddEntry(hTopPt13tev, "Top Pt (13 TeV)", "l")
# Two pads to draw also the ratio
pad1 = r.TPad("pad1", "pad1", 0, 0.3, 1, 1)
pad1.SetBottomMargin(0.02)
pad1.Draw()
pad2 = r.TPad("pad2", "pad2", 0, 0, 1, 0.3)
pad2.SetTopMargin(0.1)
pad2.SetBottomMargin(0.3)
pad2.Draw()
pad1.cd()
hTopPt13tev.GetXaxis().SetLabelSize(0)
hTopPt13tev.Draw("axis")
hTopPt.Draw("samehistE")
hTopPt13tev.Draw("histEsame")
leg.Draw("same")


pad2.cd()
hTopPtRatio = hTopPt.Clone("hTopPtRatio")
hTopPtRatio.Divide(hTopPt13tev)
# Propagate the uncertainty by hand to the ratio
for i in range(1, hTopPtRatio.GetNbinsX()+1):
    try:
        # Bad
        #hTopPtRatio.SetBinError(i, hTopPt.GetBinError(i)/hTopPt13tev.GetBinContent(i))
        # Good (using error propagation in denominator and numerator)
        hTopPtRatio.SetBinError(i, hTopPtRatio.GetBinContent(i)*r.TMath.Sqrt((hTopPt.GetBinError(i)/hTopPt.GetBinContent(i))**2 + (hTopPt13tev.GetBinError(i)/hTopPt13tev.GetBinContent(i))**2))

    except ZeroDivisionError:
        hTopPtRatio.SetBinError(i, 0)

hTopPtRatio.SetLineColor(r.kBlack)
hTopPtRatio.SetMarkerStyle(20)
hTopPtRatio.SetMarkerSize(0.8)
hTopPtRatio.SetMarkerColor(r.kBlack)
hTopPtRatio.SetMinimum(0.8)
hTopPtRatio.SetMaximum(1.2)
hTopPtRatio.GetYaxis().SetNdivisions(505)

hTopPtRatio.GetYaxis().SetTitle("Ratio")
# Center the title
hTopPtRatio.GetYaxis().CenterTitle()

# Do a linear fit to the ratio and draw it
fit = r.TF1("fit", "pol1", 0, 500)
hTopPtRatio.Fit(fit, "R")
fit.SetLineColor(r.kRed)

# Extend the line to the full range
hTopPtRatio.Draw("ep")
fit.Draw("same")
# Draw another line using parameters from the fit and extending it to the full range
fit2 = r.TF1("fit2", "pol1", 0, 1000)
fit2.SetParameters(fit.GetParameter(0), fit.GetParameter(1))
fit2.SetLineColor(r.kRed)

fit2.Draw("same")

# Add a legend for the fit with the fit parameters expressed as p1 * p_{T} + p0 (with p0 and p1 the fit parameters)
legFitTopPtRatio = r.TLegend(0.1, 0.6, 0.4, 0.9)
legFitTopPtRatio.SetBorderSize(0)
legFitTopPtRatio.SetFillStyle(0)
legFitTopPtRatio.AddEntry(fit2, "Fit: {0:.3f} + {1:.6f} ".format(fit.GetParameter(0), fit.GetParameter(1)) + "pt", "l")
legFitTopPtRatio.Draw("same")






# Draw horizontal and vertical grid lines using setgrid
pad2.SetGrid()

pad1.cd()
doSpam('#scale[1.2]{#bf{CMS}}#scale[1.0]{#it{Simulation Preliminary}}',.1, .915, .35, .955,textSize = 22)

pad2.cd()


#c.Draw()
c.SaveAs("TopPt.png")
c.SaveAs("TopPt.pdf")
c.SaveAs("TopPt.root")



c2 = r.TCanvas("c2", "c2", 800, 600)
hTopPt2.Draw()
hTopPt2_13tev.Draw("same")
c2.Draw()
c2.SaveAs("TopPt2.png")