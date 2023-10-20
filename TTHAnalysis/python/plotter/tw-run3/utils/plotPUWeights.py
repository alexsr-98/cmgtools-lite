import ROOT as r
import sys, os

path = "/beegfs/data/nanoAODv11/tw-run3/productions/2023-06-02/2022PostEE/"

rootFile = "TTto2L2Nu_TuneCP5_13p6TeV_powheg_pythia8_analysis_0.root"
friendFile = "4_scalefactors/TTto2L2Nu_TuneCP5_13p6TeV_powheg_pythia8_analysis_0_Friend.root"


f = r.TFile.Open(path + rootFile, "READ")
f2 = r.TFile.Open(path + friendFile, "READ")

t = f.Get("Events")
t2 = f2.Get("Friends")

t.AddFriend(t2)

# Now plot puWeight vs PV_npvsGood
c = r.TCanvas("c", "c", 800, 600)
c.cd()

t.Draw("puWeight:PV_npvsGood>>h(100, 0, 100, 100, 0, 5)", "", "colz")

c.SaveAs("test.png")


# Now promediate the puWeight
h = r.TH1F("h", "h", 100, 0, 100)
h.SetStats(False)
h.GetXaxis().SetTitle("Number of good PV")
h.GetYaxis().SetTitle("puWeight")
h.GetXaxis().SetTitleOffset(1.2)
h.GetYaxis().SetTitleOffset(1.2)


