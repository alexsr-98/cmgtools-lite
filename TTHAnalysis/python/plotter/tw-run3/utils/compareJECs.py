import correctionlib._core as core
from copy import deepcopy
import ROOT
import os
from array import array

path = "/mnt_pool/c3_users/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/utils/JECjsons_forComparison"
jsonWinter = "jet_jerc_Winter22.json.gz"
jsonSummer = "jet_jerc_Summer22.json.gz"

# Load the JECs

jerc_corrs_winter = core.CorrectionSet.from_file(os.path.join(path, jsonWinter))
jerc_corrs_summer = core.CorrectionSet.from_file(os.path.join(path, jsonSummer))


jer_corrector_winter = jerc_corrs_winter["JR_Winter22Run3_V1_MC_ScaleFactor_AK4PFPuppi"]
jer_corrector_summer = jerc_corrs_summer["Summer22EEPrompt22_JRV1_MC_ScaleFactor_AK4PFPuppi"]

print([input.name for input in jer_corrector_winter.inputs])
print([input.name for input in jer_corrector_summer.inputs])

# Do a scan in eta from -2.4 and 2.4, evaluate the SFs and plot them

eta = -2.4
etas = []
sf_winter = []
sf_winter_up = []
sf_winter_down = []
sf_summer = []
sf_summer_up = []
sf_summer_down = []
while eta < 2.4:
    etas.append(eta)
    sf_winter.append(jer_corrector_winter.evaluate(eta,30.,'nom'))
    sf_winter_up.append(jer_corrector_winter.evaluate(eta,30.,'up'))
    sf_winter_down.append(jer_corrector_winter.evaluate(eta,30.,'down'))
    sf_summer.append(jer_corrector_summer.evaluate(eta, 'nom'))
    sf_summer_up.append(jer_corrector_summer.evaluate(eta, 'up'))
    sf_summer_down.append(jer_corrector_summer.evaluate(eta, 'down'))
    eta += 0.01

# Plot the SFs

c = ROOT.TCanvas("c", "c", 800, 600)
c.SetGrid()
c.SetLeftMargin(0.15)
c.SetBottomMargin(0.15)

gr_winter = ROOT.TGraph(len(etas), array("d", etas), array("d", sf_winter))
gr_winter_up = ROOT.TGraph(len(etas), array("d", etas), array("d", sf_winter_up))
gr_winter_down = ROOT.TGraph(len(etas), array("d", etas), array("d", sf_winter_down))
gr_summer = ROOT.TGraph(len(etas), array("d", etas), array("d", sf_summer))
gr_summer_up = ROOT.TGraph(len(etas), array("d", etas), array("d", sf_summer_up))
gr_summer_down = ROOT.TGraph(len(etas), array("d", etas), array("d", sf_summer_down))

gr_winter.SetLineColor(ROOT.kRed)
gr_winter_up.SetLineColor(ROOT.kRed)
gr_winter_down.SetLineColor(ROOT.kRed)
gr_summer.SetLineColor(ROOT.kBlue)
gr_summer_up.SetLineColor(ROOT.kBlue)
gr_summer_down.SetLineColor(ROOT.kBlue)

gr_winter.SetLineWidth(2)
gr_winter_up.SetLineWidth(2)
gr_winter_down.SetLineWidth(2)
gr_summer.SetLineWidth(2)
gr_summer_up.SetLineWidth(2)
gr_summer_down.SetLineWidth(2)

gr_winter.SetLineStyle(1)
gr_winter_up.SetLineStyle(1)
gr_winter_down.SetLineStyle(1)
gr_summer.SetLineStyle(1)
gr_summer_up.SetLineStyle(1)
gr_summer_down.SetLineStyle(1)

# Make less transparent up and down
gr_winter_up.SetLineColorAlpha(ROOT.kRed, 0.2)
gr_winter_down.SetLineColorAlpha(ROOT.kRed, 0.2)
gr_summer_up.SetLineColorAlpha(ROOT.kBlue, 0.2)
gr_summer_down.SetLineColorAlpha(ROOT.kBlue, 0.2)

gr_winter.SetTitle("JER SFs")
gr_winter.GetXaxis().SetTitle("#eta")
gr_winter.GetYaxis().SetTitle("SF")

gr_winter.Draw("AL")
gr_winter_up.Draw("L SAME")
gr_winter_down.Draw("L SAME")
gr_summer.Draw("L SAME")
gr_summer_up.Draw("L SAME")
gr_summer_down.Draw("L SAME")

leg = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
leg.AddEntry(gr_winter, "Winter22", "L")
leg.AddEntry(gr_summer, "Summer22", "L")
leg.Draw()

c.SaveAs(path + "/jerSFs.pdf")
c.SaveAs(path + "/jerSFs.png")
