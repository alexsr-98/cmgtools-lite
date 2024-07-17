import ROOT as r
import os
import sys
import argparse
import math
from array import array

r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(0)
# Axis ticks in all sides of the plot
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)

# limit the decimals in SetPaintTextFormat
r.gStyle.SetPaintTextFormat(".3f")

'''
Theory: https://indico.cern.ch/event/731987/contributions/3252702/attachments/1785783/2908948/tutorial_trigeff.pdf
'''

#Previous, v4 of the note
#lumidict     = {"2022"       : 7.78,
#                "2022PostEE" : 20.67}
lumidict     = {"2022"       : 7.88,
                "2022PostEE" : 26.34}

# Year
spamsSize = 30

# Some dictionaries
xAxisTitles = { "trigger" : ["Electron p_{T} (GeV)","Electron p_{T} (GeV)","Muon p_{T} (GeV)","Muon #eta"],
                "trigger-ee" : ["Leading electron p_{T} (GeV)"],
                "trigger-mm" : ["Leading muon p_{T} (GeV)"]}

yAxisTitles = { "trigger" : ["Muon p_{T} (GeV)","Electron #eta","Muon #eta","Electron #eta"],
                "trigger-ee" : ["Subleading electron p_{T} (GeV)"],
                "trigger-mm" : ["Subleading muon p_{T} (GeV)"]}

savingNames = { "trigger" : ["h2D_SF_emu_lepABpt_FullError","h2D_SF_emu_lepABeta1_FullError","h2D_SF_emu_lepABeta2_FullError","h2D_SF_emu_lepABeta3_FullError"],
                "trigger-ee" : ["h2D_SF_ee_lepABpt_FullError"],
                "trigger-mm" : ["h2D_SF_mumu_lepABpt_FullError"]}

histoNames = { "trigger" : ["ept_mpt","ept_eeta","mpt_meta","meta_eeta"],
                "trigger-ee" : ["lep1pt_lep2pt"],
                "trigger-mm" : ["lep1pt_lep2pt"]}

# Dictionaries for 1D histograms
xAxisTitles1D = { "trigger" : ["Electron p_{T} (GeV)","Muon p_{T} (GeV)", "Electron #eta", "Muon #eta"],
                "trigger-ee" : ["Leading electron p_{T} (GeV)", "Subleading electron p_{T} (GeV)", "Leading electron #eta", "Subleading electron #eta"],
                "trigger-mm" : ["Leading muon p_{T} (GeV)", "Subleading muon p_{T} (GeV)", "Leading muon #eta", "Subleading muon #eta"]}

yAxisTitles1D = { "trigger" : ["Efficiency","Efficiency", "Efficiency", "Efficiency"],
                "trigger-ee" : ["Efficiency", "Efficiency", "Efficiency", "Efficiency"],
                "trigger-mm" : ["Efficiency", "Efficiency", "Efficiency", "Efficiency"]}

histoNames1D = { "trigger" : ["e_pt","m_pt", "e_eta", "m_eta"],
                "trigger-ee" : ["lep1_pt", "lep2_pt", "lep1_eta", "lep2_eta"],
                "trigger-mm" : ["lep1_pt", "lep2_pt", "lep1_eta", "lep2_eta"]}

# Dictionaries for the correlation plots
xAxisTitlesCorr = { "trigger" : ["Electron p_{T} (GeV)","Muon p_{T} (GeV)", "Electron #eta", "Muon #eta", ""],
                "trigger-ee" : ["Leading electron p_{T} (GeV)", "Subleading electron p_{T} (GeV)", "Leading electron #eta", "Subleading electron #eta", ""],
                "trigger-mm" : ["Leading muon p_{T} (GeV)", "Subleading muon p_{T} (GeV)", "Leading muon #eta", "Subleading muon #eta", ""]}

yAxisTitlesCorr = { "trigger" : ["Efficiency","Efficiency", "Efficiency", "Efficiency", "Efficiency"],
                "trigger-ee" : ["Efficiency", "Efficiency", "Efficiency", "Efficiency", "Efficiency"],
                "trigger-mm" : ["Efficiency", "Efficiency", "Efficiency", "Efficiency", "Efficiency"]}

histoNamesCorr = { "trigger" : ["e_pt","m_pt", "e_eta", "m_eta", "tot_weight"],
                "trigger-ee" : ["lep1_pt", "lep2_pt", "lep1_eta", "lep2_eta", "tot_weight"],
                "trigger-mm" : ["lep1_pt", "lep2_pt", "lep1_eta", "lep2_eta", "tot_weight"]}

def avoidUnconsistentEff(passed, total):
    '''
    In MC there is the posibility that if the eff is close to 1, some events
    from the pass with negative weights are removed.
    The result is some bins of pass are bigger than the total.

    This function redifines the histograms by taking the value passed = total
    '''
    for i in range(1, passed.GetNbinsX() + 1):
        for j in range(1, passed.GetNbinsY() + 1):
            if passed.GetBinContent(i, j) > total.GetBinContent(i, j):
                passed.SetBinContent(i, j, total.GetBinContent(i, j))

def avoidUnconsistentEff1D(passed, total):
    '''
    In MC there is the posibility that if the eff is close to 1, some events
    from the pass with negative weights are removed.
    The result is some bins of pass are bigger than the total.

    This function redifines the histograms by taking the value passed = total
    '''
    for i in range(1, passed.GetNbinsX() + 1):
        if passed.GetBinContent(i) > total.GetBinContent(i):
            passed.SetBinContent(i, total.GetBinContent(i))


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

if __name__=="__main__":
    parser = argparse.ArgumentParser(usage = "python3 plotterHelper.py [options]", description = "Helper for plotting.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inpath',   '-i', metavar = 'inpath',   dest = "inpath", required = False, default = "./temp/varplots")
    parser.add_argument('--channel',    '-c', metavar = 'channel',    dest = "channel",  required = False, default = "mu")
    parser.add_argument('--outpath',  '-o', metavar = 'outpath',  dest = "outpath",  required = False, default = "./temp/varplots")
    parser.add_argument('--year',      '-y', metavar = 'year',      dest = "year",    required = False, default = "2022")


    args     = parser.parse_args()
    inpath  = args.inpath
    channel   = args.channel # Can be: trigger, trigger-ee, trigger-mm
    outpath = args.outpath
    year     = args.year

    # Expand the outpath with the year
    outpath = outpath + "/" + year + "/"
    # Create the output directory if it doesn't exist
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    metOnlyFolder = "METTrig"
    metLepFolder = "METandLeptonTrig"
    noneOnlyFolder = "None"
    leptonOnlyFolder = "LeptonTrig"
    # Plots to study the signal efficiency (ttbar) for the single and double lepton triggers
    singleLeptonNoMETOnlyFolder = "SingleLeptonTrigNoMETcut"
    doubleLeptonNoMETOnlyFolder = "DoubleLeptonTrigNoMETcut"
    leptonNoMETOnlyFolder = "LeptonTrigNoMETcut"
    noneNoMETOnlyFolder = "NoneNoMETcut"
    # For the uncertainties
    metOnlyMore3jFolder = "METTrigMore3j"
    metLepMore3jFolder = "METandLeptonTrigMore3j"
    metOnlyLess3jFolder = "METTrigLess3j"
    metLepLess3jFolder = "METandLeptonTrigLess3j"
    metOnlyMore35vtx = "METTrigMore35vtx"
    metLepMore35vtx = "METandLeptonTrigMore35vtx"
    metOnlyLess35vtx = "METTrigLess35vtx"
    metLepLess35vtx = "METandLeptonTrigLess35vtx"

    plotsRoot = "plots-tw-{channel}.root".format(channel = channel)

    # Get the histograms
    metOnly   = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + metOnlyFolder + "/" + plotsRoot, "READ")
    metLepton = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + metLepFolder + "/" + plotsRoot, "READ")
    noneOnly  = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + noneOnlyFolder + "/" + plotsRoot, "READ")
    leptonOnly = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + leptonOnlyFolder + "/" + plotsRoot, "READ")

    try:
        singleLeptonNoMETOnly = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + singleLeptonNoMETOnlyFolder + "/" + plotsRoot, "READ")
        doubleLeptonNoMETOnly = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + doubleLeptonNoMETOnlyFolder + "/" + plotsRoot, "READ")
        leptonNoMETOnly = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + leptonNoMETOnlyFolder + "/" + plotsRoot, "READ") 
        noneNoMETOnly = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + noneNoMETOnlyFolder + "/" + plotsRoot, "READ")
    except:
        print("Not plotting single and double lepton triggers")

    # For the uncertainties
    metOnlyMore3j   = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + metOnlyMore3jFolder + "/" + plotsRoot, "READ")
    metLepMore3j = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + metLepMore3jFolder + "/" + plotsRoot, "READ")
    metOnlyLess3j   = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + metOnlyLess3jFolder + "/" + plotsRoot, "READ")
    metLepLess3j = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + metLepLess3jFolder + "/" + plotsRoot, "READ")
    metOnlyMore35vtx   = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + metOnlyMore35vtx + "/" + plotsRoot, "READ")
    metLepMore35vtx = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + metLepMore35vtx + "/" + plotsRoot, "READ")
    metOnlyLess35vtx   = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + metOnlyLess35vtx + "/" + plotsRoot, "READ")
    metLepLess35vtx = r.TFile.Open(inpath + "/" + year + "/" + channel + "/" + metLepLess35vtx + "/" + plotsRoot, "READ")

    # List of all folders based on the uncertainties that we want to compute
    uncList = {"stat": [metOnly, metLepton], "more3j" : [metOnlyMore3j, metLepMore3j], "less3j" : [metOnlyLess3j, metLepLess3j], "more35vtx" : [metOnlyMore35vtx, metLepMore35vtx], "less35vtx" : [metOnlyLess35vtx, metLepLess35vtx]}

    #uncList = {"stat": [metOnly, metLepton],"_jesUp" : [metOnly, metLepton], "_jesDown" : [metOnly, metLepton], "_unclenergyUp" : [metOnly, metLepton], "_unclenergyDown" : [metOnly, metLepton]}
    specialUncs = ["_jesUp", "_jesDown", "_unclenergyUp", "_unclenergyDown"]

    for unc in uncList: 
        for idx in range(len(xAxisTitles[channel])):
            if unc in specialUncs:
                uncName = unc
            else:
                uncName = ""
            # Get the 2D histograms data and MC
            histoName = histoNames[channel][idx]
            metOnlyData   = uncList[unc][0].Get("{name}_data".format(name = histoName))
            metOnlyMC     = uncList[unc][0].Get("{name}_ttbar{unc}".format(name = histoName, unc = uncName))

            metLeptonData = uncList[unc][1].Get("{name}_data".format(name = histoName))
            metLeptonMC   = uncList[unc][1].Get("{name}_ttbar{unc}".format(name = histoName, unc = uncName))

            avoidUnconsistentEff(metLeptonMC, metOnlyMC)
            # Efficiency histograms (ratio of metLeptonData to metOnlyData)
            effData = metLeptonData.Clone("effData")
            effData.Divide(metOnlyData)

            effMC = metLeptonMC.Clone("effMC")
            effMC.Divide(metOnlyMC)

            # Compute the scale factors (ratio of effData to effMC)
            sf = effData.Clone("sf")
            sf.Divide(effMC)

            # Plot the scale factors and efficiencies
            c = r.TCanvas("c", "c", 1400, 900)
            # Set margins
            c.SetBottomMargin(0.15)
            c.SetLeftMargin(0.2)
            c.SetRightMargin(0.15)
            c.cd()

            # Scale factors
            sf.SetTitle("Scale Factors")
            sf.GetXaxis().SetTitle(xAxisTitles[channel][idx])
            sf.GetYaxis().SetTitle(yAxisTitles[channel][idx])
            #sf.GetZaxis().SetRangeUser(0.5, 1.)
            sf.Draw("colz TEXT45")
            c.SaveAs(outpath + "/sf_" + channel + str(idx) + unc + year + ".pdf")
            c.SaveAs(outpath + "/sf_" + channel + str(idx) + unc + year + ".png")

            # Efficiencies
            effData.SetTitle("Efficiencies")
            effData.GetXaxis().SetTitle(xAxisTitles[channel][idx])
            effData.GetYaxis().SetTitle(yAxisTitles[channel][idx])
            #effData.GetZaxis().SetRangeUser(0.5, 1.5)
            effData.Draw("colz TEXT45")
            c.SaveAs(outpath + "/effData_" + channel + str(idx) + unc + year + ".pdf")
            c.SaveAs(outpath + "/effData_" + channel + str(idx) + unc + year + ".png")

            effMC.SetTitle("Efficiencies")
            effMC.GetXaxis().SetTitle(xAxisTitles[channel][idx])
            effMC.GetYaxis().SetTitle(yAxisTitles[channel][idx])
            #effMC.GetZaxis().SetRangeUser(0.5, 1.5)
            effMC.Draw("colz TEXT45")
            c.SaveAs(outpath + "/effMC_" + channel + str(idx) + unc + year + ".pdf")
            c.SaveAs(outpath + "/effMC_" + channel + str(idx) + unc + year + ".png")

            # Now plot the uncertainties of the scale factors and efficiencies
            # Get the uncertainties
            # Scale factors
            sfUnc = sf.Clone("sfUnc")
            for i in range(1, sfUnc.GetNbinsX()+1):
                for j in range(1, sfUnc.GetNbinsY()+1):
                    sfUnc.SetBinContent(i, j, sf.GetBinError(i, j))

            # Efficiencies
            effDataUnc = effData.Clone("effDataUnc")
            for i in range(1, effDataUnc.GetNbinsX()+1):
                for j in range(1, effDataUnc.GetNbinsY()+1):
                    effDataUnc.SetBinContent(i, j, effData.GetBinError(i, j))

            effMCUnc = effMC.Clone("effMCUnc")
            for i in range(1, effMCUnc.GetNbinsX()+1):
                for j in range(1, effMCUnc.GetNbinsY()+1):
                    effMCUnc.SetBinContent(i, j, effMC.GetBinError(i, j))

            # Scale factors
            sfUnc.SetTitle("Scale Factor Uncertainties")
            sfUnc.GetXaxis().SetTitle(xAxisTitles[channel][idx])
            sfUnc.GetYaxis().SetTitle(yAxisTitles[channel][idx])
            sfUnc.GetZaxis().SetRangeUser(0., 0.2)
            sfUnc.Draw("colz TEXT45")
            c.SaveAs(outpath + "/sfUnc_" + channel + str(idx) + unc + year + ".pdf")
            c.SaveAs(outpath + "/sfUnc_" + channel + str(idx) + unc + year + ".png")

            # Efficiencies
            effDataUnc.SetTitle("Efficiency Uncertainties")
            effDataUnc.GetXaxis().SetTitle(xAxisTitles[channel][idx])
            effDataUnc.GetYaxis().SetTitle(yAxisTitles[channel][idx])
            effDataUnc.GetZaxis().SetRangeUser(0., 0.2)
            effDataUnc.Draw("colz TEXT45")
            c.SaveAs(outpath + "/effDataUnc_" + channel + str(idx) + unc + year + ".pdf")
            c.SaveAs(outpath + "/effDataUnc_" + channel + str(idx) + unc + year + ".png")

            effMCUnc.SetTitle("Efficiency Uncertainties")
            effMCUnc.GetXaxis().SetTitle(xAxisTitles[channel][idx])
            effMCUnc.GetYaxis().SetTitle(yAxisTitles[channel][idx])
            effMCUnc.GetZaxis().SetRangeUser(0., 0.2)
            effMCUnc.Draw("colz TEXT45")
            c.SaveAs(outpath + "/effMCUnc_" + channel + str(idx) + unc + year + ".pdf")
            c.SaveAs(outpath + "/effMCUnc_" + channel + str(idx) + unc + year + ".png")

            # Now compute the scale factor uncertainties but using the Clopper-Pearson method
            efficiencyDataCP = r.TEfficiency(metLeptonData, metOnlyData)
            efficiencyMCCP = r.TEfficiency(metLeptonMC, metOnlyMC)

            # Create a clone of the eff histograms to store the uncertainties from TEfficiency
            effDataUncCP = effData.Clone("effDataUncCP")
            effMCUncCP = effMC.Clone("effMCUncCP")
            for i in range(1, efficiencyDataCP.GetPassedHistogram().GetNbinsX()+1):
                for j in range(1, efficiencyDataCP.GetPassedHistogram().GetNbinsY()+1):
                    dataUnc = (efficiencyDataCP.GetEfficiencyErrorLow(efficiencyDataCP.GetGlobalBin(i, j)) + efficiencyDataCP.GetEfficiencyErrorUp(efficiencyDataCP.GetGlobalBin(i, j))) / 2.
                    mcUnc = (efficiencyMCCP.GetEfficiencyErrorLow(efficiencyMCCP.GetGlobalBin(i, j)) + efficiencyMCCP.GetEfficiencyErrorUp(efficiencyMCCP.GetGlobalBin(i, j))) / 2.
                    effDataUncCP.SetBinContent(i, j, dataUnc)
                    effMCUncCP.SetBinContent(i, j, mcUnc)

            # Efficiencies
            effDataUncCP.SetTitle("Efficiency Uncertainties (CP)")
            effDataUncCP.GetXaxis().SetTitle(xAxisTitles[channel][idx])
            effDataUncCP.GetYaxis().SetTitle(yAxisTitles[channel][idx])
            effDataUncCP.GetZaxis().SetRangeUser(0., 0.2)
            effDataUncCP.Draw("colz TEXT45")
            c.SaveAs(outpath + "/effDataUncCP_" + channel + str(idx) + unc + year + ".pdf")
            c.SaveAs(outpath + "/effDataUncCP_" + channel + str(idx) + unc + year + ".png")

            effMCUncCP.SetTitle("Efficiency Uncertainties (CP)")
            effMCUncCP.GetXaxis().SetTitle(xAxisTitles[channel][idx])
            effMCUncCP.GetYaxis().SetTitle(yAxisTitles[channel][idx])
            effMCUncCP.GetZaxis().SetRangeUser(0., 0.2)
            effMCUncCP.Draw("colz TEXT45")
            c.SaveAs(outpath + "/effMCUncCP_" + channel + str(idx) + unc + year + ".pdf")
            c.SaveAs(outpath + "/effMCUncCP_" + channel + str(idx) + unc + year + ".png")

            # Propagate the uncertainties to the scale factors, we have efficiencyDataCP, efficiencyMCCP and their uncertainties
            # We apply the formula: SF = effDataCP / effMCCP and propagate the uncertainties

            # Scale factors
            sfUncCP = sf.Clone("sfUncCP")
            for i in range(1, sfUncCP.GetNbinsX()+1):
                for j in range(1, sfUncCP.GetNbinsY()+1):
                    # Check if the denominator is zero
                    if effMC.GetBinContent(i, j) == 0. or effData.GetBinContent(i, j) == 0.:
                        sfUncCP.SetBinContent(i, j, 0.)
                    else:
                        sfUncCP.SetBinContent(i, j, sf.GetBinContent(i, j) * math.sqrt((effDataUncCP.GetBinContent(i, j) / effData.GetBinContent(i, j))**2 + (effMCUncCP.GetBinContent(i, j) / effMC.GetBinContent(i, j))**2))

            # Scale factors
            sfUncCP.SetTitle("Scale Factor Uncertainties (CP)")
            sfUncCP.GetXaxis().SetTitle(xAxisTitles[channel][idx])
            sfUncCP.GetYaxis().SetTitle(yAxisTitles[channel][idx])
            #sfUncCP.GetZaxis().SetRangeUser(0., 0.2)
            sfUncCP.Draw("colz TEXT45")
            c.SaveAs(outpath + "/sfUncCP_" + channel + str(idx) + unc + year + ".pdf")
            c.SaveAs(outpath + "/sfUncCP_" + channel + str(idx) + unc + year + ".png")

            # Output root file with the histos (append if it exists)
            outroot = r.TFile.Open(outpath + "/triggerSFs.root", "UPDATE")
            # We will only write the sf but we have to update the bin errors to the ones from the CP method
            for i in range(1, sf.GetNbinsX()+1):
                for j in range(1, sf.GetNbinsY()+1):
                    sf.SetBinError(i, j, sfUncCP.GetBinContent(i, j))

            # Save the sf with savingNames[channel]
            sf.SetName(savingNames[channel][idx] + unc)
            sf.Write()
            outroot.Close()
    
    # Now collect from the outroot file the scale factors and compute the differences to evaluate the systematics
    outroot = r.TFile.Open(outpath + "/triggerSFs.root", "UPDATE")
    for idx in range(len(xAxisTitles[channel])):
        # Get the 2D histogram
        sfNominal = outroot.Get(savingNames[channel][idx] + "stat")
        # Now, for ease of display, we check if the last bin edge in x and y is 500, if it is, we clone the histo to one with last bin equal to the mean of the last two bin edges
        if sfNominal.GetXaxis().GetBinUpEdge(sfNominal.GetNbinsX()) == 500. and sfNominal.GetYaxis().GetBinUpEdge(sfNominal.GetNbinsY()) == 500.:
            # We create a new TH2D with the last bin equal to the mean of the last two bin edges
            arrayXbins = sfNominal.GetXaxis().GetXbins()
            arrayYbins = sfNominal.GetYaxis().GetXbins()
            # List x bins
            listXbins = []
            for i in range(len(arrayXbins)):
                listXbins.append(arrayXbins[i])
            # List y bins
            listYbins = []
            for i in range(len(arrayYbins)):
                listYbins.append(arrayYbins[i])

            # Replace the last bin edge with the mean of the last two bin edges
            listXbins[-1] = (listXbins[-1] + listXbins[-2]) / 2.
            listYbins[-1] = (listYbins[-1] + listYbins[-2]) / 2.
            # Convert the lists to arrays
            arrayXbins = array('d', listXbins)
            arrayYbins = array('d', listYbins)
            # Give better naming
            sfNominalNew = r.TH2D(savingNames[channel][idx] + "stat", savingNames[channel][idx] + "stat", sfNominal.GetNbinsX(), arrayXbins, sfNominal.GetNbinsY(), arrayYbins)
            # Copy the content of the old histo to the new one, also errors
            for i in range(1, sfNominal.GetNbinsX()+1):
                for j in range(1, sfNominal.GetNbinsY()+1):
                    sfNominalNew.SetBinContent(i, j, sfNominal.GetBinContent(i, j))
                    sfNominalNew.SetBinError(i, j, sfNominal.GetBinError(i, j))
            # Delete the old histo and rename the new one
            sfNominal.Delete()
            sfNominal = sfNominalNew
            
        sfNominalStat = sfNominal.Clone("sfNominalStat")
        sfNominalSys = sfNominal.Clone("sfNominalSys")
        # Assign the bin errors to the bin content
        for i in range(1, sfNominalStat.GetNbinsX()+1):
            for j in range(1, sfNominalStat.GetNbinsY()+1):
                sfNominalStat.SetBinContent(i, j, sfNominal.GetBinError(i, j))
                sfNominalSys.SetBinError(i, j, 0)
        for unc in uncList:
            # We take one uncertainty, do the difference with the nominal and then we add it in quadrature to the binError of the nominal
            sfUnc = outroot.Get(savingNames[channel][idx] + unc)
            # Compute the difference
            sfDiff = sfUnc.Clone("sfDiff")
            sfDiff.Add(sfNominal, -1.)
            # Add the difference in quadrature to the binError of the nominal
            for i in range(1, sfNominal.GetNbinsX()+1):
                for j in range(1, sfNominal.GetNbinsY()+1):
                    sfNominal.SetBinError(i, j, math.sqrt(sfNominal.GetBinError(i, j)**2 + sfDiff.GetBinContent(i, j)**2))
                    sfNominalSys.SetBinError(i, j, math.sqrt(sfNominalSys.GetBinError(i, j)**2 + sfDiff.GetBinContent(i, j)**2))
                    

        #### Now iterate over the systematic uncertainties and statistical and reassign the sfNominal errors to be stat only if stat>sys and sys+stat in quadrature if sys>stat
        ###for i in range(1, sfNominal.GetNbinsX()+1):
        ###    for j in range(1, sfNominal.GetNbinsY()+1):
        ###        if sfNominalStat.GetBinContent(i, j) > sfNominalSys.GetBinError(i, j):
        ###            sfNominal.SetBinError(i, j, sfNominalStat.GetBinContent(i, j))
        ###        else:
        ###            continue
        
        # Save the nominal with the systematics
        sfNominal.SetName(savingNames[channel][idx])
        sfNominal.Write()

        # Plot the final SF with its error (2 histograms, one for the SF and one for the error), in a separate folder
        c = r.TCanvas("c", "c", 1400, 900)
        # Set margins
        c.SetBottomMargin(0.15)
        c.SetLeftMargin(0.2)
        c.SetRightMargin(0.15)
        c.cd()

        # Scale factors
        sfNominal.SetTitle("Scale Factors")
        sfNominal.GetXaxis().SetTitle(xAxisTitles[channel][idx])
        sfNominal.GetYaxis().SetTitle(yAxisTitles[channel][idx])
        sfNominal.Draw("colz TEXT45")
        doSpam(str(lumidict[year]) + " fb^{-1} (13.6 TeV)",0.65, .913, .975, .935,textSize = spamsSize*0.98)
        # Create a folder for the final SFs
        if not os.path.exists(outpath + "/finalSFs"):
            os.makedirs(outpath + "/finalSFs")
        c.SaveAs(outpath + "/finalSFs/sf_" + channel + str(idx) + year + ".pdf")
        c.SaveAs(outpath + "/finalSFs/sf_" + channel + str(idx) + year + ".png")

        # And the uncertainty
        sfUnc = sfNominal.Clone("sfUnc")
        for i in range(1, sfUnc.GetNbinsX()+1):
            for j in range(1, sfUnc.GetNbinsY()+1):
                sfUnc.SetBinContent(i, j, sfNominal.GetBinError(i, j))
        
        sfUnc.SetTitle("Scale Factor Uncertainties")
        sfUnc.GetXaxis().SetTitle(xAxisTitles[channel][idx])
        sfUnc.GetYaxis().SetTitle(yAxisTitles[channel][idx])
        sfUnc.Draw("colz TEXT45")
        doSpam(str(lumidict[year]) + " fb^{-1} (13.6 TeV)",0.65, .913, .975, .935,textSize = spamsSize*0.98)
        c.SaveAs(outpath + "/finalSFs/sfUnc_" + channel + str(idx) + year + ".pdf")
        c.SaveAs(outpath + "/finalSFs/sfUnc_" + channel + str(idx) + year + ".png")

        # And the systematic uncertainty
        sfSysUnc = sfNominalSys.Clone("sfSysUnc")
        for i in range(1, sfSysUnc.GetNbinsX()+1):
            for j in range(1, sfSysUnc.GetNbinsY()+1):
                sfSysUnc.SetBinContent(i, j, sfNominalSys.GetBinError(i, j))

        sfSysUnc.SetTitle("Scale Factor Systematic Uncertainties")
        sfSysUnc.GetXaxis().SetTitle(xAxisTitles[channel][idx])
        sfSysUnc.GetYaxis().SetTitle(yAxisTitles[channel][idx])
        sfSysUnc.Draw("colz TEXT45")
        doSpam(str(lumidict[year]) + " fb^{-1} (13.6 TeV)",0.65, .913, .975, .935,textSize = spamsSize*0.98)
        c.SaveAs(outpath + "/finalSFs/sfSysUnc_" + channel + str(idx) + year + ".pdf")
        c.SaveAs(outpath + "/finalSFs/sfSysUnc_" + channel + str(idx) + year + ".png")

        # And the statistical uncertainty
        sfStatUnc = sfNominalStat.Clone("sfStatUnc")        
        sfStatUnc.SetTitle("Scale Factor Statistical Uncertainties")
        sfStatUnc.GetXaxis().SetTitle(xAxisTitles[channel][idx])
        sfStatUnc.GetYaxis().SetTitle(yAxisTitles[channel][idx])
        sfStatUnc.Draw("colz TEXT45")
        doSpam(str(lumidict[year]) + " fb^{-1} (13.6 TeV)",0.65, .913, .975, .935,textSize = spamsSize*0.98)
        c.SaveAs(outpath + "/finalSFs/sfStatUnc_" + channel + str(idx) + year + ".pdf")
        c.SaveAs(outpath + "/finalSFs/sfStatUnc_" + channel + str(idx) + year + ".png")
    outroot.Close()

    

    # Now plot the 1D histograms
    for idx in range(len(xAxisTitles1D[channel])):
        # Get the 1D histograms data and MC
        histoName = histoNames1D[channel][idx]
        metOnlyData   = metOnly.Get("{name}_data".format(name = histoName))
        metOnlyMC     = metOnly.Get("{name}_ttbar".format(name = histoName))

        metLeptonData = metLepton.Get("{name}_data".format(name = histoName))
        metLeptonMC   = metLepton.Get("{name}_ttbar".format(name = histoName))

        avoidUnconsistentEff1D(metLeptonMC, metOnlyMC)
        # Efficiency histograms (ratio of metLeptonData to metOnlyData)
        effData = metLeptonData.Clone("effData")
        effData.Divide(metOnlyData)

        effMC = metLeptonMC.Clone("effMC")
        effMC.Divide(metOnlyMC)

        # Compute the uncertainties of the efficiencies by using the Clopper-Pearson method
        efficiencyDataCP = r.TEfficiency(metLeptonData, metOnlyData)
        efficiencyMCCP = r.TEfficiency(metLeptonMC, metOnlyMC)

        # Create a clone of the eff histograms to store the uncertainties from TEfficiency
        effDataUncCP = effData.Clone("effDataUncCP")
        effMCUncCP = effMC.Clone("effMCUncCP")

        for i in range(1, efficiencyDataCP.GetPassedHistogram().GetNbinsX()+1):
            dataUnc = (efficiencyDataCP.GetEfficiencyErrorLow(efficiencyDataCP.GetGlobalBin(i)) + efficiencyDataCP.GetEfficiencyErrorUp(efficiencyDataCP.GetGlobalBin(i))) / 2.
            mcUnc = (efficiencyMCCP.GetEfficiencyErrorLow(efficiencyMCCP.GetGlobalBin(i)) + efficiencyMCCP.GetEfficiencyErrorUp(efficiencyMCCP.GetGlobalBin(i))) / 2.
            effDataUncCP.SetBinContent(i, dataUnc)
            effMCUncCP.SetBinContent(i, mcUnc)
            # Set the error in the original histograms
            effData.SetBinError(i, dataUnc)
            effMC.SetBinError(i, mcUnc)
        


        # Compute the scale factors (ratio of effData to effMC)
        sf = effData.Clone("sf")
        sf.Divide(effMC)

        # Compute the uncertainties of the scale factors by propagating the uncertainties of the efficiencies
        sfUncCP = sf.Clone("sfUncCP")
        for i in range(1, sfUncCP.GetNbinsX()+1):
            # Check if the denominator is zero
            if effMC.GetBinContent(i) == 0. or effData.GetBinContent(i) == 0.:
                sfUncCP.SetBinContent(i, 0.)
            else:
                sfUncCP.SetBinContent(i, sf.GetBinContent(i) * math.sqrt((effDataUncCP.GetBinContent(i) / effData.GetBinContent(i))**2 + (effMCUncCP.GetBinContent(i) / effMC.GetBinContent(i))**2))
                # Set the error in the original histograms
                sf.SetBinError(i, sfUncCP.GetBinContent(i))



        # Plot the scale factors and efficiencies
        c = r.TCanvas("c", "c", 1400, 900)
        # Set margins
        c.SetBottomMargin(0.15)
        c.SetLeftMargin(0.2)
        c.SetRightMargin(0.15)
        c.cd()



        # Divide the canvas in two pads (one for the scale factors and one for the efficiencies)
        # First pad for efficiencies with 80% of the canvas
        pad1 = r.TPad("pad1", "pad1", 0., 0.3, 1., 1.)
        pad1.SetBottomMargin(0.03)
        pad1.SetLeftMargin(0.2)
        pad1.SetRightMargin(0.15)
        # Set grid
        pad1.SetGridx()
        pad1.SetGridy()
        pad1.Draw()
        pad1.cd()


        # Efficiencies
        effData.SetTitle("")
        effData.GetXaxis().SetTitle(xAxisTitles1D[channel][idx])
        effData.GetYaxis().SetTitle(yAxisTitles1D[channel][idx])
        effData.GetYaxis().SetTitleOffset(1.2)
        effData.GetYaxis().SetTitleSize(0.05)
        effData.GetYaxis().SetLabelSize(0.05)
        effData.GetXaxis().SetTitleSize(0.05)
        effData.GetXaxis().SetLabelSize(0.05)
        effData.GetXaxis().SetTitleOffset(1.2)
        effData.GetYaxis().SetRangeUser(0.8, 1.1)
        # Remove tick labels from the x axis
        effData.GetXaxis().SetLabelOffset(999)
        # Change the color to blue
        effData.SetLineColor(r.kBlue)
        effData.SetMarkerColor(r.kBlue)
        effData.Draw("same")
        # Change the color to red
        effMC.SetLineColor(r.kRed)
        effMC.SetMarkerColor(r.kRed)
        effMC.Draw("same")
        effData.Draw("same axis")
        doSpam('#splitline{#scale[1.1]{#bf{CMS}}}{#scale[0.9]{#it{Preliminary}}}',.23, .805, .35, .845,textSize = spamsSize)
        doSpam(str(lumidict[year]) + " fb^{-1} (13.6 TeV)",0.65, .913, .975, .935,textSize = spamsSize*0.98)
        # Add the legend
        leg = r.TLegend(0.7, 0.7, 0.9, 0.9)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.05)
        leg.AddEntry(effData, "Data", "lp")
        leg.AddEntry(effMC, "MC", "lp")
        leg.Draw()

        # Second pad for scale factors with 20% of the canvas
        c.cd()
        pad2 = r.TPad("pad2", "pad2", 0., 0., 1., 0.3)
        pad2.SetTopMargin(0.05)
        pad2.SetBottomMargin(0.3)
        pad2.SetLeftMargin(0.2)
        pad2.SetRightMargin(0.15)
        # Set grid
        pad2.SetGridx()
        pad2.SetGridy()
        pad2.Draw()
        pad2.cd()
        
        # Scale factors
        sf.SetTitle("")
        sf.GetXaxis().SetTitle(xAxisTitles1D[channel][idx])
        sf.GetYaxis().SetTitle("Scale Factor")
        sf.GetYaxis().SetTitleOffset(0.5)
        sf.GetYaxis().SetTitleSize(0.05)
        sf.GetYaxis().SetLabelSize(0.05)
        sf.GetXaxis().SetTitleSize(0.05)
        sf.GetXaxis().SetLabelSize(0.05)
        sf.GetXaxis().SetTitleOffset(1.2)
        sf.GetYaxis().SetRangeUser(0.95, 1.05)
        # Center the title
        sf.GetYaxis().CenterTitle()
        # Make the axis labels and tick labels the same size as in the first pad
        sf.GetXaxis().SetLabelSize(pad1.GetHNDC() / pad2.GetHNDC() * sf.GetXaxis().GetLabelSize())
        sf.GetYaxis().SetLabelSize(pad1.GetHNDC() / pad2.GetHNDC() * sf.GetYaxis().GetLabelSize())
        sf.GetXaxis().SetTitleSize(pad1.GetHNDC() / pad2.GetHNDC() * sf.GetXaxis().GetTitleSize())
        sf.GetYaxis().SetTitleSize(pad1.GetHNDC() / pad2.GetHNDC() * sf.GetYaxis().GetTitleSize())
        # Change the ticks to match the ones in the first pad
        sf.GetXaxis().SetTickLength(0.1)
        sf.GetXaxis().SetNdivisions(510)
        sf.GetYaxis().SetTickLength(0.03)
        sf.GetYaxis().SetNdivisions(505)
        # Change the color to blue
        sf.SetLineColor(r.kBlue)
        sf.SetMarkerColor(r.kBlue)
        sf.Draw("same")
        sf.Draw("same axis")

        # Draw a line at 1
        line = r.TLine(sf.GetXaxis().GetXmin(), 1., sf.GetXaxis().GetXmax(), 1.)
        line.SetLineColor(r.kRed)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()

        # Save the canvas
        # Create a folder for the 1D histograms
        if not os.path.exists(outpath + "/1D"):
            os.makedirs(outpath + "/1D")
        c.SaveAs(outpath + "/1D/sf_" + channel + str(idx) + year + ".pdf")
        c.SaveAs(outpath + "/1D/sf_" + channel + str(idx) + year + ".png")



    # Now plot the correlation plots only for MC
    # correlation = (effMC_Dilepton * effMC_MET) / effMC_DileptonMET
    # Four folders to compute everything: LeptonTrig, METTrig, METandLeptonTrig, None

    for idx in range(len(xAxisTitlesCorr[channel])):
        # Get the 1D histograms data and MC
        histoName = histoNamesCorr[channel][idx]
        metOnlyMC     = metOnly.Get("{name}_ttbar".format(name = histoName))
        metLeptonMC   = metLepton.Get("{name}_ttbar".format(name = histoName))
        noneOnlyMC    = noneOnly.Get("{name}_ttbar".format(name = histoName))
        leptonOnlyMC  = leptonOnly.Get("{name}_ttbar".format(name = histoName))

        avoidUnconsistentEff1D(metOnlyMC, noneOnlyMC)
        avoidUnconsistentEff1D(metLeptonMC, noneOnlyMC)
        avoidUnconsistentEff1D(leptonOnlyMC, noneOnlyMC)

        # Efficiency histograms (ratio of x to noneOnlyMC)
        effMET = metOnlyMC.Clone("effMET")
        effMET.Divide(noneOnlyMC)

        effDilepton = leptonOnlyMC.Clone("effDilepton")
        effDilepton.Divide(noneOnlyMC)

        effDileptonMET = metLeptonMC.Clone("effDileptonMET")
        effDileptonMET.Divide(noneOnlyMC)

        # Compute the correlation
        correlation = effDilepton.Clone("correlation")
        correlation.Multiply(effMET)
        correlation.Divide(effDileptonMET)

        # Plot the efficiency and correlation in the same canvas (as it was done with the scale factors)
        c = r.TCanvas("c", "c", 1400, 900)
        # Set margins
        c.SetBottomMargin(0.15)
        c.SetLeftMargin(0.2)
        c.SetRightMargin(0.15)
        c.cd()

        # Divide the canvas in two pads (one for the scale factors and one for the efficiencies)
        # First pad for efficiencies with 80% of the canvas
        pad1 = r.TPad("pad1", "pad1", 0., 0.3, 1., 1.)
        pad1.SetBottomMargin(0.03)
        pad1.SetLeftMargin(0.2)
        pad1.SetRightMargin(0.15)
        # Set grid
        pad1.SetGridx()
        pad1.SetGridy()
        pad1.Draw()
        pad1.cd()

        # Efficiencies
        effMET.SetTitle("")
        effMET.GetXaxis().SetTitle(xAxisTitlesCorr[channel][idx])
        effMET.GetYaxis().SetTitle(yAxisTitlesCorr[channel][idx])
        effMET.GetYaxis().SetTitleOffset(1.2)
        effMET.GetYaxis().SetTitleSize(0.05)
        effMET.GetYaxis().SetLabelSize(0.05)
        effMET.GetXaxis().SetTitleSize(0.05)
        effMET.GetXaxis().SetLabelSize(0.05)
        effMET.GetXaxis().SetTitleOffset(1.2)
        effMET.GetYaxis().SetRangeUser(0, 1.4)
        # Remove tick labels from the x axis
        effMET.GetXaxis().SetLabelOffset(999)
        # Change the color to blue
        effMET.SetLineColor(r.kBlue)
        effMET.SetMarkerColor(r.kBlue)
        effMET.Draw("same")
        # Change the color to red
        effDilepton.SetLineColor(r.kRed)
        effDilepton.SetMarkerColor(r.kRed)
        effDilepton.Draw("same")
        # Change the color to green
        effDileptonMET.SetLineColor(r.kGreen)
        effDileptonMET.SetMarkerColor(r.kGreen)
        effDileptonMET.Draw("same")
        effMET.Draw("same axis")
        doSpam('#splitline{#scale[1.1]{#bf{CMS}}}{#scale[0.9]{#it{Preliminary}}}',.23, .805, .35, .845,textSize = spamsSize)
        doSpam(str(lumidict[year]) + " fb^{-1} (13.6 TeV)",0.65, .913, .975, .935,textSize = spamsSize*0.98)
        # Add the legend
        leg = r.TLegend(0.5, 0.7, 0.7, 0.9)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.05)
        leg.AddEntry(effMET, "MET", "lp")
        leg.AddEntry(effDilepton, "Single + double lepton", "lp")
        leg.AddEntry(effDileptonMET, "(Single + double lepton) & MET", "lp")
        leg.Draw()

        # Second pad for scale factors with 20% of the canvas
        c.cd()
        pad2 = r.TPad("pad2", "pad2", 0., 0., 1., 0.3)
        pad2.SetTopMargin(0.05)
        pad2.SetBottomMargin(0.3)
        pad2.SetLeftMargin(0.2)
        pad2.SetRightMargin(0.15)
        # Set grid
        pad2.SetGridx()
        pad2.SetGridy()
        pad2.Draw()
        pad2.cd()

        # Correlation
        correlation.SetTitle("")
        correlation.GetXaxis().SetTitle(xAxisTitlesCorr[channel][idx])
        correlation.GetYaxis().SetTitle("Correlation")
        correlation.GetYaxis().SetTitleOffset(0.5)
        correlation.GetYaxis().SetTitleSize(0.05)
        correlation.GetYaxis().SetLabelSize(0.05)
        correlation.GetXaxis().SetTitleSize(0.05)
        correlation.GetXaxis().SetLabelSize(0.05)
        correlation.GetXaxis().SetTitleOffset(1.2)
        correlation.GetYaxis().SetRangeUser(0.98, 1.02)
        # Center the title
        correlation.GetYaxis().CenterTitle()
        # Make the axis labels and tick labels the same size as in the first pad
        correlation.GetXaxis().SetLabelSize(pad1.GetHNDC() / pad2.GetHNDC() * correlation.GetXaxis().GetLabelSize())
        correlation.GetYaxis().SetLabelSize(pad1.GetHNDC() / pad2.GetHNDC() * correlation.GetYaxis().GetLabelSize())
        correlation.GetXaxis().SetTitleSize(pad1.GetHNDC() / pad2.GetHNDC() * correlation.GetXaxis().GetTitleSize())
        correlation.GetYaxis().SetTitleSize(pad1.GetHNDC() / pad2.GetHNDC() * correlation.GetYaxis().GetTitleSize())
        # Change the ticks to match the ones in the first pad
        correlation.GetXaxis().SetTickLength(0.1)
        correlation.GetXaxis().SetNdivisions(510)
        correlation.GetYaxis().SetTickLength(0.03)
        correlation.GetYaxis().SetNdivisions(505)
        # Change the color to blue
        correlation.SetLineColor(r.kBlue)
        correlation.SetMarkerColor(r.kBlue)
        correlation.Draw("same")
        correlation.Draw("same axis")

        # Draw a line at 1
        line = r.TLine(correlation.GetXaxis().GetXmin(), 1., correlation.GetXaxis().GetXmax(), 1.)
        line.SetLineColor(r.kRed)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()

        # Save the canvas
        # Create a folder for the correlation plots
        if not os.path.exists(outpath + "/corr"):
            os.makedirs(outpath + "/corr")
        c.SaveAs(outpath + "/corr/corr_" + channel + str(idx) + year + ".pdf")
        c.SaveAs(outpath + "/corr/corr_" + channel + str(idx) + year + ".png")


    # Now plot the efficiency of the single, double lepton triggers and the or of them
    for idx in range(len(xAxisTitles1D[channel])):
        # Get the 1D histograms data and MC
        histoName = histoNames1D[channel][idx]
        try:
            noneNoMETOnlyMC    = noneNoMETOnly.Get("{name}_ttbar".format(name = histoName))
            singleLeptonNoMETOnlyMC = singleLeptonNoMETOnly.Get("{name}_ttbar".format(name = histoName))
            doubleLeptonNoMETOnlyMC = doubleLeptonNoMETOnly.Get("{name}_ttbar".format(name = histoName))
            leptonNoMETOnlyMC = leptonNoMETOnly.Get("{name}_ttbar".format(name = histoName))
        except:
            print("WARNING: Not plotting single and double lepton triggers for channel " + channel + " because they don't exist")
            break
        avoidUnconsistentEff1D(singleLeptonNoMETOnlyMC, noneNoMETOnlyMC)
        avoidUnconsistentEff1D(doubleLeptonNoMETOnlyMC, noneNoMETOnlyMC)
        avoidUnconsistentEff1D(leptonNoMETOnlyMC, noneNoMETOnlyMC)

        # Efficiency histograms (ratio of x to noneNoMETOnlyMC)
        effSingleLeptonNoMET = singleLeptonNoMETOnlyMC.Clone("effSingleLeptonNoMET")
        effSingleLeptonNoMET.Divide(noneNoMETOnlyMC)

        effDoubleLeptonNoMET = doubleLeptonNoMETOnlyMC.Clone("effDoubleLeptonNoMET")
        effDoubleLeptonNoMET.Divide(noneNoMETOnlyMC)

        effLeptonNoMET = leptonNoMETOnlyMC.Clone("effLeptonNoMET")
        effLeptonNoMET.Divide(noneNoMETOnlyMC)

        # Plot the efficiency and the ratio of x to effLeptonNoMET (which is the or) in the same canvas (as it was done with the scale factors)
        c = r.TCanvas("c", "c", 1400, 900)
        # Set margins
        c.SetBottomMargin(0.15)
        c.SetLeftMargin(0.2)
        c.SetRightMargin(0.15)
        c.cd()

        # Divide the canvas in two pads (one for the scale factors and one for the efficiencies)
        # First pad for efficiencies with 80% of the canvas
        pad1 = r.TPad("pad1", "pad1", 0., 0.3, 1., 1.)
        pad1.SetBottomMargin(0.03)
        pad1.SetLeftMargin(0.2)
        pad1.SetRightMargin(0.15)
        # Set grid
        pad1.SetGridx()
        pad1.SetGridy()
        pad1.Draw()
        pad1.cd()

        # Efficiencies
        effSingleLeptonNoMET.SetTitle("")
        effSingleLeptonNoMET.GetXaxis().SetTitle(xAxisTitles1D[channel][idx])
        effSingleLeptonNoMET.GetYaxis().SetTitle(yAxisTitles1D[channel][idx])
        effSingleLeptonNoMET.GetYaxis().SetTitleOffset(1.2)
        effSingleLeptonNoMET.GetYaxis().SetTitleSize(0.05)
        effSingleLeptonNoMET.GetYaxis().SetLabelSize(0.05)
        effSingleLeptonNoMET.GetXaxis().SetTitleSize(0.05)
        effSingleLeptonNoMET.GetXaxis().SetLabelSize(0.05)
        effSingleLeptonNoMET.GetXaxis().SetTitleOffset(1.2)
        effSingleLeptonNoMET.GetYaxis().SetRangeUser(0.6, 1.2)
        # Remove tick labels from the x axis
        effSingleLeptonNoMET.GetXaxis().SetLabelOffset(999)
        # Change the color to blue
        effSingleLeptonNoMET.SetLineColor(r.kBlue)
        effSingleLeptonNoMET.SetMarkerColor(r.kBlue)
        effSingleLeptonNoMET.Draw("same")
        # Change the color to red
        effDoubleLeptonNoMET.SetLineColor(r.kRed)
        effDoubleLeptonNoMET.SetMarkerColor(r.kRed)
        effDoubleLeptonNoMET.Draw("same")
        # Change the color to green
        effLeptonNoMET.SetLineColor(r.kGreen)
        effLeptonNoMET.SetMarkerColor(r.kGreen)
        effLeptonNoMET.Draw("same")
        effSingleLeptonNoMET.Draw("same axis")
        doSpam('#splitline{#scale[1.1]{#bf{CMS}}}{#scale[0.9]{#it{Preliminary}}}',.23, .805, .35, .845,textSize = spamsSize)
        doSpam(str(lumidict[year]) + " fb^{-1} (13.6 TeV)",0.65, .913, .975, .935,textSize = spamsSize*0.98)
        # Add the legend
        leg = r.TLegend(0.5, 0.7, 0.8, 0.9)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.05)
        leg.AddEntry(effSingleLeptonNoMET, "Single lepton", "lp")
        leg.AddEntry(effDoubleLeptonNoMET, "Double lepton", "lp")
        leg.AddEntry(effLeptonNoMET, "Single + double lepton", "lp")
        leg.Draw()
        
        # Second pad for scale factors with 20% of the canvas
        c.cd()
        pad2 = r.TPad("pad2", "pad2", 0., 0., 1., 0.3)
        pad2.SetTopMargin(0.05)
        pad2.SetBottomMargin(0.3)
        pad2.SetLeftMargin(0.2)
        pad2.SetRightMargin(0.15)
        # Set grid
        pad2.SetGridx()
        pad2.SetGridy()
        pad2.Draw()
        pad2.cd()

        # Ratio of x to effLeptonNoMET
        ratioSingleLeptonNoMET = effSingleLeptonNoMET.Clone("ratioSingleLeptonNoMET")
        ratioSingleLeptonNoMET.Divide(effLeptonNoMET)

        ratioDoubleLeptonNoMET = effDoubleLeptonNoMET.Clone("ratioDoubleLeptonNoMET")
        ratioDoubleLeptonNoMET.Divide(effLeptonNoMET)

        # Ratio of x to effLeptonNoMET
        ratioDoubleLeptonNoMET.SetTitle("")
        ratioDoubleLeptonNoMET.GetXaxis().SetTitle(xAxisTitles1D[channel][idx])
        ratioDoubleLeptonNoMET.GetYaxis().SetTitle("Ratio")
        ratioDoubleLeptonNoMET.GetYaxis().SetTitleOffset(0.5)
        ratioDoubleLeptonNoMET.GetYaxis().SetTitleSize(0.05)
        ratioDoubleLeptonNoMET.GetYaxis().SetLabelSize(0.05)
        ratioDoubleLeptonNoMET.GetXaxis().SetTitleSize(0.05)
        ratioDoubleLeptonNoMET.GetXaxis().SetLabelSize(0.05)
        ratioDoubleLeptonNoMET.GetXaxis().SetTitleOffset(1.2)
        ratioDoubleLeptonNoMET.GetYaxis().SetRangeUser(0.9, 1.05)
        # Center the title
        ratioDoubleLeptonNoMET.GetYaxis().CenterTitle()
        # Make the axis labels and tick labels the same size as in the first pad
        ratioDoubleLeptonNoMET.GetXaxis().SetLabelSize(pad1.GetHNDC() / pad2.GetHNDC() * ratioDoubleLeptonNoMET.GetXaxis().GetLabelSize())
        ratioDoubleLeptonNoMET.GetYaxis().SetLabelSize(pad1.GetHNDC() / pad2.GetHNDC() * ratioDoubleLeptonNoMET.GetYaxis().GetLabelSize())
        ratioDoubleLeptonNoMET.GetXaxis().SetTitleSize(pad1.GetHNDC() / pad2.GetHNDC() * ratioDoubleLeptonNoMET.GetXaxis().GetTitleSize())
        ratioDoubleLeptonNoMET.GetYaxis().SetTitleSize(pad1.GetHNDC() / pad2.GetHNDC() * ratioDoubleLeptonNoMET.GetYaxis().GetTitleSize())
        # Change the ticks to match the ones in the first pad
        ratioDoubleLeptonNoMET.GetXaxis().SetTickLength(0.1)
        ratioDoubleLeptonNoMET.GetXaxis().SetNdivisions(510)
        ratioDoubleLeptonNoMET.GetYaxis().SetTickLength(0.03)
        ratioDoubleLeptonNoMET.GetYaxis().SetNdivisions(505)
        # Change the color to red
        ratioDoubleLeptonNoMET.SetLineColor(r.kRed)
        ratioDoubleLeptonNoMET.SetMarkerColor(r.kRed)
        ratioDoubleLeptonNoMET.Draw("same")
        ratioDoubleLeptonNoMET.Draw("same axis")

        # Change the color to blue
        ratioSingleLeptonNoMET.SetLineColor(r.kBlue)
        ratioSingleLeptonNoMET.SetMarkerColor(r.kBlue)
        ratioSingleLeptonNoMET.Draw("same")

        # Save the canvas
        # Create a folder for the correlation plots
        if not os.path.exists(outpath + "/singleVsDouble"):
            os.makedirs(outpath + "/singleVsDouble")
        c.SaveAs(outpath + "/singleVsDouble/efficiency_" + channel + str(idx) + year + ".pdf")
        c.SaveAs(outpath + "/singleVsDouble/efficiency_" + channel + str(idx) + year + ".png")









