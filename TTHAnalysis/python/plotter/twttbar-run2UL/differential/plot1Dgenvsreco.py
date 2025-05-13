import ROOT as r
import os, sys
from copy import deepcopy


textSize = 0.05616*0.4



processDictionary_Lep1_Pt = {"x_dy"            : [r.TColor.GetColor(87,144,252), "DY"],
                     "x_vvttv"         : [r.TColor.GetColor(150,74,139), "vvttv"],
                     "x_nonworz"       : [ r.TColor.GetColor(156,156,161), "non-W/Z"],
                     "x_bb4l_bkg"      : [2, "bb4l non-fiducial"],
                     "x_bb4l_partbin0" : [40, "bb4l [25, 35]"],
                     "x_bb4l_partbin1" : [41, "bb4l [35, 45]"],
                     "x_bb4l_partbin2" : [42, "bb4l [45, 55]"],
                     "x_bb4l_partbin3" : [43, "bb4l [55, 65]"],
                     "x_bb4l_partbin4" : [44, "bb4l [65, 75]"],
                     "x_bb4l_partbin5" : [45, "bb4l [75, 85]"],
                     "x_bb4l_partbin6" : [46, "bb4l [85, 95]"],
                     "x_bb4l_partbin7" : [47, "bb4l [95, 150]"],
}   

processDictionary_Jet1_Pt = {"x_dy"            : [r.TColor.GetColor(87,144,252), "DY"],
                     "x_vvttv"         : [r.TColor.GetColor(150,74,139), "vvttv"],
                     "x_nonworz"       : [ r.TColor.GetColor(156,156,161), "non-W/Z"],
                     "x_bb4l_bkg"      : [2, "bb4l non-fiducial"],
                     "x_bb4l_partbin0" : [40, "bb4l [30, 50]"],
                     "x_bb4l_partbin1" : [41, "bb4l [50, 70]"],
                     "x_bb4l_partbin2" : [42, "bb4l [70, 85]"],
                     "x_bb4l_partbin3" : [43, "bb4l [85, 110]"],
                     "x_bb4l_partbin4" : [44, "bb4l [110, 150]"],
}   
[0., 40., 60., 80., 100., 120., 140., 160., 180., 200., 220., 240., 
                       270., 310., 380., 420.],
processDictionary_minimax = {"x_dy"            : [r.TColor.GetColor(87,144,252), "DY"],
                     "x_vvttv"         : [r.TColor.GetColor(150,74,139), "vvttv"],
                     "x_nonworz"       : [ r.TColor.GetColor(156,156,161), "non-W/Z"],
                     "x_bb4l_bkg"      : [2, "bb4l non-fiducial"],
                     "x_bb4l_partbin0" : [40, "bb4l [0, 40]"],
                     "x_bb4l_partbin1" : [41, "bb4l [40, 60]"],
                     "x_bb4l_partbin2" : [42, "bb4l [60, 80]"],
                     "x_bb4l_partbin3" : [43, "bb4l [80, 100]"],
                     "x_bb4l_partbin4" : [44, "bb4l [100, 120]"],
                     "x_bb4l_partbin5" : [45, "bb4l [120, 140]"],
                     "x_bb4l_partbin6" : [46, "bb4l [140, 160]"],
                     "x_bb4l_partbin7" : [47, "bb4l [160, 180]"],
                     "x_bb4l_partbin8" : [48, "bb4l [180, 200]"],
                     "x_bb4l_partbin9" : [49, "bb4l [200, 220]"],
                     "x_bb4l_partbin10" : [30, "bb4l [220, 240]"],
                     "x_bb4l_partbin11" : [31, "bb4l [240, 270]"],
                     "x_bb4l_partbin12" : [32, "bb4l [270, 310]"],
                     "x_bb4l_partbin13" : [33, "bb4l [310, 380]"],
                     "x_bb4l_partbin14" : [38, "bb4l [380, 420]"],
}   

variablesDictionary = {"Lep1_Pt" : processDictionary_Lep1_Pt,
                       "Jet1_Pt" : processDictionary_Jet1_Pt,
                       "minimax_ATLAS" : processDictionary_minimax,
}

variablesNameDictionary = {"Lep1_Pt" : "Leading lepton #it{p}_{T} (GeV)",
                           "Jet1_Pt" : "Leading jet #it{p}_{T} (GeV)",
                           "minimax_ATLAS" : "#it{m}^{minimax} (GeV)",
}

def createFolder(folder):
    if not os.path.exists(folder): os.system("mkdir -p %s"%folder)

def createCanvas():
    c = r.TCanvas("canvas", "",  600, 600)
    topSpamSize     = 1.1
    c.SetTopMargin(c.GetTopMargin() * topSpamSize)
    c.Divide(1,2)    		        
    
    #p1 = c.GetPad(1)    		        
    #p1.SetPad(0, 0.25, 1, 1)
    #p1.SetTopMargin(0.055)
    #p1.SetBottomMargin(0.025)
    #p1.SetLeftMargin(0.16)
    #p1.SetRightMargin(0.03)
	#      		        
    #p2 = c.GetPad(2)
    #p2.SetPad(0, 0, 1, 0.25)
    #p2.SetTopMargin(0.06)
    #p2.SetBottomMargin(0.42)
    #p2.SetLeftMargin(0.16)
    #p2.SetRightMargin(0.03)

    return c

def createLegend():
    legend = r.TLegend(.3, .9, .93, .6) #0.85 for the first number in the CMGTools plotter
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetShadowColor(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(textSize)
    legend.SetNColumns(3)
    return legend

def addProcess(processDict, process, file, stackHisto, legend):
    histo = file.Get(process)
    histo.SetLineColor(1)
    histo.SetFillColor(processDict[process][0])
    stackHisto.Add(histo)
    legend.AddEntry(histo,processDict[process][1],"f")
    return histo


if __name__ == "__main__":
    variable = ["Lep1_Pt", "Jet1_Pt", "minimax_ATLAS"][0]
    rootDatacard_path = "/nfs/fanae/user/asoto/Proyectos/bb4lRun2/CMSSW_13_3_3/src/CMGTools/TTHAnalysis/python/plotter/temp_Run2_cards_diff/2025-03-31_combine_geq2jetsbjets/run2/{var}/sigextr_fit_combine/bincards/card_{var}.root".format(var = variable)
    outputPath = "./temp_tests/2025-03-31_combine_geq2jetsbjets/"


    processDictionary = variablesDictionary[variable]

    rootDatacard_file = r.TFile.Open(rootDatacard_path)


    canvas = createCanvas()
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.03)
    if variable == "minimax_ATLAS":
        canvas.SetLogy()
    legend = createLegend()

    hstack = r.THStack()


    for key in processDictionary:
        histo = addProcess(processDictionary, key, rootDatacard_file, hstack, legend)
    
    
    hstack.Draw("hist")
    hstack.GetYaxis().SetTitle("Events")
    hstack.GetXaxis().SetTitle(variablesNameDictionary[variable])
    legend.Draw("same")   

    createFolder(outputPath)
    canvas.SaveAs("%s/%s.png"%(outputPath, variable))
    canvas.SaveAs("%s/%s.pdf"%(outputPath, variable))


    rootDatacard_file.Close()
