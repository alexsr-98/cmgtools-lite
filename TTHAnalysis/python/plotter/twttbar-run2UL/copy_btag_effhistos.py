import ROOT as r
from copy import deepcopy
import os

tmpfil   = "./temp_Run2_plots/2024_09_08_btagEff_ttbar/{year}/eff/"
savepath = "{cmsswpath}/src/CMGTools/TTHAnalysis/data/TopRun2UL/btagging/".format(cmsswpath = os.environ['CMSSW_BASE'])

algodict = {"deepjet" : "DeepFlav",
            }    
wp = "M"
writethings = []
for y in ["2016apv", "2016", "2017", "2018"]:
    tmpf = r.TFile(tmpfil.format(year=y) + "/output.root", "READ")
    for algo in ["deepjet"]:
        for el in ["B", "C", "L"]:
            print("btageff_{a}_{t}_btag_{t}_ttbar".format(t = el, a = algo))
            writethings.append(deepcopy(tmpf.Get("btageff_{a}_{t}_btag_{t}_ttbar".format(t = el, a = algo)).Clone("BtagSF{t}_{a}{w}_{y}".format(t = el, a = algodict[algo], w = wp, y = y))))
    tmpf.Close()


outF = r.TFile(savepath + "/btagEffs_2024_09_08_btagEff.root", "RECREATE")
for el in writethings: el.Write()
outF.Close()
