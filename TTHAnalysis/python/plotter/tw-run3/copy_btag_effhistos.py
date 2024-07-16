import ROOT as r
from copy import deepcopy
import os

tmpfil   = "./temp_Run3_plots/2024_02_28_btagEff_tW_comb/{year}/eff/"
savepath = "{cmsswpath}/src/CMGTools/TTHAnalysis/data/TopRun3/btagging/".format(cmsswpath = os.environ['CMSSW_BASE'])

algodict = {"deepjet" : "DeepFlav",
            "particlenet"     : "PNet",
            "particletransformer" : "RobustParTAK4"
            }    
wp = "M"
writethings = []
for y in ["2022","2022PostEE"]:
    tmpf = r.TFile(tmpfil.format(year=y) + "/output.root", "READ")
    for algo in ["deepjet", "particlenet", "particletransformer"]:
        for el in ["B", "C", "L"]:
            print("btageff_{a}_{t}_btag_{t}_tw".format(t = el, a = algo))
            writethings.append(deepcopy(tmpf.Get("btageff_{a}_{t}_btag_{t}_tw".format(t = el, a = algo)).Clone("BtagSF{t}_{a}{w}_{y}".format(t = el, a = algodict[algo], w = wp, y = y))))
    tmpf.Close()


outF = r.TFile(savepath + "/btagEffs_2024_01_16_btagEff.root", "RECREATE")
for el in writethings: el.Write()
outF.Close()
