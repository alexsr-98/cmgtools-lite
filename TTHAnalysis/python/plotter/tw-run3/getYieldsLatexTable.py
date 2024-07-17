#import ROOT as r
import warnings as wr
import sys, os, copy, math
from collections import OrderedDict
from math import log10, floor

outputFolder = "/mnt_pool/c3_users/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/temp_Run3_plots/tables/"

translateDict = OrderedDict(
    [("1j1b"               , "/mnt_pool/c3_users/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/temp_Run3_plots/2024_04_26_forAN_TopPtRew13p6TeV/run3/1j1t/tot_weight.txt"),
    ("1j1b_Differential"   , "/mnt_pool/c3_users/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/temp_Run3_plots/2024_04_26_forAN_TopPtRew13p6TeV/run3/1j1t/differential/tot_weight.txt"),
    ("2j1b"                , "/mnt_pool/c3_users/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/temp_Run3_plots/2024_04_26_forAN_TopPtRew13p6TeV/run3/2j1t/tot_weight.txt"),
    ("2j2b"                , "/mnt_pool/c3_users/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/temp_Run3_plots/2024_04_26_forAN_TopPtRew13p6TeV/run3/2j2t/tot_weight.txt")]
)

translateDict_processes = OrderedDict(
    [("tW"         , "tW"),
    ("t\\bar{t}"       , "t$\\bar{\\mathrm{t}}$"),
    ("DY"          , "Drell-Yan"),
    ("VV+t\\bar{t}V"       , "VV$+\\mathrm{t}\\bar{\\mathrm{t}}$V"),
    ("Non-W/Z"     , "Non-W/Z"),
    ("DATA"        , "Data")]
)

def round_to_reference(x, y=2):
    
    if str(y)[0] == "1" or str(y)[0] == "0" and str(y)[2] == "1":
        if y>=2:
            return int(round(x, 1-int(floor(log10(y)))))
        else:
            return round(x, 1-int(floor(log10(y))))
    else:
        if y>=2:
            return int(round(x, 0-int(floor(log10(y)))))
        else:
            return round(x, 0-int(floor(log10(y))))


def getYieldLaTeXtables(RegionsandPath, ProcessesDict, OutputFolder):
    for region in RegionsandPath:
        table = "\\begin{tabular}{ccccc} \n \hline \n%25s & %10s & %10s & %10s & %10s \\\ \hline \n" %("Process", "Yields", "Stat. Unc", "Sys. Unc", "Total Unc")
        YieldsFile = open(RegionsandPath[region],"r")
        totalYields = 0
        StatUncTotalYields = 0
        SysUncTotalYields = 0
        TotUncTotalYields = 0
        for line in YieldsFile:
            for process in ProcessesDict:
                if line.split()[0] == process and process != "DATA":
                    #Asignaciones de yields e incertidumbres
                    yields = float(line.split()[1])
                    statUnc = float(line.split()[3])
                    sysUnc = float(line.split()[6])
                    totalUnc = float(line.split()[10])
                    totalYields += yields
                    StatUncTotalYields += statUnc**2
                    SysUncTotalYields += sysUnc**2
                    TotUncTotalYields += totalUnc**2
                    #Redondeo
                    yields = str(round_to_reference(yields))
                    statUnc = str(round_to_reference(statUnc))
                    sysUnc = str(round_to_reference(sysUnc))
                    totalUnc = str(round_to_reference(totalUnc))
                    table = table + "%25s & %10s & %10s & %10s & %10s \\\ \n" %(ProcessesDict[process], yields, statUnc, sysUnc, totalUnc)
                if line.split()[0] == "DATA":
                    #Asignaciones de yields e incertidumbres
                    yieldsData = float(line.split()[1])
                    statUncData = yieldsData**0.5
                else: continue
        
        #Ratio data-MC
        ratio = yieldsData/totalYields
        StatUncRatio = 1/totalYields * (statUncData**2 + ratio**2 * StatUncTotalYields)**0.5
        SysUncRatio = 1/totalYields * ratio * SysUncTotalYields**0.5
        TotalUncRatio = (StatUncRatio**2 + SysUncRatio**2)**0.5
        ratio = str(round_to_reference(ratio, TotalUncRatio))
        StatUncRatio = str(round_to_reference(StatUncRatio, StatUncRatio))
        SysUncRatio = str(round_to_reference(SysUncRatio, SysUncRatio))
        TotalUncRatio = str(round_to_reference(TotalUncRatio, TotalUncRatio))
        
        StatUncTotalYields = StatUncTotalYields**0.5
        SysUncTotalYields = SysUncTotalYields**0.5
        TotUncTotalYields = TotUncTotalYields**0.5
        #Redondeo
        totalYields = str(round_to_reference(totalYields))
        StatUncTotalYields = str(round_to_reference(StatUncTotalYields))
        SysUncTotalYields = str(round_to_reference(SysUncTotalYields))
        TotUncTotalYields = str(round_to_reference(TotUncTotalYields))
        table = table + "%25s & %10s & %10s & %10s & %10s \\\ \n" %("\\textbf{Total}", totalYields, StatUncTotalYields, SysUncTotalYields, TotUncTotalYields)
        #Data
        yieldsData = str(round_to_reference(yieldsData))
        statUncData = str(round_to_reference(statUncData))
        table = table + "%25s & %10s & %10s & %10s & %10s \\\ \hline \n" %("\\textbf{Data}", yieldsData, statUncData, "0", statUncData)
        table = table + "%25s & %10s & %10s & %10s & %10s \\\ \hline \n" %("\\textbf{Data/Pred.}", ratio, StatUncRatio, SysUncRatio, TotalUncRatio)
        table = table + "\end{tabular} \n"
        
        LatexTable = open(OutputFolder + region + "_Table.tex","w")
        LatexTable.write(table)
        LatexTable.close()
        
        
    return


if __name__ == "__main__":
    getYieldLaTeXtables(translateDict,translateDict_processes, outputFolder)
