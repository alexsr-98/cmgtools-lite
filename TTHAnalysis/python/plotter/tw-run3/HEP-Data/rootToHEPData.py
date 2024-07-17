from hepdata_lib import RootFileReader, Submission, Variable, Uncertainty, Table
import numpy as np
import os, sys

sys.path.append('./tw-run3')
from plot_postfit import dictBinsCenterRegions, dictRegionsYaxisLabels, dictRegionsXaxisLabels

### To do list
#   - Differential plots


#Create the submission obj
submission = Submission()


#### Configurations  Inclusive-------
pathOfCards = "/mnt_pool/c3_users/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/temp_Run3_cards/2024-04-26_topptrew13p6TeV_v4/run3/"
fitDiagnosticsName = "fitDiagnosticsrun3_1j1t2j1t2j2t.root"

processes = {"tw": "tW (signal)", "ttbar": "ttbar", "vvttv": "VV+ttbarV", "dy": "DY", "nonworz": "Non-W/Z", "total_background": "Total background", "data": "Data"}
channels = {"ch1": "a", "ch2": "b", "ch3": "c"}
channelsJets = {"ch1": "1j1b", "ch2": "2j1b", "ch3": "2j2b"}
unitsDict = {"ch1": "", "ch2": "", "ch3": "GeV"}
preOrPost = {"shapes_prefit": "prefit", "shapes_fit_s": "postfit"}
figure = {"ch1": "Figure_006-a.pdf", "ch2": "Figure_006-b.pdf", "ch3": "Figure_006-c.pdf"}
positionCh = {"ch1": "top left", "ch2": "top right", "ch3": "bottom"}
numberBins = {"ch1": 20, "ch2": 12, "ch3": 16}
dictBinsExtremesRegions = {
    "ch1"      : [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 12), (12, 13), (13, 14), (14, 15), (15, 16), (16, 17), (17, 18), (18, 19), (19, 20)],
    "ch2"      : [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 12)],
    #"ch3"      : [(30, 35), (40, 45), (50, 55), (60, 65), (70, 75), (80, 85), (90, 95), (100, 105), (110, 115), (120, 125), (130, 135), (140, 145), (150, 155), (160, 165), (170, 175), (180, 185)],
    "ch3"      : [(30, 40), (40, 50), (50, 60), (60, 70), (70, 80), (80, 90), (90, 100), (100, 110), (110, 120), (120, 130), (130, 140), (140, 150), (150, 160), (160, 170), (170, 180), (180, 190)],
    #"ch1"      : [1,2],
    #"ch2"      : [1,2],
    #"ch3"      : [70,150],
}
dictRegionsXaxisLabels = {
    "ch1"      : "RF discriminant",
    "ch2"      : "RF discriminant",
    "ch3"      : "Subleading jet $p_{T}$",
}

#### -------

#### Configurations  Differential-------
pathOfDiffCards = "/mnt_pool/c3_users/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/temp_cards_diff/2024-04-26_TopPtRew13p6TeV/run3/"
histosFile = "/particlefidbinOutput.root"
variablesOrdered = ["Lep1_Pt", "Lep1Lep2Jet1_Pz","Jet1_Pt", "Lep1Lep2Jet1_M", "Lep1Lep2_DPhi", "Lep1Lep2Jet1MET_Mt"]
variablesOrdered2 = ["Lep1_Pt", "Jet1_Pt", "Lep1Lep2_DPhi", "Lep1Lep2Jet1_Pz", "Lep1Lep2Jet1MET_Mt", "Lep1Lep2Jet1_M"]
variables = {"Lep1_Pt": "$Leading$ $lepton$ $p_{T}$", "Lep1Lep2Jet1_Pz": "$p_{Z}(e^{\pm}, \mu^{\pm}, j)$","Jet1_Pt": "$Jet$ $p_{T}$", "Lep1Lep2Jet1_M": "$m(e^{\pm}, \mu^{\pm}, j)$", "Lep1Lep2_DPhi": "$\Delta\phi(e^{\pm}, \mu^{\pm})/\pi$", "Lep1Lep2Jet1MET_Mt": "$m_{T}(e^{\pm}, \mu^{\pm}, j, p_{T}^{miss})$"}
unitsDiffCov = {"Lep1_Pt": "$1/GeV^{2}$", "Lep1Lep2Jet1_Pz": "$1/GeV^{2}$", "Jet1_Pt": "$1/GeV^{2}$", "Lep1Lep2Jet1_M": "$1/GeV^{2}$", "Lep1Lep2_DPhi": "1", "Lep1Lep2Jet1MET_Mt": "$1/GeV^{2}$"}
unitsDiffCov2 = {"Lep1_Pt": "$1/GeV^{2}$", "Lep1Lep2Jet1_Pz": "$1/GeV^{2}$", "Jet1_Pt": "$1/GeV^{2}$", "Lep1Lep2Jet1_M": "$1/GeV^{2}$", "Lep1Lep2_DPhi": "", "Lep1Lep2Jet1MET_Mt": "$1/GeV^{2}$"}
processesDiff = {"tw": "tW PH DR + P8", "twds": "tW PH DS + P8", "twherwig": "tW PH DR + H7", "twamcatnlo_dr": "tW aMC DR + P8", "twamcatnlo_dr2": "tW aMC DR2 + P8", "twamcatnlo_ds": "tW aMC DS + P8", "twamcatnlo_ds_runningBW": "tW aMC DS dyn. + P8"} #Here the data has the same name as the variable
figureDiff = {"Jet1_Pt": "Figure_009-c.pdf", "Lep1Lep2Jet1MET_Mt": "Figure_009-f.pdf", "Lep1Lep2Jet1_Pz": "Figure_009-b.pdf", "Lep1Lep2Jet1_M": "Figure_009-d.pdf", "Lep1Lep2_DPhi": "Figure_009-e.pdf", "Lep1_Pt": "Figure_009-a.pdf"}
unitsDiff = {"Jet1_Pt": "GeV", "Lep1Lep2Jet1MET_Mt": "GeV", "Lep1Lep2Jet1_Pz": "GeV", "Lep1Lep2Jet1_M": "GeV", "Lep1Lep2_DPhi": "", "Lep1_Pt": "GeV"}
unitsDiffY = {"Jet1_Pt": "1/GeV", "Lep1Lep2Jet1MET_Mt": "1/GeV", "Lep1Lep2Jet1_Pz": "1/GeV", "Lep1Lep2Jet1_M": "1/GeV", "Lep1Lep2_DPhi": "", "Lep1_Pt": "1/GeV"}
names = {"Jet1_Pt": "c", "Lep1Lep2Jet1MET_Mt": "f", "Lep1Lep2Jet1_Pz": "b", "Lep1Lep2Jet1_M": "d", "Lep1Lep2_DPhi": "e", "Lep1_Pt": "a"}
positionDiff = {"Lep1_Pt": "top left", "Lep1Lep2Jet1_Pz": "top right","Jet1_Pt": "middle left", "Lep1Lep2Jet1_M": "middle right", "Lep1Lep2_DPhi": "bottom left", "Lep1Lep2Jet1MET_Mt": "bottom right"}

#### -------

pathOfPaperFig = "/nfs/fanae/user/asoto/Proyectos/tW-Run3/TOP-23-008/"

def createIndependentVar(xBins, name, is_binned, units):
    indep = Variable(name, is_independent=True, is_binned=is_binned, units=units)
    
    if type(xBins) == type([1,2]):
        indep.values = xBins
    
    else:
        if is_binned:
            edges = "_edges"
        else:
            edges = ""
        indep.values = xBins["x"+edges]
    
    return indep
    
def createDependentVar(process, name, is_binned, units, nBins,uncType, symm, uncHisto=False, alsoStat = None, depLabel = "y"):
    dep = Variable(name, is_independent=False, is_binned=is_binned, units=units)
    if nBins == -1:
        nBins = len(process[depLabel])
    
    dep.values = process[depLabel][0:nBins]
    
    # Add uncertainty
    if uncType == None:
        return dep
    
    if alsoStat != None:
        uncStat = Uncertainty("Statistical uncertainty", is_symmetric = symm)
        uncStat.values = alsoStat
        dep.add_uncertainty(uncStat)    

    if uncHisto:
        unc = Uncertainty(uncType, is_symmetric = symm)
        unc.values = uncHisto
        dep.add_uncertainty(unc)
    else:
        unc = Uncertainty(uncType, is_symmetric = symm)
        unc.values = process["d"+depLabel][0:nBins]
        dep.add_uncertainty(unc)
    
    return dep


if __name__=="__main__":
    
    # Inclusive cross section results ----
    reader = RootFileReader(pathOfCards + fitDiagnosticsName)
    for folder in ["shapes_fit_s"]: #["shapes_prefit", "shapes_fit_s"]:
        for ch in channels:
            ## Table
            table = Table("Figure 6" + channels[ch])
            #table.description = preOrPost[folder].capitalize() + " distribution of the " + dictRegionsXaxisLabels[ch] + " in the " + channelsJets[ch] + "  region. The uncertainty band includes the postfit uncertainties. The bottom of each panel shows the ratios of data to the predictions after the fit (points)."
            table.description = "The distribution of the " + dictRegionsXaxisLabels[ch] + " for events in the " + channelsJets[ch] + "  region. The data (points) and the MC predictions (filled histograms) after the maximum likelihood fit are shown. The vertical bars on the points represent the statistical uncertainty in the data, and the hatched band the total uncertainty in the MC prediction. The lower panels display the ratio of the data to the sum of the MC (points) predictions after the fit, with the bands giving the corresponding uncertainties."
            table.location = "Data from Figure 6 (" + positionCh[ch] + "), located on page 17."
            #table.keywords["observables"] = ["N"]
            table.add_image(pathOfPaperFig + figure[ch])
            
            nameX, nameY, units = dictRegionsXaxisLabels[ch], dictRegionsYaxisLabels[ch] + " ", unitsDict[ch]
            table.add_variable(createIndependentVar(dictBinsExtremesRegions[ch], nameX, True, units))
    
            for process in processes:
                if process == "data":
                    dataGraph = reader.read_graph(folder + "/" + ch + "/" + process)
                    table.add_variable(createDependentVar(dataGraph, nameY + processes[process], False, "", numberBins[ch],"Statistical uncertainty", False))
                else:
                    histogram = reader.read_hist_1d(folder + "/" + ch + "/" + process)
                    table.add_variable(createDependentVar(histogram, nameY + processes[process], False, "", numberBins[ch], "Total uncertainty", True))
            ## Add table and create the yaml files
            submission.add_table(table)
    
    # ----
    table = Table("Cross section")
    table.location = "Data from page 16."
    table.description = "Observed and theoretical cross sections. In the observed, the first uncertainty is statistical, the second is the systematic, and the third the luminosity. In the theoretical, the first uncertainty is due to scale variations, the second due to the choice of PDF. The theoretical cross section has been computed at approximate third-order in QCD, with the addition of third-order corrections of soft gluon emission terms, assuming a top quark mass of 172.5 GeV and using the PDF4LHC21 PDF set."
    sigmaTheo = Variable("$\sigma_{tW}^{SM}$ ($\sqrt{s}=13.6$ TeV)", is_independent=False, is_binned=False, units="pb")
    sigmaTheoUncScale = Uncertainty("scale", is_symmetric = False)
    sigmaTheoUncPDF = Uncertainty("PDF+$\\alpha_{S}$", is_symmetric = True)
    sigmaTheo.values = [87.9]
    sigmaTheoUncScale.values = [(-1.9,2.0)]
    sigmaTheoUncPDF.values = [2.4]
    sigmaTheo.add_uncertainty(sigmaTheoUncScale)
    sigmaTheo.add_uncertainty(sigmaTheoUncPDF)
    
    sigmaExp = Variable("$\sigma_{tW}^{Obs.}$ ($\sqrt{s}=13.6$ TeV)", is_independent=False, is_binned=False, units="pb")
    sigmaExpUncStat = Uncertainty("stat", is_symmetric = True)
    sigmaExpUncSys = Uncertainty("syst", is_symmetric = False)
    sigmaExpUncLum = Uncertainty("lumi", is_symmetric = True)
    sigmaExp.values = [82.3]
    sigmaExpUncStat.values = [2.1]
    sigmaExpUncSys.values = [(-9.8,10.2)]
    sigmaExpUncLum.values = [3.3]
    sigmaExp.add_uncertainty(sigmaExpUncStat)
    sigmaExp.add_uncertainty(sigmaExpUncSys)
    sigmaExp.add_uncertainty(sigmaExpUncLum)
    
    table.add_variable(sigmaExp)
    table.add_variable(sigmaTheo)
    submission.add_table(table)
    
    # Yields
    table = Table("Table 3")
    table.description = "The number of observed and MC predicted events after the fit in the 1j1b, 2j1b, and 2j2b regions. The total uncertainties in the predictions are given."
    table.location = "Data from Table 3, located on page 16."
    col1 = Variable("Process", is_independent=True, is_binned=False, units="")
    processesList = ["tW", "ttbar", "Drell-Yan", "VV+ttbarV", "Non-W/Z", "Total", "Data"]
    col1.values = processesList
    table.add_variable(col1)
    
    processesList_Yields = {"1j1b": [8000,49200,670,460,340,58700,58635], 
                            "2j1b": [3670,42520,330,450,810,47780,47810],
                            "2j2b": [1140,33380,42,190,64,24810,34818],
                            #"Postfit/prefit" : [1.11,0.91,1.09,1.18,1.68,0.94,""],
                            }
    
    processesList_unc = {"1j1b": [300,300,60,50,70,150,0], 
                            "2j1b": [160,190,40,70,50,130,0],
                            "2j2b": [70,160,7,30,12,130,0],
                           # "Postfit/prefit" : ["","","","","","",""],
                            }
  
    channelslist = ["1j1b", "2j1b", "2j2b"]#, "Postfit/prefit"]
    for ch in channelslist:
        ## Table
        nameX, nameY, units = "", "Events", ""
        cols = Variable(ch, is_independent=False, is_binned=False, units="")
        cols.values = processesList_Yields[ch]
        uncTot = Uncertainty("Total uncertainty", is_symmetric = True)
        uncTot.values = processesList_unc[ch]
        cols.add_uncertainty(uncTot)   
        table.add_variable(cols)  
        
        
    ## Add table and create the yaml files
    submission.add_table(table)
    
    
    # Differential cross section ----
    for folder in variablesOrdered:
        processesDiff[folder] = "Data"
        reader = RootFileReader(pathOfDiffCards + folder + histosFile)
        ## Table
        table = Table("Figure 9" + names[folder])
        table.description = "Normalised fiducial differential tW production cross section as a function of the " + variables[folder] + ". The horizontal bars on the points show the bin width. Predictions from POWHEG (PH) + PYTHIA 8 (P8) DR and DS, POWHEG + HERWIG 7 (H7) DR, MADGRAPH5_aMC@NLO (aMC) + PYTHIA 8 DR, DR2, DS and DS with a dynamic factor are also shown. The grey band represents the statistical uncertainty and the orange band the total uncertainty. In the lower panels, the ratio of the predictions to the data is shown."
        table.location = "Data from Figure 9 (" + positionDiff[folder] + "), located on page 20."
        #table.keywords["observables"] = ["N"]
        table.add_image(pathOfPaperFig + figureDiff[folder])
        
        nameX, nameY, unitsX, unitsY = variables[folder], "$(1/\sigma_{fid.})d\sigma/d$" + variables[folder], unitsDiff[folder], unitsDiffY[folder]
        
    
        for process in processesDiff:
            if processesDiff[process] == "Data":
                dataGraph = reader.read_hist_1d(folder)
                table.add_variable(createIndependentVar(dataGraph, nameX, True, unitsX))
                histSysUp = reader.read_hist_1d(folder + "_totalUp")
                histSysDown = reader.read_hist_1d(folder + "_totalDown")
                histSysUpVal = histSysUp["dy"]
                histSysDownVal = histSysDown["dy"]
                histStatUp = reader.read_hist_1d(folder + "_statUp")
                histStatDown = reader.read_hist_1d(folder + "_statDown")
                histStatUpVal = histStatUp["dy"]
                histStatDownVal = histStatDown["dy"]
                for nbin in range(0,len(histSysUpVal)):
                    histSysUpVal[nbin] = abs(histSysUpVal[nbin])
                    histSysDownVal[nbin] = -abs(histSysDownVal[nbin])
                    histStatUpVal[nbin] = abs(histStatUpVal[nbin])
                    histStatDownVal[nbin] = -abs(histStatDownVal[nbin])

                histoSys = list(zip(histSysUpVal,histSysDownVal))
                histoStat = list(zip(histStatUpVal,histStatDownVal))
                table.add_variable(createDependentVar(dataGraph, nameY + " " + processesDiff[process], False, unitsY, -1,"Total uncertainty", False, histoSys, histoStat))
                
            else:
                histogram = reader.read_hist_1d(process)
                table.add_variable(createDependentVar(histogram, nameY + " " + processesDiff[process], False, unitsY, -1, None, True))
        
        del processesDiff[folder]
        ## Add table and create the yaml files
        submission.add_table(table)  
    

    
    # Create the chi2 tables for the differential measurements
    processesDiff = {"tru": "PH DR + P8", "tru_DS": "PH DS + P8", "tru_herwig": "PH DR + H7", "tru_aMC_dr": "aMC DR + P8", "tru_aMC_dr2": "aMC DR2 + P8", "tru_aMC_ds": "aMC DS + P8", "tru_aMC_ds_runn": "aMC DS dyn. + P8"}
    
    ## Table
    table = Table("Tables 4 and 5")
    table.description = "The p-values from the goodness-of-fit tests comparing the six differential cross section measurements with the predictions from POWHEG (PH) + PYTHIA 8 (P8) DR and DS,  POWHEG + HERWIG 7 (H7) DR, MADGRAPH5 aMC@NLO (aMC) + PYTHIA 8 DR, DR2, DS, and DS with a dynamic factor. The complete covariance matrix from the results and the statistical uncertainties in the predictions are taken into account."
    table.location = "Data from Tables 4 and 5 located on page 18 and 21."

    col1 = Variable("Variable", is_independent=True, is_binned=False, units="")
    variablesTable = []
    for folder in variablesOrdered2:

        #table.keywords["observables"] = ["N"]
        
        nameX, nameY, unitsX, unitsY = "", "$\chi^2$", "", ""
        variablesTable.append(variables[folder])   
    
    col1.values = variablesTable
    table.add_variable(col1)       
    
    cols_vals = {"tru":             [0.92,0.93,0.79,0.86,0.93,0.88], 
                 "tru_DS":          [0.94,0.98,0.83,0.88,0.96,0.74], 
                 "tru_herwig":      [0.92,0.94,0.77,0.85,0.93,0.92], 
                 "tru_aMC_dr":      [0.91,0.91,0.75,0.89,0.92,0.98], 
                 "tru_aMC_dr2":     [0.95,0.93,0.79,0.86,0.96,0.88], 
                 "tru_aMC_ds":      [0.92,0.94,0.79,0.86,0.94,0.94], 
                 "tru_aMC_ds_runn": [0.93,0.98,0.76,0.85,0.95,0.93]}
    for process in processesDiff:
        cols = Variable(processesDiff[process], is_independent=False, is_binned=False, units="")
        cols.values = cols_vals[process]
        table.add_variable(cols)                    
                
        
        #del processesDiff[folder]
        ## Add table and create the yaml files
      
    submission.add_table(table) 
    # ----

    ## Suplementary material
    figureprefit = {"ch1": "prefit_ch1.pdf", "ch2": "prefit_ch2.pdf", "ch3": "prefit_ch3.pdf"}
    pathOfSuplementaryFigures = "/mnt_pool/c3_users/user/asoto/Proyectos/tW-Run3/TOP-23-008/supplemental_figures/"
    reader = RootFileReader(pathOfCards + fitDiagnosticsName)
    for folder in ["shapes_prefit"]: #["shapes_prefit", "shapes_fit_s"]:
        for ch in channels:
            ## Table
            table = Table("Supplemental material: prefit " + channelsJets[ch] + " " + dictRegionsXaxisLabels[ch])
            #table.description = preOrPost[folder].capitalize() + " distribution of the " + dictRegionsXaxisLabels[ch] + " in the " + channelsJets[ch] + "  region. The uncertainty band includes the postfit uncertainties. The bottom of each panel shows the ratios of data to the predictions after the fit (points)."
            table.description = "Distribution of the " + dictRegionsXaxisLabels[ch] + " for events in the " + channelsJets[ch] + "  region. The data (points) and the MC predictions (filled histograms) before the maximum likelihood fit are shown. The vertical bars on the points represent the statistical uncertainty in the data, and the hatched band the total uncertainty in the MC prediction. The lower panels display the ratio of the data to the sum of the MC (points) predictions before the fit, with the bands giving the corresponding uncertainties."
            table.location = "Supplemental material."
            #table.keywords["observables"] = ["N"]
            table.add_image(pathOfSuplementaryFigures + figureprefit[ch])
            
            nameX, nameY, units = dictRegionsXaxisLabels[ch], dictRegionsYaxisLabels[ch] + " ", unitsDict[ch]
            table.add_variable(createIndependentVar(dictBinsExtremesRegions[ch], nameX, True, units))
    
            for process in processes:
                if process == "data":
                    dataGraph = reader.read_graph(folder + "/" + ch + "/" + process)
                    table.add_variable(createDependentVar(dataGraph, nameY + processes[process], False, "", numberBins[ch],"Statistical uncertainty", False))
                else:
                    histogram = reader.read_hist_1d(folder + "/" + ch + "/" + process)
                    table.add_variable(createDependentVar(histogram, nameY + processes[process], False, "", numberBins[ch], "Total uncertainty", True))
            ## Add table and create the yaml files
            submission.add_table(table)


    pathOfSuplementary = "/mnt_pool/c3_users/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/temp_cards_diff/2024-04-26_TopPtRew13p6TeV/run3/"
    listOfObs = ["Lep1_Pt", "Lep1Lep2Jet1_Pz", "Jet1_Pt", "Lep1Lep2Jet1_M", "Lep1Lep2_DPhi", "Lep1Lep2Jet1MET_Mt"]
    templateName = "R{Var}_"
    responseMatrices = "UnfoldingInfo.root"
    covMatrices = "CovMat_particlefidbin_v2.root"
    
    for var in listOfObs:
        reader = RootFileReader(pathOfSuplementary + var + "/" + responseMatrices)
        ## Table
        table = Table("Supplemental material: " + var + " response matrix")
        table.description = "Response matrix between detector and particle level for the " +  variables[var] + "."
        table.location = "Supplemental material."
        
        table.add_image(pathOfSuplementaryFigures + "Rnonumb_{Var}_.pdf".format(Var = var))
        
        nameX, nameY, unitsX, unitsY = variables[var] + "[particle level]", variables[var] + "[detector level]", unitsDiff[var], unitsDiff[var]
        
        histogram = reader.read_hist_2d(templateName.format(Var=var))
        table.add_variable(createIndependentVar(histogram["x_edges"], nameX, True, unitsX))
        table.add_variable(createIndependentVar(histogram["y_edges"], nameY, True, unitsY))
        table.add_variable(createDependentVar(histogram, "Response", False, "", -1, "Statistical uncertainty", True,depLabel="z"))
        submission.add_table(table)


    for var in listOfObs:
        reader = RootFileReader(pathOfSuplementary + var + "/" + covMatrices)
        ## Table
        table = Table("Supplemental material: " + var + " covariance matrix")
        table.description = "Covariance matrix including all uncertainties of the normalised differential cross section in bins of " + variables[var]  + " in units of " + unitsDiffCov[var] + "."
        table.location = "Supplemental material."
        table.add_image(pathOfSuplementaryFigures + "Cov_{Var}_fidbin.pdf".format(Var = var))
        
        nameX, nameY, unitsX, unitsY = variables[var], variables[var]+ " ", unitsDiff[var], unitsDiff[var]
        
        histogram = reader.read_hist_2d("CovMat_fidbin")
        table.add_variable(createIndependentVar(histogram["x_edges"], nameX, True, unitsX))
        table.add_variable(createIndependentVar(histogram["y_edges"], nameY, True, unitsY))
        table.add_variable(createDependentVar(histogram, "Covariance", False, unitsDiffCov2[var], -1, None, True,depLabel="z"))
        submission.add_table(table)


    
    submission.create_files("./tw-run3/HEP-Data/example_output")
    os.system("mv submission.tar.gz ./tw-run3/HEP-Data/example_output")

