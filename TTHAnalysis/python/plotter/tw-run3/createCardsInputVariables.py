import os, sys
import numpy as np
import ROOT as r
# -- Variables to be used in the MVA - Multiclass
inputVariables_1j1b = ["JetLoose1_Pt", "Lep1Lep2Jet1_Pt", "Mll", "Lep1Lep2_DPhi", "Lep1Jet1_Pt", "LepGood_pt_corrAll[0]", "Jet1_Pt","Lep1Lep2Jet1_M"]
inputVariables_1j1b_latexName = ["Leading loose jet $p_{T}$", "$p_T (e^{\pm}, \mu^{\mp}, j)$", "$m(e^{\pm}, \mu^{\mp})$", "$\Delta\phi(e^{\pm}, \mu^{\mp})$", "$p_T (\ell_1, j)$", "Leading lepton $p_{T}$", "Jet $p_{T}$", "$m(e^{\pm}, \mu^{\mp}, j)$"]

# -- Variables to be used in the MVA - Multiclass
inputVariables_2j1b = ["Mll", "Lep12Jet12_DR", "Lep1Lep2Jet1_Pt", "Lep1Jet1_DR", "Lep1Lep2_DR", "LepGood_pt_corrAll[1]", "Jet2_Pt","Lep1Lep2Jet1_C"]
inputVariables_2j1b_latexName = ["$m(e^{\pm}, \mu^{\mp})$", "$\Delta R(\ell_{12}, \ jet_{12})$", "$p_T (e^{\pm}, \mu^{\mp}, j)$", "$\Delta R(\ell_1, \ jet_1)$", "$\Delta R(e^{\pm}, \mu^{\mp})$", "Subleading lepton $p_{T}$", "Subleading jet $p_{T}$", "$C(e^{\pm}, \mu^{\mp}, j)$"]

### OPTIONS
outputFolder = "temp_Run3_cards/2024-03-04_GOFs"
year = "run3"
region = "1j1t"
runCards = False
runGOF = False
runPlots = True
###########

#inputVariables_1j1b = ["Lep1Jet1_Pt", "Lep1Lep2Jet1MET_M", "Lep1Lep2Jet1_Pt", "Lep1Lep2Jet1MET_PtOverHTtot", "Lep1_PtLep2_PtOverHTtot", "Lep1Jet1_DR", "Lep1Lep2Jet1MET_Mt", "Mll", "Jet1_Pt", "Lep1Lep2Jet1_C", "Lep1Lep2Jet1_Pz", "Lep1Lep2_DR", "Lep1Lep2_DPhi", "METgood_pt", "HTtot","JetLoose1_Pt"]
#inputVariables_2j1b = ["Jet2_Pt", "Lep1Jet1_DR", "Lep12Jet12_DR", "Lep1Jet1_Pt", "Lep1Lep2Jet1MET_M", "Lep1Lep2Jet1_Pt", "Lep1Lep2Jet1MET_PtOverHTtot", "Lep1_PtLep2_PtOverHTtot", "Lep1Lep2Jet1MET_Mt", "Mll", "Jet1_Pt", "Lep1Lep2Jet1_C", "Lep1Lep2Jet1_Pz", "Lep1Lep2_DR", "Lep1Lep2_DPhi", "METgood_pt", "HTtot"]

cardsCommand = "python3 tw-run3/cardsHelper.py -P 2023-12-23 -y {y} -o {o}/{v} -r {reg} -v '{var}' -b '{bins}' -q batch -j 16"
gofCommand = "python tw-run3/fitsHelper.py -y {y} -i {o}/{v} -r {reg}  --gofprefit -q batch --nToys 10000"
gofCommandPlot = "python tw-run3/fitsHelper.py -y {y} -i {o}/{v} -r {reg}  --gofplot --gofplotTitle {v}"

plots_1j1btxt = open("tw-run3/plots-tw/plots-tw-1j1t_MVAtrain.txt", "r")
plots_2j1btxt = open("tw-run3/plots-tw/plots-tw-2j1t_MVAtrain.txt", "r")

nbins = 20
formatedVar = "min(max({v}, {lowerBin}), {upperBin})"



# For the GOF table
gofs = {}

if region == "1j1t":
    print("Region 1j1t")
    inputVariables = inputVariables_1j1b
    inputVariables_latexName = inputVariables_1j1b_latexName
    plots = plots_1j1btxt
elif region == "2j1t":
    print("Region 2j1t")
    inputVariables = inputVariables_2j1b
    inputVariables_latexName = inputVariables_2j1b_latexName
    plots = plots_2j1btxt



# Read the 1j1b plot file and look for wich line contain each variable
for line in plots:
    if line.startswith("#"): continue
    if not line: continue
    if len(line) < 2: continue
    variable = line.split(";")[0].split(":")[1]
    # Remove empty spaces
    variable = variable.replace(" ", "")
    binning = line.split(";")[0].split(":")[2]
    # Remove empty spaces
    binningVals = binning.replace(" ", "").split(",")
    lowerBin = float(binningVals[1])
    upperBin = float(binningVals[2])
    listBins = str(list(np.linspace(lowerBin, upperBin, nbins+1)))
    for var in inputVariables:
        if var == variable:
            formatedVar = formatedVar.format(v=var, lowerBin=lowerBin, upperBin=upperBin)       
            print(var, binning)
            if runCards:
                print(cardsCommand.format(o=outputFolder, y=year, v=var, reg=region, var=var, bins=listBins))
                os.system(cardsCommand.format(o=outputFolder, y=year, v=var, reg=region, var=var, bins=listBins))
                print("")
            if runGOF:
                print(gofCommand.format(o=outputFolder, y=year, v=var, reg=region))
                os.system(gofCommand.format(o=outputFolder, y=year, v=var, reg=region))
            if runPlots:
                # We have to compute the nbins opening the root file
                cardsRootFile = r.TFile("{o}/{v}/{y}/{r}/cuts-tw-{r}.root".format(o=outputFolder,y=year,r=region,v=var), "READ")
                dataObs = cardsRootFile.Get("x_data_obs")
                realNbins = 0
                # Iterate over all the bins 
                for i in range(dataObs.GetNbinsX()):
                    if dataObs.GetBinContent(i+1) > 10:
                        realNbins += 1
                print(gofCommandPlot.format(o=outputFolder, y=year, v=var, reg=region, ndof=realNbins-1))
                os.system(gofCommandPlot.format(o=outputFolder, y=year, v=var, reg=region, ndof=realNbins-1))
                # Read the GOF txt file
                goftxt = open("{o}/{v}/{y}/GOF_{y}_{r}/plots/GOF_{y}_{r}.txt".format(o=outputFolder,y=year,r=region,v=var), "r")
                for line in goftxt:
                    gofs[var] = float(line)
                goftxt.close()


# Close the files
plots_1j1btxt.close()
plots_2j1btxt.close()

if runPlots:
    # Create the GOF latex table
    latex_table = "\\begin{tabular}{cc}\n\\hline\nVariable & p-value \\\\\n\\hline\n"
    for i, var in enumerate(inputVariables):
        latex_table += "{var} & {pvalue} \\\\\n".format(var=inputVariables_latexName[i], pvalue=gofs[var])
    latex_table += "\\end{tabular}\n"

    # Write the latex table
    latex_table_file = open("{o}/{y}_{r}_GOF.tex".format(o=outputFolder,y=year,r=region), "w")
    latex_table_file.write(latex_table)
    latex_table_file.close()