import os, sys
import numpy as np

# -- Variables to be used in the MVA - Multiclass
inputVariables_1j1b = ["JetLoose1_Pt", "Lep1Lep2Jet1_Pt", "Mll", "Lep1Lep2_DPhi", "Lep1Jet1_Pt", "Lep1Lep2_DR", "LepGood_pt_corrAll[0]", "Jet1_Pt"]

# -- Variables to be used in the MVA - Multiclass
inputVariables_2j1b = ["Mll", "Lep12Jet12_DR", "Lep1Lep2Jet1_Pt", "Lep1Jet1_DR", "Lep1Lep2_DR", "LepGood_pt_corrAll[1]", "Jet2_Pt"]



#inputVariables_1j1b = ["Lep1Jet1_Pt", "Lep1Lep2Jet1MET_M", "Lep1Lep2Jet1_Pt", "Lep1Lep2Jet1MET_PtOverHTtot", "Lep1_PtLep2_PtOverHTtot", "Lep1Jet1_DR", "Lep1Lep2Jet1MET_Mt", "Mll", "Jet1_Pt", "Lep1Lep2Jet1_C", "Lep1Lep2Jet1_Pz", "Lep1Lep2_DR", "Lep1Lep2_DPhi", "METgood_pt", "HTtot","JetLoose1_Pt"]
#inputVariables_2j1b = ["Jet2_Pt", "Lep1Jet1_DR", "Lep12Jet12_DR", "Lep1Jet1_Pt", "Lep1Lep2Jet1MET_M", "Lep1Lep2Jet1_Pt", "Lep1Lep2Jet1MET_PtOverHTtot", "Lep1_PtLep2_PtOverHTtot", "Lep1Lep2Jet1MET_Mt", "Mll", "Jet1_Pt", "Lep1Lep2Jet1_C", "Lep1Lep2Jet1_Pz", "Lep1Lep2_DR", "Lep1Lep2_DPhi", "METgood_pt", "HTtot"]

cardsCommand = "python3 tw-run3/cardsHelper.py -P 2023-06-02 -y 2022PostEE -o ./temp_Run3_cards/2023-10-03_GOFs/{v} -r {reg} -v '{var}' -b '{bins}' -q batch -j 16"
gofCommand = "python tw-run3/fitsHelper.py -y 2022PostEE -i temp_Run3_cards/2023-10-03_GOFs/{v} -r {reg}  --gofprefit -q batch --nToys 10000"
gofCommandPlot = "python tw-run3/fitsHelper.py -y 2022PostEE -i temp_Run3_cards/2023-10-03_GOFs/{v} -r {reg}  --gofplot --gofplotTitle {v}"

plots_1j1btxt = open("tw-run3/plots-tw/plots-tw-1j1t.txt", "r")
plots_2j1btxt = open("tw-run3/plots-tw/plots-tw-2j1t.txt", "r")

nbins = 20
formatedVar = "min(max({v}, {lowerBin}), {upperBin})"

### OPTIONS
region = "1j1t"
runCards = False
runGOF = False
runPlots = True
###########

if region == "1j1t":
    print("Region 1j1t")
    inputVariables = inputVariables_1j1b
    plots = plots_1j1btxt
elif region == "2j1t":
    print("Region 2j1t")
    inputVariables = inputVariables_2j1b
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
                print(cardsCommand.format(v=var, reg=region, var=var, bins=listBins))
                os.system(cardsCommand.format(v=var, reg=region, var=var, bins=listBins))
                print("")
            if runGOF:
                print(gofCommand.format(v=var, reg=region))
                os.system(gofCommand.format(v=var, reg=region))
            if runPlots:
                print(gofCommandPlot.format(v=var, reg=region))
                os.system(gofCommandPlot.format(v=var, reg=region))

# Read the 2j1b plot file and look for wich line contain each variable


# Close the files
plots_1j1btxt.close()
plots_2j1btxt.close()
