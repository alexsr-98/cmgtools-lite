'''
Library of all the variables so everything is centralised and there are not confusions
'''

import ROOT as r
import warnings as wr
import os, sys
from copy import deepcopy
from array import array

# === ESSENTIAL PARAMETERS OF THE ANALYSIS. CHANGING THIS APPLIES TO EVERYTHING. ===
nuncs       = 3         # Number of uncs. shown in the relative uncertainty plots
diffControlReg = "3j2t" # Control region used in the differential signal extraction step

asimov      = True     # Use of Asimov dataset or data
doxsec      = True      # Show events or diff. cross section in final results
doPre       = True      # Show or not show the "Preliminary" in the plots
#doPre       = False      # Show or not show the "Preliminary" in the plots
doSym       = True      # Symmetrise the uncertainties or not
#doSym       = False      # Symmetrise the uncertainties or not
doReg       = False     # Regularise distributions (general setting, can be overrrided for each variable separately)
doArea      = False     # Apply area constraint in unfolding (general setting, can be overrrided for each variable separately)
#onlyTotal   = True      # Only show total unc. line in the differential relative unc. plots.
onlyTotal   = False      # Only show total unc. line in the differential relative unc. plots.

vetolist = ["plots", "control", "tables", "response"]
arXivtext = "arXiv:2208.06485"

# === OTHER IMPORTANT DEFINITIONS ===
LumiDict     = {"2022"       : 7.78,
                "2022PostEE" : 20.67}
TotalLumi   = LumiDict["2022"] + LumiDict["2022PostEE"] # In femtobarns

plotlimits   = tuple([float(i) for i in "0.00, 0.25, 1.00, 1.00".split(',')]) # xlow, ylow, xup, yup
ratiolimits  = tuple([float(i) for i in "0.00, 0.00, 1.00, 0.25".split(',')]) # xlow, ylow, xup, yup
margins      = "0.06, 0.1, 0.04, 0.1" # top, bottom, right, left
marginsratio = "0.03, 0.4, 0.04, 0.1" # top, bottom, right, left
legpos       = (0.82, 0.65, 0.93, 0.93)

if asimov: labellegend = 'Asimov dataset'
else:      labellegend = 'Data'


#def GetLastFolder(stpth):
#    savefolders   = next(os.walk(stpth))[1]
#    saveyears     = list(map(int, [i[6:]  for i in savefolders]))
#    savefolders   = [i for i in savefolders if int(i[6:]) == max(saveyears)]
#    savemonths    = list(map(int, [i[3:5] for i in savefolders]))
#    savefolders   = [i for i in savefolders if int(i[3:5]) == max(savemonths)]
#    savedays      = list(map(int, [i[:2]  for i in savefolders]))
#    savefolders   = [i for i in savefolders if int(i[:2]) == max(savedays)]
#    return (stpth + savefolders[0] + "/")

def SetUpWarnings():
    wr.simplefilter("always", UserWarning)
    wr.filterwarnings(action = 'ignore', category = RuntimeWarning, message = 'TClass::Init:0: RuntimeWarning: no dictionary for class')
    return

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

#def GiveMeTheExpNamesWOJER(inl):
#    l = []
#    for el in inl:
#        if "JER" in el or "Down" in el: continue
#        l.append(el.replace("Up", ""))
#    return ",".join(l)


def giveMeOneComparison(thef, name, scalevalue, iV, part = False, normfid = False, normbin = False):
    #print name
    outH = deepcopy(thef.Get("x_" + name).Clone(name))
    outH.Scale(scalevalue)
    outH.SetMarkerSize(1.25)
    outH.SetMarkerColor(comparisonColourDict[name])
    outH.SetMarkerStyle(markersdict[name])
    outH.SetLineColor(0)
    outH.SetLineWidth(0)

    binstr = "bins_detector" if not part else "bins_particle"
    if normfid:
        fiduval = None; fiduunc = None
        htmp = deepcopy(outH.Clone(outH.GetName() + "_fid"))
        fiduval = sum([htmp.GetBinContent(i) for i in range(1, htmp.GetNbinsX() + 1)])
        fiduunc = r.TMath.Sqrt(sum([htmp.GetBinError(i)**2 for i in range(1, htmp.GetNbinsX() + 1)]))

        for bin in range(1, htmp.GetNbinsX() + 1):
            htmp.SetBinContent(bin, fiduval)
            htmp.SetBinError(  bin, fiduunc)

        thegooduncs = []
        for bin in range(1, outH.GetNbinsX() + 1):
            thegooduncs.append(r.TMath.Sqrt(outH.GetBinError(bin)**2 / fiduval**2 +
                                            outH.GetBinContent(bin)**2 * fiduunc**2 / fiduval**4 -
                                            2 * (outH.GetBinContent(bin) * outH.GetBinError(bin)**2) / fiduval**3 )
            )

        outH.Divide(htmp)

        for iB in range(1, outH.GetNbinsX() + 1):
            outH.SetBinError(iB, thegooduncs[iB - 1])
        del htmp

    if normbin:
        outH.Scale(1, "width")

    if spacingdict[name] != 0.:
        goodoutH = deepcopy(outH.Clone(name))
        smallestw = 9999
        for i in range(1, len(varList[iV][binstr])):
            thedif = varList[iV][binstr][i] - varList[iV][binstr][i - 1]
            if thedif < smallestw:
                smallestw = thedif

        newbins = array("d", [el + spacingdict[name]*smallestw/2 for el in varList[iV][binstr]])
        goodoutH.SetBins(len(varList[iV][binstr]) - 1, newbins)
        for iB in range(1, goodoutH.GetNbinsX() + 1):
            goodoutH.SetBinContent(iB, outH.GetBinContent(iB))
            goodoutH.SetBinError(iB, outH.GetBinError(iB))
        del outH
        return goodoutH
    else:
        return outH


def parseRelUncs(path, verbose = False):
    if not os.path.isfile(path):
        raise RuntimeError("FATAL: given systematic table/relative uncs. result file, {f}, does not exist.".format(f = path))
    theres = {}
    with open(path, "r") as res:
        for line in res:
            tmpline = line.strip().replace(" ", "").split(":")
            if len(tmpline) > 2:
                tmpup = float(tmpline[1]); tmpdn = float(tmpline[2]); tmpnam = tmpline[0]
                theres[tmpnam] = {}
                if tmpdn == 0. and tmpup != 0.:
                    tmpdn = tmpup
                    theres[tmpnam]["err"] = "down"
                    if verbose:
                        print("\t- WARNING: for file {f}, unc. source {u} effect estimation has failed for the \'down\' variation. This will be symmetrised.".format(f = path,
                                                                                                                                                                     u = tmpnam))
                elif tmpup == 0. and tmpdn != 0.:
                    tmpup = tmpdn
                    theres[tmpnam]["err"] = "up"
                    if verbose:
                        print("\t- WARNING: for file {f}, unc. source {u} effect estimation has failed for the \'up\' variation. This will be symmetrised.".format(f = path,
                                                                                                                                                                   u = tmpnam))
                elif tmpup == 0. and tmpdn == 0.:
                    raise RuntimeError("FATAL: for file {f}, unc. source {u} effect estimation has failed for both variations".format(f = path,
                                                                                                                                      u = tmpnam))
                
                theres[tmpnam]["up"]   = tmpup
                theres[tmpnam]["down"] = tmpdn
            else:
                theres["munom"] = float(tmpline[1])

    return theres


#def prepareCards(thef):
#    tmpdict = {}
#
#    #### Get the histos
#    tmpf = r.TFile(thef, "READ")
#    for key in tmpf.GetListOfKeys():
#        tmpdict[key.GetName()] = deepcopy(tmpf.Get(key.GetName()).Clone(key.GetName()))
#
#    tmpf.Close(); del tmpf
#
#    #### Process them
#    # tW PH+P8 DR + ttbar PH+P8
#    tmpdict["twttbardr"] = deepcopy(tmpdict["ttbar"].Clone("twttbardr"))
#    tmpdict["twttbardr"].Add(tmpdict["twdr"])
#
#    # tW PH+P8 DS + ttbar PH+P8
#    tmpdict["twttbards"] = deepcopy(tmpdict["ttbar"].Clone("twttbards"))
#    tmpdict["twttbards"].Add(tmpdict["twds"])
#
#    # tW PH+H7 DR + ttbar PH+H7
#    tmpdict["twttbarherwig"] = deepcopy(tmpdict["ttbarherwig"].Clone("twttbarherwig"))
#    tmpdict["twttbarherwig"].Add(tmpdict["twherwig"])
#
#    # tW aMC+P8 DR + ttbar aMC+P8
#    tmpdict["twttbaramc_dr"] = deepcopy(tmpdict["ttbaramc"].Clone("twttbaramc_dr"))
#    tmpdict["twttbaramc_dr"].Add(tmpdict["twamc_dr"])
#
#    # tW aMC+P8 DR2 + ttbar aMC+P8
#    tmpdict["twttbaramc_dr2"] = deepcopy(tmpdict["ttbaramc"].Clone("twttbaramc_dr2"))
#    tmpdict["twttbaramc_dr2"].Add(tmpdict["twamc_dr2"])
#
#    # tW aMC+P8 DS + ttbar aMC+P8
#    tmpdict["twttbaramc_ds"] = deepcopy(tmpdict["ttbaramc"].Clone("twttbaramc_ds"))
#    tmpdict["twttbaramc_ds"].Add(tmpdict["twamc_ds"])
#
#    # tW aMC+P8 DS IS + ttbar aMC+P8
#    tmpdict["twttbaramc_ds_is"] = deepcopy(tmpdict["ttbaramc"].Clone("twttbaramc_ds_is"))
#    tmpdict["twttbaramc_ds_is"].Add(tmpdict["twamc_ds_is"])
#
#    # tW aMC+P8 DS runningBW + ttbar aMC+P8
#    tmpdict["twttbaramc_ds_runningBW"] = deepcopy(tmpdict["ttbaramc"].Clone("twttbaramc_ds_runningBW"))
#    tmpdict["twttbaramc_ds_runningBW"].Add(tmpdict["twamc_ds_runningBW"])
#
#    # tW aMC+P8 DS IS runningBW + ttbar aMC+P8
#    tmpdict["twttbaramc_ds_is_runningBW"] = deepcopy(tmpdict["ttbaramc"].Clone("twttbaramc_ds_is_runningBW"))
#    tmpdict["twttbaramc_ds_is_runningBW"].Add(tmpdict["twamc_ds_is_runningBW"])
#
#
#    #### Delete used histograms
#    del tmpdict["ttbar"]
#    del tmpdict["ttbarherwig"]
#    del tmpdict["ttbaramc"]
#    del tmpdict["twdr"]
#    del tmpdict["twds"]
#    del tmpdict["twherwig"]
#    del tmpdict["twamc_dr"]
#    del tmpdict["twamc_dr2"]
#    del tmpdict["twamc_ds"]
#    del tmpdict["twamc_ds_runningBW"]
#    del tmpdict["twamc_ds_is"]
#    del tmpdict["twamc_ds_is_runningBW"]
#
#
#    #### Rename current file
#    os.system("mv " + thef + " " + thef.replace(".root", "_old.root"))
#
#
#    #### Save new file
#    tmpf = r.TFile(thef, "RECREATE")
#    for key in tmpdict: tmpdict[key].Write()
#    tmpf.Close(); del tmpf
#    return


def getInfoFromFD(path, fdpath, unc, iY, iV, isNom = False):
    signalname = "x_tw"
    ints = []
    card = r.TFile.Open(path + "/particle.root", "READ")
    nparticlebins = card.Get(signalname).GetNbinsX()
    scaleval = 1
    thelumi = TotalLumi if iY == "run3" else LumiDict[iY]
    if doxsec: scaleval = 1/thelumi/1000.
    for iB in range(1, nparticlebins + 1):
        ints.append(card.Get(signalname).GetBinContent(iB))
    card.Close(); del card

    tfile  = r.TFile.Open(fdpath + "/fitDiagnostics{u}.root".format(u = unc))
    tfile2 = r.TFile.Open((fdpath + "/higgsCombine{u}.FitDiagnostics.mH125.root".format(u = unc)) if os.path.isfile(fdpath + "/higgsCombine{u}.FitDiagnostics.mH125.root".format(u = unc)) else (fdpath + "/higgsCombine{u}.FitDiagnostics.mH120.root".format(u = unc)))
    poiList = r.RooArgList('poiList')
    w       = tfile2.Get('w')
    fitResult = tfile.Get('fit_s')
    for i in range(nparticlebins):
        var = w.var('r_tW_%d'%i)
        poiList.add(var)
    cov   = deepcopy(fitResult.reducedCovarianceMatrix(poiList))

    results = {}
    count = 0
    for var in fitResult.floatParsFinal():
        count += 1
        results[var.GetName()] = [ var.getVal(), abs(var.getErrorLo()), var.getErrorHi(), var.getError() ]
        # if "r_tW" in var.GetName(): print var.getVal(), abs(var.getErrorLo()), var.getErrorHi(), var.getError()
        if count == fitResult.floatParsFinal().getSize(): break

    tfile.Close();  del tfile
    tfile2.Close(); del tfile2

    outUp   = r.TH1D('{u}Up'.format(u = unc) if not isNom else "hNom",
                        '', nparticlebins, array('d', varList[iV]['bins_particle']))
    outDown = deepcopy(outUp.Clone('{u}Down'.format(u = unc)))

    hCov = r.TH2D('CovMat_{u}'.format(u = unc if not isNom else ""), '',
                  nparticlebins, -0.5, nparticlebins - 0.5,
                  nparticlebins, -0.5, nparticlebins - 0.5)
    for i in range(1, nparticlebins + 1):
        results['r_tW_%d'%(i-1)][0] *= ints[i-1] * scaleval
        results['r_tW_%d'%(i-1)][1] *= ints[i-1] * scaleval
        results['r_tW_%d'%(i-1)][2] *= ints[i-1] * scaleval
        results['r_tW_%d'%(i-1)][3] *= ints[i-1] * scaleval

        # ...and these here are asymm.:
        # if not isNom:
        #     outUp.  SetBinContent(i,   results['r_tW_%d'%(i-1)][0]) # asym
        #     outUp.  SetBinError  (i,   results['r_tW_%d'%(i-1)][2]) # asym
        #     outDown.SetBinContent(i,   results['r_tW_%d'%(i-1)][0]) # asym
        #     outDown.SetBinError  (i,   results['r_tW_%d'%(i-1)][1]) # asym
        # else:
        #     outUp.  SetBinContent(i,   results['r_tW_%d'%(i-1)][0]) # asym
        #     outUp.  SetBinError  (i,   mean([results['r_tW_%d'%(i-1)][1], results['r_tW_%d'%(i-1)][2]])) # asym

        ##### FORCED SYMM
        outUp.  SetBinContent(i,   results['r_tW_%d'%(i-1)][0]) # sym
        outUp.  SetBinError  (i,   results['r_tW_%d'%(i-1)][3]) # sym
        outDown.SetBinContent(i,   results['r_tW_%d'%(i-1)][0]) # sym
        outDown.SetBinError  (i,   results['r_tW_%d'%(i-1)][3]) # sym

        # print outUp.GetBinContent(iB), outUp.GetBinError(iB), scaleval, ints[i-1]

        for j in range(1, nparticlebins + 1):
            normx = ints[i - 1] * scaleval
            normy = ints[j - 1] * scaleval

            cov[i-1][j-1] = cov[i-1][j-1] * normx * normy
            hCov.SetBinContent( hCov.GetBin(i,j), cov[i-1][j-1] )
    
    if isNom:
        return outUp, hCov
    else:
        return outUp, outDown, hCov


def getActualUnc(D):
    #print D
    rawD = deepcopy(D)
    for i,iU in enumerate(individual_list):
        # print "\n", iU
        if (iU + "Up") not in D and (iU + "Down") not in D: raise RuntimeError("FATAL: uncertainty group from the individual list is not available in the relative unc. information needed to be conveyed to the differential plots.")

        if i == 0:
            key = ""
        else:
            key = individual_list[i - 1]

        # if key == "":
        #     print "AQUI"
        for iB in range(1, D[""].GetNbinsX() + 1):
            # print D[key + ("Up"   if key != "" else "")].GetBinError(iB), D[iU + "Up"].GetBinError(iB)
            D[iU + "Up"].  SetBinError(iB, r.TMath.Sqrt(rawD[key + ("Up"   if key != "" else "")].GetBinError(iB)**2 - D[iU + "Up"].  GetBinError(iB)**2))
            D[iU + "Down"].SetBinError(iB, r.TMath.Sqrt(rawD[key + ("Down" if key != "" else "")].GetBinError(iB)**2 - D[iU + "Down"].GetBinError(iB)**2))

            # print D[iU + "Up"].GetBinError(iB)
    return D

def getActualCovMat(D):
    rawD = deepcopy(D)
    for i,iU in enumerate(individual_list):
        #print "\n", iU
        if iU not in D: raise RuntimeError("FATAL: uncertainty group from the individual list is not available in the relative unc. information needed to be conveyed to the differential plots.")

        if i == 0:
            key = ""
        else:
            key = individual_list[i - 1]

        # if key == "":
        #     print "AQUI"

        # for iB in range(1, D[iU].GetNbinsX() + 1):
        #     print D[iU].GetBinContent(iB, iB)
        # print " "
        # for iB in range(1, D[iU].GetNbinsX() + 1):
        #     print rawD[key].GetBinContent(iB, iB)
        D[iU].Add(rawD[key], -1)
        D[iU].Scale(-1)
        # print " "
        #for iB in range(1, D[iU].GetNbinsX() + 1):
            #print D[iU].GetBinContent(iB, iB)
        # sys.exit()
    return D


def getInfoForRelUncs(path, iY, iV, verbose = False):
    #if not os.path.isfile(path):
        #raise RuntimeError("FATAL: given systematic table/relative uncs. result file, {f}, does not exist.".format(f = path))
    theres = {}; covdict = {}
    nompath = path + "/sigextr_fit_combine/fitdiagnostics/fitDiagnostics{y}_{v}.root".format(y = iY, v = iV)
    glodir  = path + "/sigextr_fit_combine/relativeuncs_fd/global/"
    inddir  = path + "/sigextr_fit_combine/relativeuncs_fd/individual/"

    # get the existent variations
    indvars = []
    for el in os.listdir(inddir):
        if "fitDiagnostics" in el:
            indvars.append(el.replace("fitDiagnostics", "").replace(".root", ""))

    # First, nominal
    theres[""], covdict[""] = deepcopy(getInfoFromFD(path,
                                       path + "/sigextr_fit_combine/fitdiagnostics",
                                       iY + "_" + iV,
                                       iY, iV,
                                       True))

    # Then, globals (stat)
    theres["statUp"], theres["statDown"], covdict["stat"] = deepcopy(getInfoFromFD(path,
                                                                    glodir,
                                                                    "systematics",
                                                                    iY, iV))

    theres["statUp"].SetName("statUp")
    theres["statDown"].SetName("statDown")
    covdict["stat"].SetName("CovMat_stat")

    # Then, syst
    hCovar_syst = deepcopy(covdict[""].Clone("CovMat_syst"))
    hCovar_syst.Add(covdict["stat"], -1)
    covdict["syst"] = hCovar_syst
    hsystup   = deepcopy(theres[""].Clone("systUp"))
    hsystdown = deepcopy(theres[""].Clone("systDown"))
    for iB in range(1, hsystup.GetNbinsX() + 1):
        hsystup.  SetBinError(iB, r.TMath.Sqrt(theres[""].GetBinError(iB)**2 - theres["statUp"].GetBinError(iB)**2))
        hsystdown.SetBinError(iB, r.TMath.Sqrt(theres[""].GetBinError(iB)**2 - theres["statDown"].GetBinError(iB)**2))
        # print covdict[""].GetBinContent(iB, iB), covdict["stat"].GetBinContent(iB, iB), covdict["syst"].GetBinContent(iB, iB)

    theres["systUp"]   = hsystup
    theres["systDown"] = hsystdown

    # Then, the individuals. First, we get the raw info.

    for iU in indvars:
        theres[iU + "Up"], theres[iU + "Down"], covdict[iU] = deepcopy(getInfoFromFD(path, inddir, iU, iY, iV))

    # And now, we process it to obtain the actual uncertainty.
    theres  = getActualUnc(theres)
    covdict = getActualCovMat(covdict)

    return theres, covdict


def getAlternateUncsHistos(nom, induncs, glouncs):
    thehistos = {}; theUncs = []
    thehistos[""] = deepcopy(nom)
    nbs = nom.GetNbinsX()

    for iU in induncs[1]:
        if iU == "" or iU == "munom": continue
        theUncs.append(iU)
    for iU in glouncs[1]:
        if iU == "" or iU == "munom": continue
        theUncs.append(iU)

    for iU in theUncs:
        tmpnam = iU.replace("inc_", "").replace("systematics", "syst")
        thehistos[tmpnam + "Down"] = deepcopy(nom.Clone(tmpnam + "Down"))
        thehistos[tmpnam + "Up"]   = deepcopy(nom.Clone(tmpnam + "Up"))
        for iB in range(1, nbs + 1):
            if "syst" == tmpnam or "stat" == tmpnam:
                thed = glouncs
            else:
                thed = induncs
            
            thehistos[tmpnam + "Down"].SetBinError(iB,  thed[iB][iU]["down"] / glouncs[iB]["munom"]  * nom.GetBinContent(iB))
            thehistos[tmpnam + "Up"].  SetBinError(iB, (thed[iB][iU]["up"]   / glouncs[iB]["munom"]) * nom.GetBinContent(iB))
            #if tmpnam == "syst":
                #print thehistos[tmpnam + "Down"].GetBinError(iB), thehistos[tmpnam + "Up"].GetBinError(iB)

    thehistos["totalUp"] = deepcopy(nom.Clone("totalUp")); thehistos["totalDown"] = deepcopy(nom.Clone("totalDown"))

    for iB in range(1, nom.GetNbinsX() + 1):
        thehistos["totalUp"].  SetBinError(iB, r.TMath.Sqrt(thehistos["statUp"].  GetBinError(iB)**2 + thehistos["systUp"].  GetBinError(iB)**2))
        thehistos["totalDown"].SetBinError(iB, r.TMath.Sqrt(thehistos["statDown"].GetBinError(iB)**2 + thehistos["systDown"].GetBinError(iB)**2))

        themax = max([thehistos["totalUp"].GetBinError(iB), thehistos["totalDown"].GetBinError(iB)])
        thehistos[""].SetBinError(iB, themax if themax < thehistos[""].GetBinContent(iB) else thehistos[""].GetBinContent(iB))

    return thehistos


def getAConservativeMaximum(h):
    theval = 0.
    for iB in range(1, h.GetNbinsX() + 1):
        newmax = (h.GetBinContent(iB) + h.GetBinError(iB)) * 1.2
        if theval < newmax: theval = newmax
    return theval


def confirm(message = "Do you wish to continue?"):
    """
    Ask user to enter y(es) or n(o) (case-insensitive).
    :return: True if the answer is Y.
    :rtype: bool
    """
    answer = ""
    while answer not in ["y", "n", "yes", "no"]:
        answer = input(message + " [Y/N]\n").lower()
    return answer[0] == "y"



# === DICTIONARIES ===
# var           := name of the variable in the tree to make cards
# var_response  := name of the variable in the response matrix without the M
varList = {}

#varList['LCurve'] = {
#    'xaxis' : '\\log_{10}\\left(\\mathscr{L}_{1}\\right)',
#    'yaxis' : '\\log_{10}\\left(\\mathscr{L}_{2}\\,/\\tau^{2}\\right)',
#}
#
#varList['LCurvature'] = {
#    'xaxis' : '\\log_{10}\\left(\\tau\\right)',
#    'yaxis' : '\\mathscr{C}\\left(\\log_{10}\\left(\\tau\\right)\\right)',
#}

#varList['LogTauCurv'] = varList['LCurvature']

varList['Names'] = {
    #'Variables'   : ["Lep1_Pt", "Jet1_Pt", "Lep1Lep2_DPhi", "Lep1Lep2Jet1_Pz", "Lep1Lep2Jet1MET_Mt", "Lep1Lep2Jet1_M"],
    'Variables'   : ["Lep1_Pt", "Jet1_Pt", "Lep1Lep2_DPhi", "Lep1Lep2Jet1_Pz", "Lep1Lep2Jet1MET_Mt", "Lep1Lep2Jet1_M", "Fiducial"],
    #'Variables'   : ["Lep1_Pt", "Jet1_Pt", "Lep1Lep2_DPhi", "Lep1Lep2Jet1_Pz", "Lep1Lep2Jet1MET_Mt", "Lep1Lep2Jet1_M"],
#    'ExpSysts'    : ["JESUp", "JESDown", "JERUp", "ElecEffUp", "ElecEffDown", "MuonEffUp", #   DO NOT MOVE THE FIRST THREE TO OTHER
#                     "MuonEffDown", "TrigUp", "TrigDown", "PUUp", "PUDown", "BtagUp",      # POSITION: it will affect the calculus
#                     #"BtagDown", "MistagUp", "MistagDown"],                                # of the response matrices.
#                     "BtagDown", "MistagUp", "MistagDown", "TopPtUp", "TopPtDown"],        # of the response matrices.  WITH TOP PT REW
#    'ttbarSysts'  : ["ttbarMEUp", "ttbarMEDown", "pdfUp", "pdfDown", "hDampUp", "hDampDown", "UEUp", "UEDown"],
#    'specialSysts': ["JERDown", "DSDown"],
#    'colorSysts'  : ["GluonMoveCRTuneerdON", "PowhegerdON", "QCDbasedCRTuneerdON", "GluonMoveCRTune"],
#    'NormSysts'   : ["ttbarUp", "ttbarDown", "Non-WorZUp", "Non-WorZDown", "DYUp", "DYDown", "VVttbarVUp", "VVttbarVDown"],
}


varList['Lep1Lep2Jet1MET_Mt'] = {
    #'xaxis'       : 'm_{T}(\\ell_{1}, \\ell_{2},\\slash{E}_{T}, j) (GeV)',
    'xaxis'       : '#it{m}_{T}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{p}_{T}^{miss}, #it{j}) (GeV)',
    'printname'   : "\\transmassvar (\GeV)",
    'printnamenodim':"\\transmassvar",
    'mathprintname': "\\transmassvar",
    'yaxis_particle'       : 'd#sigma/d(#it{m}_{T}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{p}_{T}^{miss}, #it{j})) (pb)',
    # 'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(#it{m}_{T}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{p}_{T}^{miss}, #it{j})) (adim.)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(#it{m}_{T}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{p}_{T}^{miss}, #it{j}))',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(#it{m}_{T}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{p}_{T}^{miss}, #it{j})) (1/GeV)',
    'yaxisnorm'   : 'd#sigma/d(#it{m}_{T}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{p}_{T}^{miss}, #it{j})) (pb/GeV)',
    # 'yaxis_unc'   : 'Relative uncertainty (adim.)',
    'yaxis_unc'   : 'Relative uncertainty',

    'bins_particle' : [100., 220., 280., 325., 400., 500.], # propuesta (5 bins)
    'bins_detector' : [100., 220., 250., 260., 270., 280., 290., 300., 325., 400., 500.], # propuesta (5 bins) SELECCIONADA

#    'bins_particle' : [100., 220., 280., 325., 400., 500.], # propuesta (5 bins)
#    'bins_detector' : [100., 220., 280., 325., 400., 500.],

    #'descbinning' : [0., 800.], ## TEMPORAL
    #'ndescbins'   : 20, ## TEMPORAL
    "var_detector"  : 'min(max(Lep1Lep2Jet1MET_Mt, 101.), 499.)',
#    "var_detector"  : 'min(max(Lep1Lep2Jet1MET_Mt, 101.), 599.)',
    'var_response'  : 'Lep1Lep2Jet1MET_Mt',
    'var_particle'  : 'min(max(DressLep1Lep2Jet1MET_Mt, 101.), 499.)',
    #'var_particle'  : 'min(max(DressLep1Lep2Jet1MET_Mt, 101.), 599.)',
    'legpos'      : (0.51, 0.55, 0.71, 0.93),
    #'legposdesc'  : (0.15, 0.425, 0.35, 0.81),
    "legposdesc"  : (0.65, 0.55, 0.85, 0.93),
    "maxdesc"  : 2600,
    'legpos_detectoras':"BL",
    'legpos_detector' : "BL",
    'legpos_particlefid'  : "TL",
    "legpos_particle"   : (.18, .26, .36, .03),
    #"legpos_particle"  : "BL",
    'legpos_particleas': "TL",
    'legpos_detectorunc' : "TL",
    'legpos_particleunc'  : "TC",
    'legpos_particlefidunc'  : "TL",
    #'legpos_particlefidbin':(.175, .73, .27, .49),
    #'legpos_particlefidbin': (.52, .9, .72, .65),
    # 'legpos_particlefidbin': (.66, .9, .8, .6),
    # 'legpos_particlefidbin': (.6, .9, .75, .55),
    #'legpos_particlefidbin': (.63, .88, .78, .53),
    'legpos_particlefidbin': (.51, .88, .78, .43),
    'legpos_particlefidbinunc':"TL",
    'resptxtsize' : 0.9,
    'txtsize_covdetector': 0.75,
    'txtsize_covparticle': 1.5,
    'txtangle_covparticle': 45,
    "txtsize_covparticlefidbin": 1.4,
    'txtangle_covparticlefidbin': 35,
    "yaxisuplimitunf": 0.176,
    "yaxismax_particlefid" : 1.6,
#    "yaxismax_particlefidbin": 0.007,
    "yaxismax_particlefidbin": 0.0105,
    #"yaxismax_particlefidbinunc" : 1.8,
    "yaxismax_particlefidbinunc" : 1.3,
    "yaxismax_unf" : 2,
    #"yaxismax_particlebinunc": 2.0,
    "yaxismax_particlebinunc": 2.1,
    'legpos_particlebin': "TL",
    "legpos_particlebinunc" : "TC",
    "yaxismax_particlebin": 0.0038,
    #"yaxismax_ratio_fidbin" : 4.0,
#    "yaxismax_ratio_fidbin" : 2.5,
    "yaxismax_ratio_fidbin" : 1.5,
    "yaxismax_ratio_norm" : 5.0,
    # "particle_smoothing" : {"jes_FlavorQCD" : {-1 : "fitOrder=1, symmAfterFit=1",
    #                                            4  : "NOM",},
    #                         "jes_RelativeBal" : "fitOrder=1, symmAfterFit=1",
    #                         "jes_RelativeSample" : "fitOrder=1, symmAfterFit=1",
    #                         "jes_RelativeSample_2016" : {-1 : "fitOrder=1, symmAfterFit=1",
    #                                                      0  : "NOM",},
    #                         "jes_RelativeSample_2017" : {1 : "fitOrder=1, symmAfterFit=1",
    #                                                      2 : "fitOrder=1, symmAfterFit=1",},
    #                         "jes_RelativeSample_2018" : "fitOrder=1, symmAfterFit=1",
    #                         "jes_Absolute" : "fitOrder=1, symmAfterFit=1",
    #                         "jes_Absolute_2016" : "fitOrder=1, symmAfterFit=1",
    #                         "jes_Absolute_2017" : "fitOrder=1, symmAfterFit=1",
    #                         "jes_Absolute_2018" : "fitOrder=1, symmAfterFit=1",
    #                         "unclenergy"        : {1 : "fitOrder=4, symmAfterFit=1",
    #                                                2 : "fitOrder=4, symmAfterFit=1",
    #                                                3 : "fitOrder=2, symmAfterFit=1",
    #                                                4 : "EstimateFromXbins=1",},
    # },
}


varList['Lep1Lep2Jet1_M'] = {
    #'xaxis'       : 'm(#ell_{1}, #ell_{2}, j) (GeV)',
    'xaxis'       : '#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j}) (GeV)',
    'printname'   : '\\invmassvar (\\GeV)',
    'printnamenodim':'\\invmassvar',
    'mathprintname': '\\invmassvar',
    'yaxis_particle'       : 'd#sigma/d(#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (pb)',
    # 'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (adim.)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j}))',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (1/GeV)',
    'yaxisnorm'   : 'd#sigma/d(#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (pb/GeV)',
    # 'yaxis_unc'   : 'Relative uncertainty (adim.)',
    'yaxis_unc'   : 'Relative uncertainty',
    'bins_particle' : [50., 130., 170., 200., 250., 325., 400.], # propuesta (6 bins) BUENA
    'bins_detector' : [50., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 250., 325., 400.], # propuesta (6 bins) SELECCIONADA BUENA

#    'bins_particle'  : [50., 130., 170., 200., 250., 325., 400.],
#    'bins_detector' :  [50., 130., 170., 200., 250., 325., 400.],

    'descbinning' : [0., 400.],
    'ndescbins'   : 16,
    'legpos'      : (0.75, 0.55, 0.95, 0.93),
    #'legposdesc'  : (0.11, 0.435, 0.31, 0.82),
    #'legposdesc'  : (0.65, 0.55, 0.85, 0.93),
    'legposdesc'  : (0.7, 0.55, 0.90, 0.93),
    "var_detector"         : 'min(max(Lep1Lep2Jet1_M, 51.), 399.)',
    'var_response'         : 'Lep1Lep2Jet1_M',
    'var_particle'         : 'min(max(DressLep1Lep2Jet1_M, 51.), 399.)',
    'legpos_detectorunc' : "TC",
    'legpos_particlefidunc'  : "TC",
    'legpos_particleunc'  : "TL",
    'legpos_particlefidbinunc': (.45, .615, .63, .9),
    # 'legpos_particlefidbin': (.66, .9, .8, .6),
    # 'legpos_particlefidbin': (.6, .9, .75, .55),
    #'legpos_particlefidbin': (.63, .88, .78, .53),
    'legpos_particlefidbin': (.51, .88, .78, .43),
    #'legpos_particlefidbin': (.52, .9, .72, .65),
    "legpos_particle"   : (.52, .9, .72, .65),
    'txtsize_covdetector': 0.5,
    'txtsize_covparticle': 1.5,
    'txtangle_covparticle': 35,
    'txtsize_covparticlefidbin': 1.3,
    'txtangle_covparticlefidbin': 35,
    "yaxisuplimitunf": 0.15,
    "yaxismax_particlefidunc" : 0.7,
    "yaxismax_particlefidbinunc" : 1.3,
#    "yaxismax_particlefidbin" : 0.009,
    "yaxismax_particlefidbin" : 0.01,
    "yaxismax_unf" : 2.,
    "yaxismax_particlebin": 0.004,
    "legpos_particlebinunc" : (.18, .5, .31, .785),
    "legpos_particlefidbinunc": (.4, .615, .58, .9),
#    "yaxismax_ratio_fidbin" : 2.5,
    "yaxismax_ratio_fidbin" : 1.5,
    # "particle_smoothing" : {"jes_BBEC1"       : {-1 : "fitOrder=1, symmAfterFit=1",
    #                                              0  : "NOM",
    #                                              5  : "NOM",},
    #                         "jes_FlavorQCD"   : "fitOrder=1, symmAfterFit=1",
    #                         "jes_RelativeBal" : {-1 : "fitOrder=1, symmAfterFit=1",
    #                                              5  : "NOM",},
    #                         "jes_RelativeSample_2018" : {-1 : "fitOrder=1, symmAfterFit=1",
    #                                                      0  : "fitOrder=3, symmAfterFit=1",
    #                                                      5  : "NOM",},
    #                         "jes_Absolute" : {5  : "EstimateFromXbins=1",},
    #                         "jes_Absolute_2016" : "fitOrder=1, symmAfterFit=1",
    #                         "jes_Absolute_2017" : "fitOrder=1, symmAfterFit=1",
    #                         "jes_Absolute_2018" : "fitOrder=1, symmAfterFit=1",
    #                         "jer_2016" : {-1 : "fitOrder=2, symmAfterFit=1",
    #                                       0  : "fitOrder=3, symmAfterFit=1",
    #                                       3  : "fitOrder=1, symmAfterFit=1",
    #                                       4  : "fitOrder=1, symmAfterFit=1",
    #                                       5  : "fitOrder=1, symmAfterFit=1",},
    #                         "jer_2017" : {-1 : "fitOrder=2, symmAfterFit=1",
    #                                       0  : "EstimateFromXbins=1",
    #                                       4  : "fitOrder=1, symmAfterFit=1",
    #                                       5  : "fitOrder=1, symmAfterFit=1",},
    #                         "jer_2018" : {-1 : "fitOrder=2, symmAfterFit=1",
    #                                       0  : "fitOrder=3, symmAfterFit=1",
    #                                       4  : "fitOrder=1, symmAfterFit=1",
    #                                       5  : "fitOrder=1, symmAfterFit=1",},
    #                         "ds" : {2 : "fitOrder=1, symmAfterFit=1"},
    #                        },
}


varList['Lep1Lep2Jet1_M_control'] = {
    #'xaxis'       : 'm(#ell_{1}, #ell_{2}, j) (GeV)',
    'xaxis'       : '#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j}) (GeV)',
    'printname'   : '\\invmassvar (\\GeV)',
    'printnamenodim':'\\invmassvar',
    'mathprintname': '\\invmassvar',
    'yaxis_particle'       : 'd#sigma/d(#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (pb)',
    # 'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (adim.)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j}))',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (1/GeV)',
    'yaxisnorm'   : 'd#sigma/d(#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (pb/GeV)',
    # 'yaxis_unc'   : 'Relative uncertainty (adim.)',
    'yaxis_unc'   : 'Relative uncertainty',


    'bins_particle'  : [50., 130., 170., 200., 250., 325., 400.],
    'bins_detector' : [0., 120., 180., 240., 300., 360., 420., 480., 540., 600.],

    'descbinning' : [0., 600.],
    'ndescbins'   : 16,
    'legpos'      : (0.75, 0.55, 0.95, 0.93),
    'legposdesc'  : (0.7, 0.55, 0.90, 0.93),
    "var_detector"         : 'min(max(Lep1Lep2Jet1_M, 0.), 599.)',
    'var_response'         : 'Lep1Lep2Jet1_M',
    'var_particle'         : 'min(max(DressLep1Lep2Jet1_M, 0.), 599.)',
    'legpos_detectorunc' : "TC",
    'legpos_particlefidunc'  : "TC",
    'legpos_particleunc'  : "TL",
    'legpos_particlefidbinunc': (.45, .615, .63, .9),
#    'legpos_particlefidbin': (.52, .9, .72, .65),
    'legpos_particlefidbin': (.51, .88, .78, .43),
    "legpos_particle"   : (.52, .9, .72, .65),
    'txtsize_covdetector': 0.5,
    'txtsize_covparticle': 1.5,
    'txtangle_covparticle': 35,
    'txtsize_covparticlefidbin': 1.3,
    'txtangle_covparticlefidbin': 35,
    "yaxisuplimitunf": 0.15,
    "yaxismax_particlefidunc" : 0.7,
    "yaxismax_particlefidbinunc" : 1.3,
    "yaxismax_unf" : 1.55,
    "yaxismax_particlebin": 0.004,
    "legpos_particlebinunc" : (.18, .5, .31, .785),
    "legpos_particlefidbinunc": (.4, .615, .58, .9),
    "yaxismax_ratio_fidbin" : 2.5,
}


varList['Jet1_Pt'] = {
    'printname'     : 'Leading jet \\pt (\\GeV)',
    'printnamenodim': 'Leading jet \\pt',
    'mathprintname' : '\\text{Leading jet }\\pt',
    'xaxis'         : 'Leading jet #it{p}_{T} (GeV)',
    'yaxis_particle': 'd#sigma/d(leading jet #it{p}_{T}) (pb)',
    'yaxisfidbin'   : '(1/#sigma_{fid.})d#sigma/d(leading jet  #it{p}_{T}) (1/GeV)',
    'yaxisbin'      : '(1/#sigma_{fid.})d#sigma/d(leading jet  #it{p}_{T}) (pb/GeV)',
    'yaxis_unc'     : 'Relative uncertainty',
    'bins_particle' : [30., 50., 70., 85., 110., 150.], # propuesta (5 bins)
    'bins_detector' : [30., 40., 45., 50., 55., 60., 65., 70., 85., 110., 150.], # propuesta (5 bins) SELECCTIONADA

#    'bins_particle' : [30., 45., 60., 80., 100., 150.],
#    'bins_detector' : [30., 35., 40., 45., 50., 55., 60., 70., 80., 90., 100., 125., 150.],

#    'bins_particle' : [30., 50., 70., 85., 110., 150.],
#    'bins_detector' : [30., 50., 70., 85., 110., 150.],

    'descbinning' : [10., 150.],
    'ndescbins'   : 14,
    'legpos'      : (0.7, 0.55, 0.90, 0.93),
    'legposdesc'  : (0.7, 0.55, 0.90, 0.93),
    'maxdesc'     : 3600,
    "legpos_particle"  : "TC",
    "var_detector": 'min(max(Jet1_Pt, 31.), 149)',
    'var_response': 'Jet1_Pt',
    'var_particle': 'min(max(DressJet1_Pt, 31.), 149)',
    'legpos_detectorunc' : "TL",
    'legpos_particleunc'  : "TL",
    'legpos_particlefidunc'  : "TL",
    'legpos_particlefidbinunc':"TL",
    # 'legpos_particlefidbin': (.66, .9, .8, .6),
    # 'legpos_particlefidbin': (.6, .9, .75, .55),
    #'legpos_particlefidbin': (.63, .88, .78, .53),
    'legpos_particlefidbin': (.51, .88, .78, .43),
    #'legpos_particlefidbin':"TC",
    "legpos_particlebinunc" : "TL",
    'resptxtsize' : 1.5,
    'txtsize_covparticle': 2,
    'txtsize_covdetector': 1.25,
    "txtsize_covparticlefidbin": 1.8,
    'txtangle_covparticle': 45,
    "txtangle_covparticlefidbin": 45,
    "yaxisuplimitunf": 0.20,
    "yaxismax_particlefidbin": 0.02,
    "yaxismax_particlefid" : 1.8,
    #"yaxismax_particlefidbinunc" : 1.2,
    "yaxismax_particlefidbinunc" : 1.3,
    "yaxismax_unf" : 1,
    "yaxismax_particlebinunc": 1.8,
    "yaxismax_particlebin": 0.010,
    "yaxismax_ratio_norm" : 3.5,
#    "yaxismax_ratio_fidbin" : 2.5,
    "yaxismax_ratio_fidbin" : 1.5,
}

varList['Lep1_Pt'] = {
    #'xaxis'       : 'p_{T}(\\ell_{1}) (GeV)',
    'printname'      : 'Leading lepton \\pt (\\GeV)',
    'printnamenodim' : 'Leading lepton \\pt',
    'mathprintname'  :'\\text{Leading lepton }\\pt',
    'xaxis'          : 'Leading lepton #it{p}_{T} (GeV)',
    'yaxis_particle' : 'd#sigma/d(leading lepton #it{p}_{T}) (pb)',
    # 'yaxisfid'       : '(1/#sigma_{fid.})d#sigma/d(leading lepton #it{p}_{T}) (adim.)',
    'yaxisfid'       : '(1/#sigma_{fid.})d#sigma/d(leading lepton #it{p}_{T})',
    'yaxisfidbin'    : '(1/#sigma_{fid.})d#sigma/d(leading lepton #it{p}_{T}) (1/GeV)',
    'yaxisbin'       : 'd#sigma/d(leading lepton #it{p}_{T}) (pb/GeV)',
    # 'yaxis_unc'      : 'Relative uncertainty (adim.)',
    'yaxis_unc'      : 'Relative uncertainty',

    #### TOP-19-003 8 bines
    #'bins_particle' : [25., 50., 60., 70., 80., 90., 105., 125., 150.],
    #'bins_detector' : [25., 35., 50., 55., 60., 65.,  70.,  75.,  80., 85., 90., 95., 100., 105., 115., 125., 150.],
    #'bins_particle' : [25., 50., 60., 70., 80., 90., 105., 150.],
    #'bins_detector' : [25., 35., 50., 55., 60., 65.,  70.,  75.,  80., 85., 90., 95., 100., 105., 115., 125., 150.],

    ### dev reluncs
    # 'bins_particle' : [25., 50., 60., 70., 80., 90., 100., 150.],
    # 'bins_detector' : [25., 35., 50., 55., 60., 65.,  70.,  75.,  80., 85., 90., 95., 100., 105., 115., 125., 150.],

    #'bins_particle' : [25., 40., 50., 60., 70., 80., 90., 100., 250.],
    #'bins_detector' : [25., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 115., 125., 135., 150., 175., 200., 225., 250.],


    ### Nuevo intento 2022-02-23 EL QUE TOI USANDO
    'bins_particle' : [25., 35., 45., 55., 65., 75., 85., 95., 150.],
    'bins_detector' : [25., 35., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 115., 125., 150.],

    ### Nuevo intento 2022-02-24
    #'bins_particle' : [25., 35., 45., 55., 65., 75., 85., 125., 150.],
    #'bins_detector' : [25., 35., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 115., 125., 150.],

    ### Nuevo intento 2022-02-24v2
    #'bins_particle' : [25., 45., 55., 65., 75., 80., 125., 150.],
    #'bins_detector' : [25., 35., 45., 50., 55., 60., 65.,  70.,  75.,  80., 90., 100., 125., 150.],

    ### Prueba
#    'bins_particle' : [25., 35., 45., 55., 65., 75., 85., 95., 150.],
#    'bins_detector' : [25., 35., 45., 55., 65., 75., 85., 95., 150.],

    'descbinning' : [10., 150.],   ### TOP-19-003 8 bines
    #'descbinning' : [10., 250.],
    'ndescbins'   : 15,
    'legpos'      : (0.7, 0.55, 0.90, 0.93),
    'legposdesc'  : (0.7, 0.55, 0.90, 0.93),
    "maxdesc"     : 3300,
    "var_detector": 'min(max(LepGood_pt_corrAll[0], 26.), 149.)',   ### TOP-19-003 8 bines
    #"var_detector": 'min(max(LepGood_pt_corrAll[0], 26.), 249.)',
    'var_response': 'Lep1_Pt',
    'var_particle': 'min(max(GenDressedLepton_pt[iDressSelLep[0]], 26.), 149.)',   ### TOP-19-003 8 bines
    #'var_particle': 'min(max(GenDressedLepton_pt[iDressSelLep[0]], 26.), 249.)',
    'legpos_detectorunc' : "TL",
    'legpos_particlefidunc'  : "TL",
    #"legpos_particle"   : (.18, .3, .32, .05),
    "legpos_particle"      : "TC",
    #'legpos_particlefidbin': "TC",
    #'legpos_particlefidbin': (.18, .3, .32, .05),
    # 'legpos_particlefidbin': (.66, .9, .8, .6),
#    'legpos_particlefidbin': (.63, .88, .78, .53),
    'legpos_particlefidbin': (.51, .88, .78, .43),
    'legpos_particleunc'   : "TL",
    'legpos_particlefidbinunc': "TL",
    'resptxtsize'  : 1.5,
    'txtsize_covdetector': 1.2,
    'txtsize_covparticle': 1.35,
    'txtangle_covparticle': 42.50,
    "txtsize_covparticlefidbin": 1.2,
    "txtangle_covparticlefidbin": 35,
    "yaxisuplimitunf": 0.2,
    "yaxismax_particlefidbin": 0.03,
    "yaxismax_particlefid" : 1.1,
    "yaxismax_particlefidbinunc" : 1.3,
    #"yaxismax_particlefidbinunc" : 0.5,
    "yaxismax_unf" : 2.0,
    "yaxismax_particlebin": 0.0085,
    "legpos_particlebinunc" : "TL",
#    "yaxismax_ratio_fidbin" : 2.5,
   "yaxismax_ratio_fidbin" : 1.08,
}

varList['Fiducial'] = {
    'xaxis'         : 'a.u.',
    'yaxis_particle': 'd#sigma (pb)',
    'bins_particle' : [25., 150.],
    'bins_detector' : [25., 150.],
    #'bins_detector' : [25., 80., 150.],
    "var_detector"  : 'min(LepGood_pt_corrAll[0], 149.)',
    #"var_detector"  : 'min(Jet1_Pt, 149.)',
    'var_response'  : 'Fiducial',
    'var_particle'  : 'min(GenDressedLepton_pt[iDressSelLep[0]], 149.)',
    #'var_particle'  : 'min(DressJet1_Pt, 149.)',
}


varList['Lep1Lep2Jet1_Pz'] = {
    'xaxis'       : '#it{p}_{Z} (#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j}) (GeV)',
    'printname'   : '\\pzvar (\\GeV)',
    'printnamenodim':'\\pzvar',
    'mathprintname': '\\pzvar',
    'yaxis_particle'       : 'd#sigma/d(#it{p}_{Z} (#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (pb)',
    # 'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(#it{p}_{Z} (#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (adim.)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(#it{p}_{Z} (#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j}))',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(#it{p}_{Z} (#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (1/GeV)',
    'yaxisnorm'   : 'd#sigma/d(#it{p}_{Z} (#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (pb/GeV)',
    # 'yaxis_unc'   : 'Relative uncertainty (adim.)',
    'yaxis_unc'   : 'Relative uncertainty',
    'bins_particle'  : [0., 50., 100., 150., 200., 250., 300., 350., 450.],                        # propuesta (8 bins)
    'bins_detector' : [0., 20., 40., 60., 80., 100., 125., 150., 175., 200., 225., 250., 275., 300., 350., 400., 450.], # propuesta (8 bins) SELECCIONADA

#    'bins_particle' : [0., 50., 100., 150., 200., 250., 300., 350., 450.],
#    'bins_detector' : [0., 50., 100., 150., 200., 250., 300., 350., 450.],

    'descbinning' : [0., 450.],
    'ndescbins'   : 18,
    "var_detector": 'min(max(abs(Lep1Lep2Jet1_Pz), 1.), 449.)',
    'var_response': 'Lep1Lep2Jet1_Pz',
    'var_particle': 'min(max(abs(DressLep1Lep2Jet1MET_Pz), 1.), 449.)',
    'legpos'      : (0.51, 0.55, 0.71, 0.93),
    #'legposdesc'  : (0.57, 0.55, 0.78, 0.93),
    'legposdesc'  : (0.65, 0.55, 0.85, 0.93),
    #'legpos_detector':  "TC",
    'legpos_detector' : (.5, .685, .725, .93),
    'legpos_detectoras':"BL",
    "legpos_particleunc" : "TL",
    'legpos_detectorunc' : "TL",
    'legpos_particlefidunc'  : "TL",
    'legpos_particlefidbinunc': "TL",
    # 'legpos_particlefidbin': (.66, .9, .8, .6),
    # 'legpos_particlefidbin': (.6, .9, .75, .55),
    #'legpos_particlefidbin': (.63, .88, .78, .53),
    'legpos_particlefidbin': (.51, .88, .78, .43),
    #'legpos_particlefidbin':"TC",
    "legpos_particle"  : "TC",
    'resptxtsize' : 1.7,
    'txtsize_covparticle': 1.4,
    'txtangle_covparticle': 45,
    'txtsize_covdetector': 1.4,
    "txtsize_covparticlefidbin": 1.2,
    "txtangle_covparticlefidbin": 45,
    "yaxisuplimitunf": 0.1,
    "yaxismax_particlefidbin": 0.005,
    "yaxismax_particlefid" : 1,
    #"yaxismax_particlefidbinunc" : 1,
    "yaxismax_particlefidbinunc" : 1.3,
    "yaxismax_unf" : 2,
    "yaxismax_particlebin" : 0.002,
    "legpos_particlebinunc" : "TC",
#    "yaxismax_ratio_fidbin" : 2.5,
    "yaxismax_ratio_fidbin" : 1.5,
    # "particle_smoothing" : {"jes_BBEC1" : {-1 : "fitOrder=1, symmAfterFit=1",
    #                                        0  : "NOM",
    #                                        1  : "NOM",
    #                                        7  : "NOM",},
    #                         "jes_FlavorQCD" : "fitOrder=1, symmAfterFit=1",
    #                         "jes_RelativeBal" : {-1 : "fitOrder=1, symmAfterFit=1",
    #                                              0  : "NOM",
    #                                              1  : "NOM",},
    #                         "jes_RelativeSample_2018" : "fitOrder=1, symmAfterFit=1",
    #                         "jes_Absolute" : "fitOrder=1, symmAfterFit=1",
    #                         "jes_Absolute_2018" : {-1 : "fitOrder=1, symmAfterFit=1",
    #                                                0  : "NOM",
    #                                                7  : "NOM",},
    #                         "jer_2016" : {-1 : "fitOrder=2, symmAfterFit=1",
    #                                       0  : "NOM",
    #                                       7  : "NOM",},
    #                         "jer_2017" : {-1 : "fitOrder=2, symmAfterFit=1",
    #                                       0  : "NOM",
    #                                       7  : "NOM",},
    #                         "jer_2018" : {-1 : "fitOrder=2, symmAfterFit=1",
    #                                       0  : "NOM",
    #                                       7  : "NOM",},
    #                         "fsr" : {-1 : "fitOrder=1, symmAfterFit=1",
    #                                  0  : "NOM",
    #                                  2  : "fitOrder=2, symmAfterFit=1",
    #                                  7  : "NOM",},
    #                         "colour_rec_cr1" : {-1 : "fitOrder=1, symmAfterFit=1",
    #                                  0  : "NOM",
    #                                  1  : "NOM",
    #                                  2  : "NOM",
    #                                  3  : "NOM",
    #                                  6  : "fitOrder=2, symmAfterFit=1",
    #                                  7  : "NOM",},
    #                         "colour_rec_cr2" : {-1 : "fitOrder=1, symmAfterFit=1",
    #                                  0  : "NOM",
    #                                  1  : "NOM",
    #                                  2  : "NOM",
    #                                  3  : "NOM",
    #                                  4  : "NOM",
    #                                  5  : "fitOrder=2, symmAfterFit=1",
    #                                  6  : "NOM",},
    #                         "colour_rec_erdon" : {-1 : "fitOrder=1, symmAfterFit=1",
    #                                  1  : "NOM",
    #                                  2  : "NOM",
    #                                  3  : "NOM",
    #                                  5  : "NOM",
    #                                  6  : "NOM",
    #                                  7  : "NOM",},
    #                         "ds" : {-1 : "fitOrder=1, symmAfterFit=1",
    #                                  2  : "NOM",
    #                                  3  : "NOM",
    #                                  4  : "NOM",
    #                                  7  : "NOM",},
    # },
}

varList['Lep1Lep2_DPhi'] = {
    'xaxis'         : "#Delta#it{#varphi}(#it{e}^{#pm}, #it{#mu}^{#mp})/#it{#pi}",
    'printname'     : "$\\deltaPhiVar/ \\pi$",
    'printnamenodim': "$\\deltaPhiVar/ \\pi$",
    'mathprintname' : "\\deltaPhiVar/ \\pi",
    'yaxis_particle': "d#sigma/d(#Delta#it{#varphi}(#it{e}^{#pm}, #it{#mu}^{#mp})/#it{#pi}) (pb)",
    'yaxisfidbin'   : '(1/#sigma_{fid.})d#sigma/d(#Delta#it{#varphi}(#it{e}^{#pm}, #it{#mu}^{#mp})/#it{#pi})',
    'yaxisbin'      : 'd#sigma/d(#Delta#it{#varphi}(#it{e}^{#pm}, #it{#mu}^{#mp})/#it{#pi}) (pb)',
    'yaxis_unc'     : 'Relative uncertainty',
    'bins_particle' : [0., .17, .33, .5, .67, .83, 1.0],                        # propuesta (6 bins)
    'bins_detector' : [0., 0.083, .17, 0.25, 0.33, 0.417, .5, 0.583, 0.67, 0.75, 0.83, 0.917, 1.0], # propuesta (6 bins)

    "var_detector"  : 'Lep1Lep2_DPhi',
    'var_response'  : 'Lep1Lep2_DPhi',
    'var_particle'  : 'DressLep1Lep2_DPhi',
    #'var_particle'     : 'abs(TGenDPhiLL)',
    #'legpos'      : (0.82, 0.14, 0.93, 0.47),
    'legpos'      : (0.15, 0.425, 0.35, 0.81),
    #'legposdesc'  : (0.15, 0.425, 0.35, 0.81),
    'legposdesc'  : (0.15, 0.52, 0.35, 0.9215),
    'maxdesc'     : 2000,
    'legpos_detector' : "TL",
    'legpos_particlefid'  : "BR",
    #'legpos_particle'  : "BC",
    "legpos_particle"   : (.43, .255, .63, .04),
    'legpos_detectoras':"TL",
    'legpos_particleas': "TL",
    # 'legpos_particlefidbin': (.65, .45, .8, .1),
#    'legpos_particlefidbin': (.63, .45, .78, .1),
    'legpos_particlefidbin': (.51, .5, .78, .05),
    #'legpos_particlefidbin':"BC",
    'legpos_particlefidbinunc': "TC",
    "legpos_particleunc"  : "TL",
    'resptxtsize' : 2,
    'txtsize_covparticle': 1.75,
    'txtangle_covparticle': 45,
    'txtsize_covdetector': 1.75,
    "txtsize_covparticlefidbin": 1.6,
    'txtangle_covparticlefidbin': 45,
    "equalbinsunf" : True,
    "equalbinsfol" : True,
#    "yaxismax_detector" : 3100,
    "yaxismax_particlefid" : 0.4,
    #"yaxismax_particlefidbinunc" : 0.35,
#    "yaxismax_particlefidbin" : 1.6,
    "yaxismax_particlefidbin" : 1.55,
    "yaxismax_unf" : 1,
    #'legpos_particlebin'  :"BR",
    'legpos_particlebin'        : (.4, .35, .78, .03),
    "yaxismax_ratio_bin"        : 1.5,
    "yaxismax_particlebin"      : 6.5,
    "yaxismax_ratio_fidbin"     : 1.04,
    "yaxismax_particlebinunc"   : 0.1,
    "yaxismax_particlefidbinunc": 1.3,
}

# Profiling things
#systMap = {
#    'btagging_2016'       : True,
#    'btagging_1718'       : True,
#    'mistagging_2016'     : True,
#    'mistagging_1718'     : True,
#    'muonen_2016'         : True,
#    'muonen_2017'         : True,
#    'muonen_2018'         : True,
#    'elecidsf'            : True,
#    'elecrecosf'          : True,
#    'muonidsf_stat_2016'  : True,
#    'muonidsf_stat_2017'  : True,
#    'muonidsf_stat_2018'  : True,
#    'muonidsf_syst'       : True,
#    'muonisosf_stat_2016' : True,
#    'muonisosf_stat_2017' : True,
#    'muonisosf_stat_2018' : True,
#    'muonisosf_syst'      : True,
#    'pileup'              : True,
#    'prefiring_2016'      : True,
#    'prefiring_2017'      : True,
#    'jes'                 : True,
#    'jer_2016'            : True,
#    'jer_2017'            : True,
#    'jer_2018'            : True,
#    'triggereff_2016'     : True,
#    'triggereff_2017'     : True,
#    'triggereff_2018'     : True,
#    'ttbar_matching' : {'tw'      : False,
#                        'ttbar'   : True,
#                        'dy'      : False,
#                        'nonworz' : False,
#                        'vvttv'   : False},
#    'ue'           : {  'tw'      : False,
#                        'ttbar'   : True,
#                        'dy'      : False,
#                        'nonworz' : False,
#                        'vvttv'   : False},
#    'isr_ttbar' : {     'tw'      : False,
#                        'ttbar'   : True,
#                        'dy'      : False,
#                        'nonworz' : False,
#                        'vvttv'   : False},
#    'isr_tw' : {        'tw'       : True,
#                        'ttbar'   : False,
#                        'dy'      : False,
#                        'nonworz' : False,
#                        'vvttv'   : False},
#    'fsr_ttbar' : {     'tw'      : False,
#                        'ttbar'   : True,
#                        'dy'      : False,
#                        'nonworz' : False,
#                        'vvttv'   : False},
#    'fsr_tw' : {        'tw'       : True,
#                        'ttbar'   : False,
#                        'dy'      : False,
#                        'nonworz' : False,
#                        'vvttv'   : False},
#    'ttbar_scales' : {  'tw'      : False,
#                        'ttbar'   : True,
#                        'dy'      : False,
#                        'nonworz' : False,
#                        'vvttv'   : False},
#    'tw_scales' : {    'tw'       : True,
#                        'ttbar'   : False,
#                        'dy'      : False,
#                        'nonworz' : False,
#                        'vvttv'   : False},
#    'topptrew'   : {    'tw'      : False,
#                        'ttbar'   : True,
#                        'dy'      : False,
#                        'nonworz' : False,
#                        'vvttv'   : False},
#    'colour_rec' : {    'tw'      : False,
#                        'ttbar'   : True,
#                        'dy'      : False,
#                        'nonworz' : False,
#                        'vvttv'   : False},
#}


UncsColourMap = {
    'fit'                 : r.kPink-7,
    'btagging_corr'       : r.TColor.GetColor("#b2df8a"),
    'btagging_2016'       : r.TColor.GetColor("#b2df8a"),
    'btagging_2017'       : r.TColor.GetColor("#b2df8a"),
    'btagging_2018'       : r.TColor.GetColor("#b2df8a"),
    'mistagging_corr'     : r.TColor.GetColor("#ff7f00"),
    'mistagging_2016'     : r.TColor.GetColor("#ff7f00"),
    'mistagging_2017'     : r.TColor.GetColor("#ff7f00"),
    'mistagging_2018'     : r.TColor.GetColor("#ff7f00"),
    'muonen_2016'         : r.kYellow-4,
    'muonen_2017'         : r.kYellow-2,
    'muonen_2018'         : r.kYellow-12,
    'elecidsf'            : r.kGray+7,
    'elecrecosf'          : r.kGray+12,
    'muonidsf_stat_2016'  : r.kBlue+5,
    'muonidsf_stat_2017'  : r.kBlue+8,
    'muonidsf_stat_2018'  : r.kBlue+11,
    'muonidsf_syst'       : r.kBlue+17,
    'muonisosf_stat_2016' : r.kGreen+6,
    'muonisosf_stat_2017' : r.kGreen+8,
    'muonisosf_stat_2018' : r.kGreen+10,
    'muonisosf_syst'      : r.kGreen+12,
    'pileup'              : r.TColor.GetColor("#fb9a99"),
    'prefiring_2016'      : r.kGreen+1,
    'prefiring_2017'      : r.kBlue+3,
    'jes'                 : r.kPink+1,
    'jes_hf'              : r.kPink+1,
    'jes_bbec1_2016'      : r.kPink+1,
    'jes_bbec1_2017'      : r.kPink+1,
    'jes_bbec1_2018'      : r.kPink+1,
    'jes_flavorqcd'       : r.kPink+1,
    'jes_relativesample_2022'   : r.kPink+1,
    'jes_relativesample_2022postee'   : r.kPink+1,
    'jes_ec2'             : r.kPink+1,
    'jes_hf_2016'         : r.kPink+1,
    'jes_hf_2017'         : r.kPink+1,
    'jes_hf_2018'         : r.kPink+1,
    'jes_relativebal'     : r.kPink+1,
    'jes_absolute_2016'   : r.kPink+1,
    'jes_absolute_2017'   : r.kPink+1,
    'jes_absolute_2018'   : r.kPink+1,
    'jes_bbec1'           : r.kPink+1,
    'jes_EC2_2016'        : r.kPink+1,
    'jes_ec2_2017'        : r.kPink+1,
    'jes_ec2_2018'        : r.kPink+1,
    'jes_absolute'        : r.kPink+1,
    'jer_2016'            : r.TColor.GetColor("#1f77b4"),
    'jer_2017'            : r.TColor.GetColor("#a6cee3"),
    'jer_2018'            : r.TColor.GetColor("#33a02c"),
    'triggereff_2016'     : r.kOrange-6,
    'triggereff_2017'     : r.kSpring-9,
    'triggereff_2018'     : r.kOrange-9,
    'isr_ttbar'           : r.TColor.GetColor("#cab2d6"),
    'isr_tw'              : r.kAzure-6,
    'fsr'                 : r.kGray+2,
    'fsr_ttbar'           : r.kGray+2,
    'fsr_tw'              : r.TColor.GetColor("#fdbf6f"),
    'topptrew'            : r.kPink-5,
    'ue'                  : r.kTeal-7,
    'ttbar_matching'      : r.kTeal,
    'scales'              : r.kGreen+4,
    'ttbar_scales'        : r.kGreen+4,
    'tw_scales'           : r.kYellow-6,
    'colour_rec'          : r.kViolet-2,
    'ttbar_norm'          : r.kBlue,
    'vvttv_norm'          : r.kAzure-9,
    'nonworz_norm'        : r.kMagenta-4,
    'dy_norm'             : r.kMagenta,
    'mtop'                : r.kMagenta-3,
    'pdfhessian'          : r.kPink-7,
}
UncsColourMap["colour"] = UncsColourMap["colour_rec"]


UncGroupsColourMap = {
    'btag'                : r.TColor.GetColor("#b2df8a"),
    'mc_stat'             : r.kRed,
    'mistag'              : r.kYellow-2,
    'pdf'                 : r.kYellow-12,
    'elec'                : r.kGray+7,
    'lumi'                : r.kGray+12,
    'muon'                : r.kBlue+5,
    'pileup'              : r.TColor.GetColor("#fb9a99"),
    'prefiring'           : r.kGreen+1,
    'jes'                 : r.kPink+1,
    'experimental'        : r.kPink+1,
    'jer'                 : r.TColor.GetColor("#1f77b4"),
    'trigger'             : r.kOrange-6,
    'isr'                 : r.TColor.GetColor("#cab2d6"),
    'fsr'                 : r.kGray+2,
    'modelling'           : r.kTeal-7,
    'toppt'               : r.kPink-5,
    'ds'                  : r.kPink-5,
    'ue'                  : r.kTeal-7,
    'matching'            : r.kTeal,
    'ttbar_scales'        : r.kPink+4,
    'modellingttbar'      : r.kPink+4,
    'tw_scales'           : r.kGreen-4,
    'modellingtw'         : r.kGreen-4,
    'colour'              : r.kViolet-2,
    'ttbar_norm'          : r.kBlue,
    'normalisation'       : r.kBlue,
    'vvttv_norm'          : r.kAzure-9,
    'nonworz_norm'        : r.kMagenta-4,
    'dy_norm'             : r.kMagenta,
    'mtop'                : r.kMagenta-3,
}

#ColorMapList = [
#    r.TColor.GetColor("#a6cee3"),
#    r.TColor.GetColor("#1f77b4"),
#    r.TColor.GetColor("#b2df8a"),
#    r.TColor.GetColor("#33a02c"),
#    r.TColor.GetColor("#e31a1c"),
#    r.TColor.GetColor("#fdbf6f"),
#    r.TColor.GetColor("#ff7f00"),
#    r.TColor.GetColor("#6a3d9a"),
#    r.kOrange-6,
#    r.TColor.GetColor("#cab2d6"),
#    r.kAzure-6,
#    r.kMagenta,
#    r.kGray+2,
#    r.kBlue,
#    r.kRed,
#    r.kGreen+4,
#    r.kYellow-6,
#    r.kTeal-7,
#    r.kViolet-2,
#    r.kPink+1,
#    r.kSpring-9,
#    r.TColor.GetColor("#fb9a99"),
#    r.kGreen+1,
#    r.kBlue+3,
#    r.kPink-5,
#    r.kTeal,
#]


ColourMapForProcesses = {
    "tw"       : 798,
    "ttbar"    : 633,
    "nonworz"  : 413,
    "vvttv"    : 390,
    "dy"       : 852
}


# Plots (llaves en minuscula)
SysNameTranslator = {
    'fit'                 : "Fit",
    'asimov'              : "Asimov",
    'Nominal'             : "Nominal",
    'btag'                : "B-tag.",
    'btagging_corr'       : "B-tag.",
    'btagging_2016'       : "B-tag. (16)",
    'btagging_2017'       : "B-tag. (17)",
    'btagging_2018'       : "B-tag. (18)",
    "lumi_2016"           : "Lumi. uncorr. (16)",
    "lumi_2017"           : "Lumi. uncorr. (17)",
    "lumi_2018"           : "Lumi. uncorr. (18)",
    "lumi_corr"           : "Luminosity corr.",
    "lumi"                : "Luminosity",
    'mistagging_corr'     : "Mistag.",
    'mistagging_2016'     : "Mistag. (16)",
    'mistagging_2017'     : "Mistag. (17)",
    'mistagging_2018'     : "Mistag. (18)",
    'muonen_2016'         : "Muon en. corr. (16)",
    'muonen_2017'         : "Muon en. corr. (17)",
    'muonen_2018'         : "Muon en. corr. (18)",
    'elecidsf'            : "Electron ID eff.",
    'elecrecosf'          : "Electron reco. eff.",
    "elec"                : "Electron effs.",
    'muonidsf_stat_2016'  : "Muon ID eff. (16, stat.)",
    'muonidsf_stat_2017'  : "Muon ID eff. (17, stat.)",
    'muonidsf_stat_2018'  : "Muon ID eff. (18, stat.)",
    'muonidsf_syst'       : "Muon ID eff. (syst.)",
    'muonisosf_stat_2016' : "Muon ISO eff. (16, stat.)",
    'muonisosf_stat_2017' : "Muon ISO eff. (17, stat.)",
    'muonisosf_stat_2018' : "Muon ISO eff. (18, stat.)",
    'muonisosf_syst'      : "Muon ISO eff. (syst.)",
    'pileup'              : "Pile-up",
    'prefiring_2016'      : "Prefiring (16)",
    'prefiring_2017'      : "Prefiring (17)",
    'jes'                 : "JES",
    'jes_hf'              : "JES",
    'jes_bbec1_2016'      : "JES",
    'jes_bbec1_2017'      : "JES",
    'jes_bbec1_2018'      : "JES",
    'jes_flavorqcd'       : "JES",
    'jes_relativesample_2022'   : "JES",
    'jes_relativesample_2022postee'   : "JES",
    'jes_ec2'             : "JES",
    'jes_hf_2016'         : "JES",
    'jes_hf_2017'         : "JES",
    'jes_hf_2018'         : "JES",
    'jes_relativebal'     : "JES",
    'jes_absolute_2016'   : "JES",
    'jes_absolute_2017'   : "JES",
    'jes_absolute_2018'   : "JES",
    'jes_bbec1'           : "JES",
    'jes_ec2_2016'        : "JES",
    'jes_ec2_2017'        : "JES",
    'jes_ec2_2018'        : "JES",
    'jes_absolute'        : "JES",
    'jer'                 : "JER",
    'jer_2016'            : "JER (16)",
    'jer_2017'            : "JER (17)",
    'jer_2018'            : "JER (18)",
    'triggereff_2016'     : "Trigger eff. (16)",
    'triggereff_2017'     : "Trigger eff. (17)",
    'triggereff_2018'     : "Trigger eff. (18)",
    'isr_ttbar'           : "ISR (t#bar{t})",
    'isr_tw'              : "ISR (tW)",
    'fsr_ttbar'           : "FSR (t#bar{t})",
    'fsr_tw'              : "FSR (tW)",
    "fsr"                 : "FSR",
    'topptrew'            : "Top p_{T} rew.",
    'ue'                  : "UE",
    'ttbar_matching'      : "ME/PS matching (t#bar{t})",
    'scales'              : "#mu_{R}/#mu_{F}",
    'ttbar_scales'        : "t#bar{t} #mu_{R}/#mu_{F}",
    'tw_scales'           : "tW #mu_{R}/#mu_{F}",
    'colour_rec'          : "Colour rec.",
    'ttbar_norm'          : "t#bar{t} norm.",
    'vvttv_norm'          : "VVt#bar{t}V norm.",
    'nonworz_norm'        : "Non-W/Z norm.",
    'dy_norm'             : "DY norm.",
    'mtop'                : "Top mass",
    "pdfhessian"          : "PDF + #alpha_{S}",
    "mc_stat"             : "MC stat.",
    "normalisation"       : "Normalisation",
    "modelling"           : "Modelling",
    "modellingtw"         : "Modelling (tW)",
    "modellingttbar"      : "Modelling (t#bar{t})",
    "experimental"        : "Experimental",
}
SysNameTranslator["colour"]   = SysNameTranslator["colour_rec"]
SysNameTranslator["matching"] = SysNameTranslator["ttbar_matching"]
SysNameTranslator["pdf"]      = SysNameTranslator["pdfhessian"]
SysNameTranslator["toppt"]    = SysNameTranslator["topptrew"]
SysNameTranslator["muon"]     = "Muon efficiencies"

# LaTeX puro (tablas)
PrintSysNameTranslator = {
    'asimov'              : "Asimov",
    'Nominal'             : "Nominal",
    #'lumi_corr'           : "Luminosity corr. (16, 17, 18)",
    #'lumi_corr1718'       : "Luminosity corr. (17, 18)",
    #'lumi_2016'           : "Luminosity uncorr. (2016)",
    #'lumi_2017'           : "Luminosity uncorr. (2017)",
    'lumi'           : "Luminosity",
    #'lumi_2018'           : "Luminosity uncorr. (2018)",
    'btagging'            : "B-tagging",
    #'btagging_2016'       : "B-tagging (2016)",
    #'btagging_1718'       : "B-tagging (2017, 2018)",
    #'mistagging'          : "Mistagging",
    #'mistagging_2016'     : "Mistagging (2016)",
    #'mistagging_1718'     : "Mistagging (2017, 2018)",
    "btagging_corr"       : "B-tagging corr.",
    'btagging_2022'       : "B-tagging uncorr. (2022)",
    'btagging_2022PostEE' : "B-tagging uncorr. (2022PostEE)",
    "mistagging_corr"     : "Mistagging corr.",
    'mistagging_2022'     : "Mistagging uncorr. (2022)",
    'mistagging_2022PostEE': "Mistagging uncorr. (2022PostEE)",
    #'muonen_2016'         : "Muon en. corr. (2016)",
    #'muonen_2017'         : "Muon en. corr. (2017)",
    #'muonen_2018'         : "Muon en. corr. (2018)",
    'elecidsf'            : "Electron ID eff.",
    'elecrecosf'          : "Electron reco. eff.",
    'muonidsf_stat_2022'  : "Muon ID eff. (2022, stat.)",
    'muonidsf_stat_2022PostEE'  : "Muon ID eff. (2022PostEE, stat.)",
    'muonidsf_stat'  : "Muon ID eff. (2022PostEE, stat.)",
    'muonidsf_syst'       : "Muon ID eff. (syst.)",
    'muonisosf_stat_2022' : "Muon ISO eff. (2022, stat.)",
    'muonisosf_stat_2022PostEE' : "Muon ISO eff. (2022PostEE, stat.)",
    'muonisosf_stat' : "Muon ISO eff. (2022PostEE, stat.)",
    'muonisosf_syst'      : "Muon ISO eff. (syst.)",
    'pileup'              : "Pile-up",
    #'prefiring_2016'      : "L1 ECAL prefiring (2016)",
    #'prefiring_2017'      : "L1 ECAL prefiring (2017)",
    'jes'                   : "JES",
    "jes_HF"                  : "JES - HF corr.",
    "jes_HF_2022"             : "JES - HF uncorr. (2022)",
    "jes_HF_2022PostEE"             : "JES - HF uncorr. (2022PostEE)",
    "jes_BBEC1"               : "JES - BBEC1 corr.",
    "jes_BBEC1_2022"          : "JES - BBEC1 uncorr. (2022)",
    "jes_BBEC1_2022PostEE"          : "JES - BBEC1 uncorr. (2022PostEE)",
    "jes_FlavorQCD"           : "JES - Flavour QCD",
    "jes_RelativeSample_2022" : "JES - Relative sample (2022)",
    "jes_RelativeSample_2022PostEE" : "JES - Relative sample (2022PostEE)",
    "jes_EC2"                 : "JES - EC2 corr.",
    "jes_EC2_2022"            : "JES - EC2 uncorr. (2022)",
    "jes_EC2_2022PostEE"            : "JES - EC2 uncorr. (2022PostEE)",
    "jes_RelativeBal"         : "JES - RelativeBal",
    "jes_Absolute"            : "JES - Absolute corr.",
    "jes_Absolute_2022"       : "JES - Absolute corr. (2022)",
    "jes_Absolute_2022PostEE"       : "JES - Absolute corr. (2022PostEE)",
    'jer'            : "JER",
    'unclenergy'    : "Uncl. energy",
    #'jer_2017'            : "JER (2017)",
    #'jer_2018'            : "JER (2018)",
    'triggereff'     : "Trigger eff.",
    #'triggereff_2017'     : "Trigger eff. (2017)",
    #'triggereff_2018'     : "Trigger eff. (2018)",
    "pdfhessian"          : "PDF + $\\alpha_{S}$",
    'isr_ttbar'           : "ISR (\\ttbar)",
    'isr_tw'              : "ISR (\\tw)",
    'fsr'                 : "FSR (\\ttbar, \\tw)",
    'fsr_ttbar'           : "FSR (\\ttbar)",
    'fsr_tw'              : "FSR (\\tw)",
    'topptrew'            : "Top \pt rew.",
    'ue'                  : "UE",
    "mtop"                : "$m_{top}$",
    "ds"                  : "tW DR/DS",
    'ttbar_matching'      : "ME/PS matching (\\ttbar)",
    'ttbar_scales'        : "\\ttbar $\\mu_{R}$/$\\mu_{F}$",
    'tw_scales'           : "\\tw $\\mu_{R}$/$\\mu_{F}$",
    'colour_rec'          : "Colour rec.",
    'colour_rec_cr1'      : "Colour rec. (CR1)",
    'colour_rec_cr2'      : "Colour rec. (CR2)",
    'colour_rec_erdon'    : "Colour rec. (MPI with ERD on)",
}


# Details for transcribing to table
ProcessNameTranslator = {
    "ttbar"        : "\\ttbar",
    "nonworz"      : "Non W/Z",
    "dy"           : "Drell-Yan",
    "vvttv"        : "VV and \\ttbar V",
    "tw"           : "tW (signal)",
    "data"         : "Data",
    "total"        : "Postfit unc.",
}
ProcessNameTranslator["tW"]           = ProcessNameTranslator["tw"]
ProcessNameTranslator["t#bar{t}"]     = ProcessNameTranslator["ttbar"]
ProcessNameTranslator["Non-WorZ"]     = ProcessNameTranslator["nonworz"]
ProcessNameTranslator["DY"]           = ProcessNameTranslator["dy"]
ProcessNameTranslator["VVttbarV"]     = ProcessNameTranslator["vvttv"]
ProcessNameTranslator["VV+t#bar{t}V"] = ProcessNameTranslator["vvttv"]


GOFTranslator = {
    "bb4l"        : "PH b\\bar{b}l^{+}l^{-}\\nu\\bar{\\nu}",
    "DR"          : "PH DR + P8",
    "DS"          : "PH DS + P8",
    "Herwig"      : "PH DR + H7",
    "aMC_dr"      : "aMC DR + P8",
    "aMC_dr2"     : "aMC DR2 + P8",
    "aMC_ds"      : "aMC DS + P8",
    "aMC_ds_runn" : "aMC DS dyn. + P8",
}


TableDict = {
    "detector" : {
            "prefix"     : "$d\\sigma/d\\left(",
            "dimensions" : "(pb)",
        },
    "particle" : {
            "prefix"     : "$d\\sigma/d\\left(",
            "dimensions" : "(pb)",
        },
    "particlefidbin" : {
            "prefix"     : "$(1/\\sigma_{\\rm{fid}})d\\sigma/d\\left(",
            "dimensions" : "(1/\\GeV)",
        },
}

#### NOTE: it is important to note that they are all normalisation uncertainties (ahora ya no jiji)
#ProfileSysts = {"dy_norm", "nonworz_norm", "ttbar_norm", "vvttv_norm",
#                "lumi_2016", "lumi_2017", "lumi_2018",
#                "lumi_BBD", "lumi_DB", "lumi_LS", "lumi_BCC", "lumi_GS", "lumi_XY"}

#ModifiedProfileSysts = {"dy_norm", "nonworz_norm", "ttbar_norm", "vvttv_norm",
#                        "lumi_2016", "lumi_2017", "lumi_2018",
#                        "lumi_BBD", "lumi_DB", "lumi_LS", "lumi_BCC", "lumi_GS", "lumi_XY",
#    'btagging_2016'  ,
#    'btagging_1718'  ,
#    'mistagging_2016',
#    'mistagging_1718',
#    'muonen_2016',
#    'muonen_2017',
#    'muonen_2018',
#    'elecidsf'  ,
#    'elecrecosf',
#    'muonidsf_stat_2016' ,
#    'muonidsf_stat_2017' ,
#    'muonidsf_stat_2018' ,
#    'muonidsf_syst'      ,
#    'muonisosf_stat_2016',
#    'muonisosf_stat_2017',
#    'muonisosf_stat_2018',
#    'muonisosf_syst'     ,
#    'pileup'        ,
#    'prefiring_2016',
#    'prefiring_2017',
#    'jes'       ,
#    'jer_2016'  ,
#    'jer_2017'  ,
#    'jer_2018'  ,
#    'triggereff_2016',
#    'triggereff_2017',
#    'triggereff_2018',
#    #'tw_scales',
#    #'ttbar_scales',
#    #'topptrew',
#    #'isr_ttbar',
#    #'isr_tw',
#    #'fsr_ttbar',
#    #'fsr_tw',
#    }

#ModifiedProfileSystsThatAreNotPresentAllYears = {
#    "lumi_2016"           : ["2016"],
#    "lumi_2017"           : ["2017"],
#    "lumi_2018"           : ["2018"],
#    "lumi_BBD"            : ["2016", "2017"],
#    "lumi_DB"             : ["2016", "2017"],
#    "lumi_LS"             : ["2017", "2018"],
#    "lumi_BCC"            : ["2017", "2018"],
#    "lumi_GS"             : ["2016", "2017"],
#    "lumi_XY"             : ["2016", "2017", "2018"],
#    'btagging_2016'       : ["2016"],
#    'btagging_1718'       : ["2017", "2018"],
#    'mistagging_2016'     : ["2016"],
#    'mistagging_1718'     : ["2017", "2018"],
#    'muonen_2016'         : ["2016"],
#    'muonen_2017'         : ["2017"],
#    'muonen_2018'         : ["2018"],
#    'elecidsf'            : ["2016", "2017", "2018"],
#    'elecrecosf'          : ["2016", "2017", "2018"],
#    'muonidsf_stat_2016'  : ["2016"],
#    'muonidsf_stat_2017'  : ["2017"],
#    'muonidsf_stat_2018'  : ["2018"],
#    'muonidsf_syst'       : ["2016", "2017", "2018"],
#    'muonisosf_stat_2016' : ["2016"],
#    'muonisosf_stat_2017' : ["2017"],
#    'muonisosf_stat_2018' : ["2018"],
#    'muonisosf_syst'      : ["2016", "2017", "2018"],
#    'pileup'              : ["2016", "2017", "2018"],
#    'prefiring_2016'      : ["2016"],
#    'prefiring_2017'      : ["2017"],
#    'jes'                 : ["2016", "2017", "2018"],
#    'jer_2016'            : ["2016"],
#    'jer_2017'            : ["2017"],
#    'jer_2018'            : ["2018"],
#    'triggereff_2016'     : ["2016"],
#    'triggereff_2017'     : ["2017"],
#    'triggereff_2018'     : ["2018"],
#    #'tw_scales'           : ["2016", "2017", "2018"],
#    #'ttbar_scales'        : ["2016", "2017", "2018"],
#    #'topptrew'            : ["2016", "2017", "2018"],
#    #'isr_ttbar'           : ["2016", "2017", "2018"],
#    #'isr_tw'              : ["2016", "2017", "2018"],
#    #'fsr_ttbar'           : ["2016", "2017", "2018"],
#    #'fsr_tw'              : ["2016", "2017", "2018"],
#}


#ProfileSystsThatAreNotPresentAllYears = {
#    "lumi_2016" : ["2016"],
#    "lumi_2017" : ["2017"],
#    "lumi_2018" : ["2018"],
#    "lumi_BBD"  : ["2016", "2017"],
#    "lumi_DB"   : ["2016", "2017"],
#    "lumi_LS"   : ["2017", "2018"],
#    "lumi_BCC"  : ["2017", "2018"],
#    "lumi_GS"   : ["2016", "2017"],
#    "lumi_XY"   : ["2016", "2017", "2018"],
#}


ProcessesNames = {"tw", "ttbar", "dy", "nonworz", "vvttv"}


ControlRegionVariableDict = {"3j2t" : "Lep1Lep2Jet1_M_control",
}


#WeightsToChange = {
#    'btagging_2016'      : "bTagWeight_btag",
#    'btagging_1718'      : "bTagWeight_btag",
#    'mistagging_2016'    : "bTagWeight_mistag",
#    'mistagging_1718'    : "bTagWeight_mistag",
#    'elecidsf'           : "ElecIDSF",
#    'elecrecosf'         : "ElecRECOSF",
#    'muonidsf_stat_2016' : "MuonIDSF_stat",
#    'muonidsf_stat_2017' : "MuonIDSF_stat",
#    'muonidsf_stat_2018' : "MuonIDSF_stat",
#    'muonidsf_syst'      : "MuonIDSF_syst",
#    'muonisosf_stat_2016': "MuonISOSF_stat",
#    'muonisosf_stat_2017': "MuonISOSF_stat",
#    'muonisosf_stat_2018': "MuonISOSF_stat",
#    'muonisosf_syst'     : "MuonISOSF_syst",
#    'pileup'             : "puWeight",
#    'prefiring_2016'     : "PrefireWeight",
#    'prefiring_2017'     : "PrefireWeight",
#    'triggereff_2016'    : "TrigSF",
#    'triggereff_2017'    : "TrigSF",
#    'triggereff_2018'    : "TrigSF",
#}

#UncertaintiesToEnvelope = {"ttbar_scales", "tw_scales", "colour"}
#UncertaintiesToEnvelope = {"ttbar_scales", "tw_scales"}
#UncertaintiesToEnvelope = {"colour"}
UncertaintiesToEnvelope = {}

VariatedSamplesUncertaintySources = { "ue", "colour_rec_erdon", "colour_rec_cr1", "colour_rec_cr2", "tw_matching", "mtop", "ds"} # Necesario para splitear los uncertainties en particle bins

nuisanceColours = {
    'Gaussian'          : 1,
    'Poisson'           : 8,
    'AsymmetricGaussian': 9,
    'Unconstrained'     : 39,
    'Unrecognised'      : 2
}

lo_colour = {
    'default': 38,
    'hesse'  : r.kOrange - 3,
    'robust' : r.kGreen + 1
}
hi_colour = {
    'default': 46,
    'hesse'  : r.kBlue,
    'robust' : r.kAzure - 5
}

coloursForDiffWithCombine = {
    'bkg': r.kOrange + 10, # bkg proc
    0    : r.kOrange + 10, # bkg proc
    1    : r.kOrange,
    2    : r.kOrange - 1,
    3    : r.kOrange - 2,
    4    : r.kOrange - 3,
    5    : r.kOrange - 4,
    6    : r.kOrange - 5,
    7    : r.kOrange - 6,
    8    : r.kOrange - 7,
    9    : r.kOrange - 8,
    10   : r.kOrange - 9,
}

coloursForToys = {
    0    : r.kOrange +10,
    1    : r.kOrange,
    2    : r.kGreen,
    3    : r.kGreen + 3,
    4    : r.kAzure,
    5    : r.kViolet,
    6    : r.kViolet + 10,
    7    : r.kPink + 6,
    8    : r.kRed - 4,
    9    : r.kOrange - 8,
}

# individual_list = ['fsr','pileup','dy_norm','colour','jes','matching','ttbar_norm','elec','nonworz_norm','btag','ue','toppt','tw_scales',
#                    'vvttv_norm','mc_stat','isr','ttbar_scales','lumi','pdf','mistag','ds','jer','trigger','prefiring','muon','mtop']

#individual_list = ['ttbar_scales','btag','jes','ttbar_norm','toppt','matching','mtop','elec','pileup','lumi','colour','dy_norm','mc_stat',
#                   'vvttv_norm','tw_scales','ue','isr','jer','nonworz_norm','fsr','ds','pdf','mistag','trigger','prefiring','muon']

individual_list = ["modelling", "experimental", "normalisation", "mc_stat"]
global_list     = ['systematics']

markersdict = {"tw"                         : 2,
               "twds"                     : 4,
               "tru_herwig"                 : 26,
               "twamcatnlo_dr"              : 32,
               "twamcatnlo_dr2"             : 24,
               "twamcatnlo_ds"              : 25,
               "twamcatnlo_ds_runningBW"    : 28,
               "twamcatnlo_ds_is"              : 30,
               "twamcatnlo_ds_is_runningBW" : 5,
                "tru"                : 2,
               "tru_DS"             : 4,
               "tru_herwig"         : 26,
               "tru_aMC_dr"         : 32,
               "tru_aMC_dr2"        : 24,
               "tru_aMC_ds"         : 25,
               "tru_aMC_ds_runn"    : 28,
               "tru_aMC_ds_IS"      : 30,
               "tru_aMC_ds_IS_runn" : 5,}

spacingdict = {"tw"                          : -0.9,
               "twds"                      : -0.6,
               "tru_herwig"                  : -0.3,
               "twamcatnlo_dr"               : +0.0,
               "twamcatnlo_dr2"              : +0.3,
               "twamcatnlo_ds"               : +0.6,
               "twamcatnlo_ds_runningBW"     : +0.9,
               "twamcatnlo_ds_is"               : +0.93,
               "twamcatnlo_ds_is_runningBW"  : +0.97,
               "tru"                : -0.9,
               "tru_DS"             : -0.6,
               "tru_herwig"         : -0.3,
               "tru_aMC_dr"         : +0.0,
               "tru_aMC_dr2"        : +0.3,
               "tru_aMC_ds"         : +0.6,
               "tru_aMC_ds_runn"    : +0.9,
               "tru_aMC_ds_IS"      : +0.93,
               "tru_aMC_ds_IS_runn" : +0.97,}

"""
comparisonColourDict = {"bb4l"                      : 880,
#                        "twttbardr"                 : r.kRed,
                        "tw"                 : r.kRed,
                        "twttbards"                 : r.kGreen+3,
                        "twttbarherwig"             : r.kMagenta,
                        "twttbaramc_dr"             : r.kAzure,
                        "twttbaramc_dr2"            : r.kAzure - 1,
                        "twttbaramc_ds"             : r.kCyan,
                        "twttbaramc_ds_runningBW"   : r.kPink - 9,
                        "twttbaramc_ds_is"          : r.kViolet - 4,
                        "twttbaramc_ds_is_runningBW": r.kViolet - 1}"""
comparisonColourDict = {"tw"                        : r.kRed,
                        "twds"                    : r.kGreen+3,
                        "tru_herwig"                : r.kMagenta,
                        "twamcatnlo_dr"             : r.kAzure,
                        "twamcatnlo_dr2"            : r.kAzure - 1,
                        "twamcatnlo_ds"             : r.kOrange + 3,
                        "twamcatnlo_ds_runningBW"   : r.kPink - 9,
                        "twamcatnlo_ds_is"             : r.kViolet - 4,
                        "twamcatnlo_ds_is_runningBW": r.kViolet - 1,
                        "tru"                : r.kRed,
                        "tru_DS"             : r.kGreen+3,
                        "tru_herwig"         : r.kMagenta,
                        "tru_aMC_dr"         : r.kAzure,
                        "tru_aMC_dr2"        : r.kAzure - 1,
                        "tru_aMC_ds"         : r.kOrange + 3,
                        "tru_aMC_ds_runn"    : r.kPink - 9,
                        "tru_aMC_ds_IS"      : r.kViolet - 4,
                        "tru_aMC_ds_IS_runn" : r.kViolet - 1}