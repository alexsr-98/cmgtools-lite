import ROOT as r
import math, sys
from copy  import deepcopy
from array import array
import varList as vl
import numpy as np
###############################################################################

def quadSum(elements):
    return math.sqrt( sum( [x*x for x in elements]))


def GetMaxUnc(nominal, uncUp, uncDown):
    return max([x*x for x in [nominal - uncUp, nominal - uncDown]])


def GetSymUnc(nominal, uncUp, uncDown):
    return vl.mean([abs(nominal - uncUp), abs(nominal - uncDown)])**2 # Applying the AVERAGE of the uncs.
    #return max    ([abs(nominal - uncUp), abs(nominal - uncDown)])**2 # Applying the MAXIMUM of the uncs.

def histogram_to_vector(hist):
    """Convert TH1 histogram to a 1D NumPy array."""
    return np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)])

def computeCovarianceFromResults(histograms, covnomi):
    # Convert each histogram to a vector
    hist_vectors = [histogram_to_vector(histograms[hist]) for hist in histograms]

    # Concatenate vectors into a single 2D array
    data_matrix = np.vstack(hist_vectors)

    # Compute the covariance matrix
    covariance_matrix = np.cov(data_matrix, rowvar=False)    

    n_bins = covariance_matrix.shape[0]
    for i in range(n_bins):
        for j in range(n_bins):
            covnomi.SetBinContent(i+1, j+1, covariance_matrix[i, j])
    return covnomi


def propagateQuantity(nom, varDict, case = 0):
    tot = 0.
    if case == 0: ### Symmetric case
        for key in varDict:
            if 'Down' in key: continue
            if not 'Up' in key:
                raise RuntimeError('Use case is not supported yet. Syst %s is not supported'%key)

            if key.replace('Up', 'Down') in varDict:
                tot = tot + GetSymUnc(nom, varDict[key], varDict[key.replace('Up', 'Down')])
            else:
                tot = tot + (nom - varDict[key])**2

    elif case > 0:
        for key in varDict: 
            if 'Down' in key: continue
            if not 'Up' in key and key not in vl.varList['Names']['colorSysts']:
                raise RuntimeError('Use case is not supported yet. Syst %s is not supported'%key)

            if key.replace('Up','Down') in varDict:
                if (nom - varDict[key] > 0 and nom - varDict[key.replace('Up', 'Down')] > 0):
                    continue
                elif (nom - varDict[key] < 0 and nom - varDict[key.replace('Up', 'Down')] < 0):
                    tot = tot + GetMaxUnc(nom, varDict[key], varDict[key.replace('Up', 'Down')])
                    #if 'fsr' in key or 'FSR' in key or 'isr' in key or 'ISR' in key:
                        #tot = tot + GetMaxUnc(nom, varDict[key], varDict[key.replace('Up', 'Down')])/2
                    #else:
                        #tot = tot + GetMaxUnc(nom, varDict[key], varDict[key.replace('Up', 'Down')])
                else:
                    if (nom - varDict[key] < 0):
                        #if 'fsr' in key or 'FSR' in key or 'isr' in key or 'ISR' in key:
                            #tot = tot + map(lambda x : x*x, [nom - varDict[key]])[0]/2
                        #else:
                            #tot = tot + map(lambda x : x*x, [nom - varDict[key]])[0]
                        tot = tot + [x*x for x in [nom - varDict[key]]][0]
                    else:
                        #if 'fsr' in key or 'FSR' in key or 'isr' in key or 'ISR' in key:
                            #tot = tot + map(lambda x : x*x, [nom - varDict[key.replace('Up', 'Down')]])[0]/2
                        #else:
                            #tot = tot + map(lambda x : x*x, [nom - varDict[key.replace('Up', 'Down')]])[0]
                        tot = tot + [x*x for x in [nom - varDict[key.replace('Up', 'Down')]]][0]
    else:
        for key in varDict: 
            if 'Down' in key: continue
            if not 'Up' in key and key not in vl.varList['Names']['colorSysts']:
                raise RuntimeError('Use case is not supported yet. Syst %s is not supported'%key)

            if key.replace('Up','Down') in varDict:
                if (nom - varDict[key] < 0 and nom - varDict[key.replace('Up', 'Down')] < 0): 
                    continue
                elif (nom - varDict[key] > 0 and nom - varDict[key.replace('Up', 'Down')] > 0):
                    #if 'fsr' in key or 'FSR' in key or 'isr' in key or 'ISR' in key:
                        #tot = tot + GetMaxUnc(nom, varDict[key], varDict[key.replace('Up', 'Down')])/2
                    #else:
                        #tot = tot + GetMaxUnc(nom, varDict[key], varDict[key.replace('Up', 'Down')])
                    tot = tot + GetMaxUnc(nom, varDict[key], varDict[key.replace('Up', 'Down')])
                else:
                    if (nom - varDict[key] > 0):
                        #if 'fsr' in key or 'FSR' in key or 'isr' in key or 'ISR' in key:
                            #tot = tot + map(lambda x : x*x, [nom - varDict[key]])[0]/2
                        #else:
                            #tot = tot + map(lambda x : x*x, [nom - varDict[key]])[0]
                        tot = tot + [x*x for x in [nom - varDict[key]]][0]
                    else:
                        #if 'fsr' in key or 'FSR' in key or 'isr' in key or 'ISR' in key:
                            #tot = tot + map(lambda x : x*x, [nom - varDict[key.replace('Up', 'Down')]])[0]/2
                        #else:
                            #tot = tot + map(lambda x : x*x, [nom - varDict[key.replace('Up', 'Down')]])[0]
                        tot = tot + [x*x for x in [nom - varDict[key.replace('Up', 'Down')]]][0]
    return math.sqrt(tot)


def propagateHisto(varDict, doSym = False):
    outUp   = deepcopy(varDict[""].Clone(varDict[""].GetName() + 'uncUp'))
    outDown = deepcopy(varDict[""].Clone(varDict[""].GetName() + 'uncDown'))

    for bin in range(1, varDict[""].GetNbinsX() + 1):
        err     = outUp.GetBinError(bin)    # <==  Fit unc. taken here or statistical unc. taken here
        cont    = outUp.GetBinContent(bin)
        tmpDict = dict([(key, histo.GetBinContent(bin)) for (key, histo) in varDict.items() if key != ""])

        envelopesUp = []; envelopesDown = []
        for iU in vl.UncertaintiesToEnvelope:
            tmpuncUp    = 0.
            tmpuncDown  = 0.
            for key,histo in varDict.items():
                if (iU + "_") in key:
                    tmpunc  = tmpDict[key] - cont
                    if tmpunc > tmpuncUp:
                        tmpuncUp    = tmpunc
                    elif tmpunc < tmpuncDown:
                        tmpuncDown  = tmpunc
                    del tmpDict[key]
            envelopesUp.append(tmpuncUp)
            envelopesDown.append(tmpuncDown)

        if doSym:
            envelopes = [vl.mean([envelopesUp[i], abs(envelopesDown[i])]) for i in range(len(envelopesUp)) ]
            totsymunc = quadSum([propagateQuantity(cont, tmpDict, 0), err] + envelopes)
            outUp.SetBinError(bin, totsymunc); outDown.SetBinError(bin, totsymunc)
        else:
            outUp.SetBinError(  bin, quadSum([propagateQuantity(cont, tmpDict, +1), err] + envelopesUp))
            outDown.SetBinError(bin, quadSum([propagateQuantity(cont, tmpDict, -1), err] + envelopesDown))

        #print "bin", bin, "inc", totsymunc, varDict[""].GetBinContent(bin), varDict["jesUp"].GetBinContent(bin) - varDict[""].GetBinContent(bin), varDict["jesDown"].GetBinContent(bin) - varDict[""].GetBinContent(bin)
    return [outUp, outDown]


def maxRelativeError(histo):
    try:
        return max( [histo.GetBinError(bin) / histo.GetBinContent(bin) for bin in range(1, histo.GetNbinsX()+1)])
    except ZeroDivisionError:
        raise RuntimeError('There is (at least) one bin with zero contents')


def relativeErrorHist(histo):
    relErr = histo.Clone( histo.GetName() + 'relErr')
    for bin in range(1, histo.GetNbinsX()+1):
        relErr.SetBinContent(bin, histo.GetBinError(bin) / histo.GetBinContent(bin) )
        relErr.SetBinError(bin, 0)
    return relErr

def MeanUncertaintyHisto(histo):
    try:
        return vl.mean( [histo.GetBinError(bin) / histo.GetBinContent(bin) for bin in range(1, histo.GetNbinsX()+1)])
    except ZeroDivisionError:
        raise RuntimeError('There is (at least) one bin with zero contents')

def maximumHisto(histo1, histo2):
    maxHist = histo1.Clone( histo1.GetName() + 'max' ) 
    for bin in range(1, maxHist.GetNbinsX()+1):
        maxHist.SetBinContent( bin, max(histo1.GetBinContent(bin), histo2.GetBinContent(bin)))
    return maxHist


def getUncList(varDict, doFit = True, doSym = False):
    tmpDict = deepcopy(varDict)
    medList = [("Fit" if doFit else "Statistical", tmpDict[""])]
    outList = []

    #histUp  = [[] for i in range(nom.GetNbinsX())]
    #histDown= [[] for i in range(nom.GetNbinsX())]
    #for var in tmpDict:
        #if var == 'asimov': continue
        #hist    = nom.Clone(nom.GetName() + var)
        #variat  = 0.
        #for bin in range(1, nom.GetNbinsX() + 1):
            #variat = nom.GetBinContent(bin) - tmpDict[var].GetBinContent(bin)
            ##if 'fsr' in var or 'FSR' in var or 'isr' in var or 'ISR' in var: variat/(2**(1/2))
            #if doEnv and var in vl.varList['Names']['colorSysts']:
                #if (variat >= 0):
                    #if abs(variat) >= nom.GetBinContent(bin): histDown[bin - 1].append(nom.GetBinContent(bin))
                    #else:                                     histDown[bin - 1].append(abs(variat))
                #else: histUp[bin - 1].append(abs(variat))
            #else:
                #if abs(variat) >= nom.GetBinContent(bin): hist.SetBinError(bin, nom.GetBinContent(bin))
                #else:                                     hist.SetBinError(bin, abs(variat))
        #if not doEnv or var not in vl.varList['Names']['colorSysts']:
            #medList.append( (var,hist) )

    #if doEnv:
        #finalhistUp     = nom.Clone(nom.GetName() + 'ColorRUp')
        #finalhistDown   = nom.Clone(nom.GetName() + 'ColorRDown')
        #for bin in range(1, nom.GetNbinsX() + 1):
            #finalhistUp.SetBinError(bin, max(histUp[bin - 1] + [0]))
            #finalhistDown.SetBinError(bin, max(histDown[bin - 1] + [0]))
        #medList.append(('ColorRUp', finalhistUp))
        #medList.append(('ColorRDown', finalhistDown))

    if not doSym:
        if len(vl.UncertaintiesToEnvelope):
            for iUE in vl.UncertaintiesToEnvelope:
                histUp   = deepcopy(tmpDict[""].Clone(tmpDict[""].GetName() + iUE + "Up"))
                histDown = deepcopy(tmpDict[""].Clone(tmpDict[""].GetName() + iUE + "Down"))
                init = False
                for var in varDict:
                    if iUE in var:
                        for iB in range(1, tmpDict[var].GetNbinsX() + 1):
                            if not init:
                                histUp.SetBinError(iB, 0)
                                histDown.SetBinError(iB, 0)
                            tmpunc  = tmpDict[var].GetBinError(iB) - tmpDict[""].GetBinError(iB)
                            if tmpunc > histUp.GetBinError(iB):
                                histUp.SetBinError(iB, tmpunc)
                            elif abs(tmpunc) > histDown.GetBinError(iB):
                                if abs(tmpunc) > histDown.GetBinContent(iB):
                                    histDown.SetBinError(iB, histDown.GetBinContent(iB))
                                else:
                                    histDown.SetBinError(iB, abs(tmpunc))

                        if not init: init = True
                        del tmpDict[var]

                if not init:
                    raise RuntimeError("FATAL: an envelope for the " + iUE + " unc. was required, but not a single uncertainty source compatible with such scaffold was found.")

                medList.append( (iUE + "Up", deepcopy(histUp)) )
                medList.append( (iUE + "Down", deepcopy(histDown)) )


        for var in tmpDict:
            if var == '': continue
            hist   = deepcopy(tmpDict[""].Clone(tmpDict[""].GetName() + var))

            for bin in range(1, tmpDict[""].GetNbinsX() + 1):
                variat = tmpDict[""].GetBinContent(bin) - tmpDict[var].GetBinContent(bin)

                if abs(variat) >= tmpDict[""].GetBinContent(bin) and "Down" in var:
                    hist.SetBinError(bin, tmpDict[""].GetBinContent(bin))
                else:
                    hist.SetBinError(bin, abs(variat))

            medList.append( (var, hist) )

        medList.sort(key = lambda x : MeanUncertaintyHisto(x[1]), reverse = True)
        #medList.sort(key = lambda x : maxRelativeError(x[1]), reverse = True)
        #medList.sort(key = lambda x : relativeErrorHist(x[1]).GetBinContent(1), reverse = True)

        #for el in medList:
            #print el

        #sys.exit()

        #### NOTE: this does not support one sided variations, except from statistical ones: ALL SYST. VARS.
        #          SHOULD HAVE UP AND DOWN VARIATIONS!

        for key in medList:
            done = False
            for el in outList:
                if key[0].replace('Up', '').replace('Down', '') == el[0]:
                    done = True
            if done: continue

            if 'Down' in key[0]:
                for key2 in medList:
                    if key2[0] == key[0].replace('Down','Up'):
                        hist = maximumHisto(key[1], key2[1])
                        outList.append( (key[0].replace('Down',''), hist) )
                        break
            elif 'Up' in key[0]:
                for key2 in medList:
                    if key2[0] == key[0].replace('Up','Down'):
                        hist = maximumHisto(key[1], key2[1])
                        outList.append( (key[0].replace('Up', ''), hist) )
            else: # We expect only the fit unc. or statistical uncertainty to arrive here
                outList.append( (key[0].replace('_',''), key[1]) )

    else: #### Symmetric uncertainties
        if len(vl.UncertaintiesToEnvelope):
            for iUE in vl.UncertaintiesToEnvelope:
                histUp   = deepcopy(tmpDict[""].Clone(tmpDict[""].GetName() + iUE))
                histDown = deepcopy(tmpDict[""].Clone(tmpDict[""].GetName() + iUE + "Down"))
                init = False
                for var in varDict:
                    #print var
                    if iUE in var:
                        for iB in range(1, tmpDict[var].GetNbinsX() + 1):
                            if not init:
                                histUp.SetBinError(iB, 0)
                                histDown.SetBinError(iB, 0)
                            tmpunc  = tmpDict[var].GetBinError(iB) - tmpDict[""].GetBinError(iB)
                            if tmpunc > histUp.GetBinError(iB):
                                histUp.SetBinError(iB, tmpunc)
                            elif abs(tmpunc) > histDown.GetBinError(iB):
                                histDown.SetBinError(iB, abs(tmpunc))

                        if not init: init = True
                        del tmpDict[var]

                if not init:
                    raise RuntimeError("FATAL: an envelope for the " + iUE + " unc. was required, but not a single uncertainty source compatible with such scaffold was found.")

                for iB in range(1, tmpDict[""].GetNbinsX() + 1):
                    theaverage = vl.mean([histUp.GetBinError(iB), histDown.GetBinError(iB)])
                    histUp.SetBinError(iB, theaverage if theaverage < histUp.GetBinContent(iB) else histUp.GetBinContent(iB))

                medList.append( (iUE, deepcopy(histUp)) )

        for var in tmpDict:
            if "Up" in var:
                hist    = deepcopy(tmpDict[""].Clone(tmpDict[""].GetName() + var.replace("Up", "")))
                variat  = 0.
                for bin in range(1, tmpDict[""].GetNbinsX() + 1):
                    if var.replace("Up", "Down") in tmpDict:
                        variat = math.sqrt(GetSymUnc(tmpDict[""].GetBinContent(bin), tmpDict[var].GetBinContent(bin), tmpDict[var.replace("Up", "Down")].GetBinContent(bin)))
                    else:
                        raise RuntimeError("FATAL: unexpectedly, a one-sided uncertainty that is not the nominal one has been found, " + var + ", after calculating all the applicable envelopes.")

                    if abs(variat) >= tmpDict[""].GetBinContent(bin): hist.SetBinError(bin, tmpDict[""].GetBinContent(bin))
                    else:                                             hist.SetBinError(bin, abs(variat))

                    #print var, hist.GetBinContent(bin), hist.GetBinError(bin)

                medList.append( (var.replace("Up", ""), hist) )

        medList.sort(key = lambda x : MeanUncertaintyHisto(x[1]), reverse = True)
        #medList.sort(key = lambda x : maxRelativeError(x[1]), reverse = True)
        #medList.sort(key = lambda x : relativeErrorHist(x[1]).GetBinContent(1), reverse = True)

        #for el in medList:
            #print el

        #sys.exit()

        #### NOTE: this does not support one sided variations, except from statistical ones: ALL SYST. VARS.
        #          SHOULD HAVE UP AND DOWN VARIATIONS!

        #for key in medList:
            #done = False
            #for el in outList:
                #if key[0].replace('Up', '').replace('Down', '') == el[0]:
                    #done = True
            #if done: continue

            #if 'Down' in key[0]:
                #for key2 in medList:
                    #if key2[0] == key[0].replace('Down','Up'):
                        #hist = maximumHisto(key[1], key2[1])
                        #outList.append( (key[0].replace('Down',''), hist) )
                        #break
            #elif 'Up' in key[0]:
                #for key2 in medList:
                    #if key2[0] == key[0].replace('Up','Down'):
                        #hist = maximumHisto(key[1], key2[1])
                        #outList.append( (key[0].replace('Up', ''), hist) )
            #else: # We expect only the fit unc. or statistical uncertainty to arrive here
                #outList.append( (key[0].replace('_',''), key[1]) )

        outList = medList

    return [(x[0], relativeErrorHist(x[1])) for x in outList]



def getUncListv2(varDict, doSym = False):
    tmpDict = deepcopy(varDict)
    medList = []
    outList = []

    # doSym = False

    #histUp  = [[] for i in range(nom.GetNbinsX())]
    #histDown= [[] for i in range(nom.GetNbinsX())]
    #for var in tmpDict:
        #if var == 'asimov': continue
        #hist    = nom.Clone(nom.GetName() + var)
        #variat  = 0.
        #for bin in range(1, nom.GetNbinsX() + 1):
            #variat = nom.GetBinContent(bin) - tmpDict[var].GetBinContent(bin)
            ##if 'fsr' in var or 'FSR' in var or 'isr' in var or 'ISR' in var: variat/(2**(1/2))
            #if doEnv and var in vl.varList['Names']['colorSysts']:
                #if (variat >= 0):
                    #if abs(variat) >= nom.GetBinContent(bin): histDown[bin - 1].append(nom.GetBinContent(bin))
                    #else:                                     histDown[bin - 1].append(abs(variat))
                #else: histUp[bin - 1].append(abs(variat))
            #else:
                #if abs(variat) >= nom.GetBinContent(bin): hist.SetBinError(bin, nom.GetBinContent(bin))
                #else:                                     hist.SetBinError(bin, abs(variat))
        #if not doEnv or var not in vl.varList['Names']['colorSysts']:
            #medList.append( (var,hist) )

    #if doEnv:
        #finalhistUp     = nom.Clone(nom.GetName() + 'ColorRUp')
        #finalhistDown   = nom.Clone(nom.GetName() + 'ColorRDown')
        #for bin in range(1, nom.GetNbinsX() + 1):
            #finalhistUp.SetBinError(bin, max(histUp[bin - 1] + [0]))
            #finalhistDown.SetBinError(bin, max(histDown[bin - 1] + [0]))
        #medList.append(('ColorRUp', finalhistUp))
        #medList.append(('ColorRDown', finalhistDown))

    if not doSym:
        for var in tmpDict:
            if var == '': continue
            hist   = deepcopy(tmpDict[""].Clone(tmpDict[""].GetName() + var))

            #for bin in range(1, tmpDict[""].GetNbinsX() + 1):
                #variat = tmpDict[""].GetBinContent(bin) - tmpDict[var].GetBinContent(bin)

                #if abs(variat) >= tmpDict[""].GetBinContent(bin) and "Down" in var:
                    #hist.SetBinError(bin, tmpDict[""].GetBinContent(bin))
                #else:
                    #hist.SetBinError(bin, abs(variat))
            for bin in range(1, tmpDict[""].GetNbinsX() + 1):
                variat = tmpDict[var].GetBinError(bin)
                hist.SetBinError(bin, abs(variat))

            medList.append( (var, hist) )

        medList.sort(key = lambda x : MeanUncertaintyHisto(x[1]), reverse = True)
        #medList.sort(key = lambda x : maxRelativeError(x[1]), reverse = True)
        #medList.sort(key = lambda x : relativeErrorHist(x[1]).GetBinContent(1), reverse = True)

        #for el in medList:
            #print el

        #sys.exit()

        #### NOTE: this does not support one sided variations, except from statistical ones: ALL SYST. VARS.
        #          SHOULD HAVE UP AND DOWN VARIATIONS!

        for key in medList:
            done = False
            for el in outList:
                if key[0].replace('Up', '').replace('Down', '') == el[0]:
                    done = True
            if done: continue

            if 'Down' in key[0]:
                for key2 in medList:
                    if key2[0] == key[0].replace('Down','Up'):
                        hist = maximumHisto(key[1], key2[1])
                        outList.append( (key[0].replace('Down',''), hist) )
                        break
            elif 'Up' in key[0]:
                for key2 in medList:
                    if key2[0] == key[0].replace('Up','Down'):
                        hist = maximumHisto(key[1], key2[1])
                        outList.append( (key[0].replace('Up', ''), hist) )

    else: #### Symmetric uncertainties
        for var in tmpDict:
            # print var
            if "Up" in var:
                hist    = deepcopy(tmpDict[""].Clone(tmpDict[""].GetName() + var.replace("Up", "")))
                variat  = 0.
                for bin in range(1, tmpDict[""].GetNbinsX() + 1):
                    if var.replace("Up", "Down") in tmpDict:
                        variat = math.sqrt(GetSymUnc(tmpDict[""].GetBinContent(bin), tmpDict[var].GetBinContent(bin) + tmpDict[var].GetBinError(bin), tmpDict[var.replace("Up", "Down")].GetBinContent(bin) - tmpDict[var.replace("Up", "Down")].GetBinError(bin)))
                    else:
                        raise RuntimeError("FATAL: unexpectedly, a one-sided uncertainty that is not the nominal one has been found, " + var + ", after calculating all the applicable envelopes.")

                    if abs(variat) >= tmpDict[""].GetBinContent(bin): hist.SetBinError(bin, tmpDict[""].GetBinContent(bin))
                    else:                                             hist.SetBinError(bin, abs(variat))

                    #print var, hist.GetBinContent(bin), hist.GetBinError(bin)

                medList.append( (var.replace("Up", ""), hist) )

        medList.sort(key = lambda x : MeanUncertaintyHisto(x[1]), reverse = True)
        #medList.sort(key = lambda x : maxRelativeError(x[1]), reverse = True)
        #medList.sort(key = lambda x : relativeErrorHist(x[1]).GetBinContent(1), reverse = True)

        #for el in medList:
            #print el

        #sys.exit()

        #### NOTE: this does not support one sided variations, except from statistical ones: ALL SYST. VARS.
        #          SHOULD HAVE UP AND DOWN VARIATIONS!

        #for key in medList:
            #done = False
            #for el in outList:
                #if key[0].replace('Up', '').replace('Down', '') == el[0]:
                    #done = True
            #if done: continue

            #if 'Down' in key[0]:
                #for key2 in medList:
                    #if key2[0] == key[0].replace('Down','Up'):
                        #hist = maximumHisto(key[1], key2[1])
                        #outList.append( (key[0].replace('Down',''), hist) )
                        #break
            #elif 'Up' in key[0]:
                #for key2 in medList:
                    #if key2[0] == key[0].replace('Up','Down'):
                        #hist = maximumHisto(key[1], key2[1])
                        #outList.append( (key[0].replace('Up', ''), hist) )
            #else: # We expect only the fit unc. or statistical uncertainty to arrive here
                #outList.append( (key[0].replace('_',''), key[1]) )

        outList = medList

    return [(x[0], relativeErrorHist(x[1])) for x in outList]


def SetTheStatsUncs(histo):
    ''' Function that sets the uncs. of the given histogram as the expecteds
    from its contents.'''
    for bin in range(1, histo.GetNbinsX() + 1):
        histo.SetBinError(bin, math.sqrt(histo.GetBinContent(bin)))
    
    return


##def SetTheUncsFromHere(histo, hlist, SetStatUncs = False):
##    ''' Function that returns to you a histogram with the quadratic sum
##    of the uncertainties corresponding to the bin uncertainties of the 
##    histograms that are given in a list, in addition to the ones from the
##    given histogram.'''
##    
##    tmperr = 0
##    for i in range(1, histo.GetNbinsX() + 1):
##        if SetStatUncs: tmperr = quadSum([math.sqrt(hel.GetBinContent(bin)) for hel in hlist] + math.sqrt(histo.GetBinContent(bin)))
##        else:           tmperr = quadSum([hel.GetBinError(bin) for hel in hlist] + histo.GetBinError(bin))
##        
##        histo.SetBinError(bin, tmperr)
##    
##    return


##def getCovarianceFromVar(nom, var, name, year = "2016", ty = "detector", doCorr = False):
##    nbins   = nom.GetXaxis().GetNbins()
##    binning = array('f', vl.varList[name]['bins_detector'] if ty == "detector" else vl.varList[name]['bins_particle'])
##    if ty == "detector" and vl.doxsec:
##        thelumi = vl.TotalLumi if year == "run2" else vl.LumiDict[int(year)]
##        scaleval = 1/thelumi/1000 if vl.doxsec else 1
##        var.Scale(scaleval)
##    cov     = r.TH2D(var.GetName().replace("data_", "").replace(name+"_", ''), '', nbins, binning, nbins, binning)
##    for x in range(nbins):
##        for y in range(nbins):
##            bin = cov.GetBin(x + 1, y + 1)
##            cov.SetBinContent(bin, (nom.GetBinContent(x + 1) - var.GetBinContent(x + 1)) * (nom.GetBinContent(y + 1) - var.GetBinContent(y + 1)))
##
##    return cov
##
##def getCovarianceFromVarv2(nom, var, name):
##    nbins   = nom.GetXaxis().GetNbins()
##    binning = array('f', vl.varList[name]['bins_particle'])
##    cov     = r.TH2D(var.GetName().replace("data_", "").replace(name+"_", ''), '', nbins, binning, nbins, binning)
##    for x in range(nbins):
##        for y in range(nbins):
##            bin = cov.GetBin(x + 1, y + 1)
##            cov.SetBinContent(bin, (nom.GetBinContent(x + 1) - var.GetBinContent(x + 1)) * (nom.GetBinContent(y + 1) - var.GetBinContent(y + 1)))
##
##    return cov

def drawTheRelUncPlot(listWithHistos, thedict, thePlot, yaxismax = "auto", doSym = False, doFit = False):
    # Calculate the order
    uncList = getUncList(thedict, doFit, doSym)

    #for el in uncList:
    #    print el

    #sys.exit()
    # Calculate the total uncertainty histogram and the total systematic one
    incmax = []; incsyst = []
    maxinctot = 0
    hincmax   = deepcopy(thedict[""].Clone('hincmax'))
    hincmax.SetLineColor(r.kBlack)
    hincmax.SetLineWidth( 2 )
    hincmax.SetFillColorAlpha(r.kBlue, 0.)

    hincsyst  = deepcopy(thedict[""].Clone('hincsyst'))
    hincsyst.SetLineColor(r.kBlack)
    hincsyst.SetLineWidth( 2 )
    hincsyst.SetLineStyle( 3 )
    hincsyst.SetFillColorAlpha(r.kBlue, 0.)

    for bin in range(1, listWithHistos[0].GetNbinsX() + 1):
        if listWithHistos[1].GetBinError(bin) > listWithHistos[0].GetBinContent(bin):
            incmax.append(max([listWithHistos[0].GetBinError(bin), listWithHistos[0].GetBinContent(bin)]))
        else:
            incmax.append(max([listWithHistos[0].GetBinError(bin), listWithHistos[1].GetBinError(bin)]))
        try:
            if math.sqrt(listWithHistos[1].GetBinError(bin)**2 - thedict[""].GetBinError(bin)**2) > listWithHistos[0].GetBinContent(bin):
                incsyst.append(max([math.sqrt(listWithHistos[0].GetBinError(bin)**2 - thedict[""].GetBinError(bin)**2),
                                    listWithHistos[0].GetBinContent(bin)]))
            else:
                incsyst.append(max([math.sqrt(listWithHistos[0].GetBinError(bin)**2 - thedict[""].GetBinError(bin)**2),
                                    math.sqrt(listWithHistos[1].GetBinError(bin)**2 - thedict[""].GetBinError(bin)**2)]))
        except ValueError:
            #print listWithHistos[0].GetBinContent(bin), listWithHistos[0].GetBinError(bin), listWithHistos[1].GetBinContent(bin), listWithHistos[1].GetBinError(bin), thedict[""].GetBinContent(bin), thedict[""].GetBinError(bin), listWithHistos[1].GetBinError(bin)**2 - thedict[""].GetBinError(bin)**2
            #raise RuntimeError("FATAL: the provided nominal/fit only unc. is larger that the total sum of all of them for bin {b}".format(b = bin))
            print("TEMPORAL FIX")
            incsyst.append(max([thedict[""].GetBinError(bin),
                                listWithHistos[0].GetBinContent(bin)]))

        hincmax.SetBinContent(bin, incmax[bin-1] / hincmax.GetBinContent(bin))
        hincmax.SetBinError(bin, 0.)

        if (hincmax.GetBinContent(bin) > maxinctot): maxinctot = hincmax.GetBinContent(bin)

        hincsyst.SetBinContent(bin, incsyst[bin-1] / hincsyst.GetBinContent(bin))
        hincsyst.SetBinError(bin, 0.)

        #print "bin:", bin, "incmax:", incmax[-1]/thedict[""].GetBinContent(bin)*100, "incsyst:", incsyst[-1]/thedict[""].GetBinContent(bin)*100, "theinc:", max([listWithHistos[0].GetBinError(bin), listWithHistos[1].GetBinError(bin)])/thedict[""].GetBinContent(bin)*100

    #sys.exit()

    # Set maximum of the y axis
    if yaxismax == "auto":
        if (maxinctot >= 0.9):
            if maxinctot >= 5:
                hincmax.GetYaxis().SetRangeUser(0, 5)
            else:
                hincmax.GetYaxis().SetRangeUser(0, maxinctot + 0.1)

        else:
            hincmax.GetYaxis().SetRangeUser(0, 0.9)
    else:
        hincmax.GetYaxis().SetRangeUser(0, yaxismax)


    # Begin drawing; first the total, systematic and statistical lines
    thePlot.addHisto(hincmax,  'hist',      'Total',      'L')
    hincstat = None
    if not vl.onlyTotal:
        thePlot.addHisto(hincsyst, 'hist,same', 'Systematic' if not doFit else "Not fit", 'L')

        for i in range(len(uncList)):
            if "Stat" in uncList[i][0]: # The zero element is the name, the 1st is the histogram
                uncList[i][1].SetLineColor(r.kBlack)
                uncList[i][1].SetLineStyle( 2 )
                uncList[i][1].SetLineWidth( 2 )

                hincstat = deepcopy(uncList[i][1].Clone("hincstat"))

                uncList[i][1].SetFillColorAlpha(r.kBlue, 0.)
                thePlot.addHisto(uncList[i][1], 'hist,same', "Statistical", 'L')
            if "Fit" in uncList[i][0]:
                uncList[i][1].SetLineColor(r.kBlack)
                uncList[i][1].SetLineStyle( 2 )

                hincstat = deepcopy(uncList[i][1].Clone("hincfit"))

                uncList[i][1].SetFillColorAlpha(r.kBlue, 0.)
                thePlot.addHisto(uncList[i][1], 'hist,same', vl.SysNameTranslator[uncList[i][0].lower()], 'L')

        # And now, the rest of the lines are up to vl.nuncs of the most dominant systematic sources
        iS = 0
        plottedsysts = 0
        maxnumofuncs = vl.nuncs

        while plottedsysts < maxnumofuncs and plottedsysts < len(uncList):
            if "Stat" in uncList[iS][0] or "Fit" in uncList[iS][0]:
                #print "holi", iS
                iS += 1
                if len(uncList) == 1: break
                continue

            #if "lumi" in uncList[iS][0].lower():
            #    uncList[iS][1].SetLineColor(r.kBlack)
            #    uncList[iS][1].SetLineStyle( 4 )
            else:
                #uncList[iS][1].SetLineColor( vl.ColorMapList[iS] )
                uncList[iS][1].SetLineColor( vl.UncsColourMap[uncList[iS][0].lower().replace("resp_", "")] )
                uncList[iS][1].SetLineWidth( 2 )

            uncList[iS][1].SetFillColorAlpha(r.kBlue, 0.)
            thePlot.addHisto(uncList[iS][1], 'H,same', vl.SysNameTranslator[uncList[iS][0].lower().replace("resp_", "")] + (" (resp.)" if "resp" in uncList[iS][0].lower() else ""), 'L')
            plottedsysts += 1

            #print iS, uncList[iS][0]
            iS += 1


    return (uncList, hincstat, hincsyst, hincmax)


# For the differential with combine 
def drawTheRelUncPlotv2(listWithHistos, thedict, thePlot, yaxismax = "auto", doSym = False):
    # Calculate the order
    subdict = deepcopy(thedict)
    del subdict["statUp"]; del subdict["statDown"]; del subdict["systUp"]; del subdict["systDown"]
    if "totalUp" in subdict:
        del subdict["totalUp"]; del subdict["totalDown"]
    #print subdict

    uncList = getUncListv2(subdict, doSym)

    #for el in uncList:
    #    print el

    #sys.exit()
    # Calculate the total uncertainty histogram and the total systematic one
    maxinctot = 0.
    hincmax   = deepcopy(thedict[""].Clone('hincmax'))
    hincmax.SetLineColor(r.kBlack)
    hincmax.SetLineWidth( 2 )
    hincmax.SetFillColorAlpha(r.kBlue, 0.)

    hincsyst  = deepcopy(thedict[""].Clone('hincsyst'))
    hincsyst.SetLineColor(r.kBlack)
    hincsyst.SetLineWidth( 2 )
    hincsyst.SetLineStyle( 3 )
    hincsyst.SetFillColorAlpha(r.kBlue, 0.)

    hincstat = deepcopy(thedict[""].Clone('hincstat'))
    hincstat.SetLineColor(r.kBlack)
    hincstat.SetLineStyle( 2 )
    hincstat.SetFillColorAlpha(r.kBlue, 0.)

    for bin in range(1, listWithHistos[0].GetNbinsX() + 1):
        hincmax.SetBinContent(bin, max([listWithHistos[0].GetBinError(bin),
                                        listWithHistos[1].GetBinError(bin)])/thedict[""].GetBinContent(bin))

        if (hincmax.GetBinContent(bin) > maxinctot): maxinctot = hincmax.GetBinContent(bin)

        hincsyst.SetBinContent(bin, max([thedict["systUp"].GetBinError(bin),
                                         thedict["systDown"].GetBinError(bin)]) / thedict["systUp"].GetBinContent(bin))

        hincstat.SetBinContent(bin, max([thedict["statUp"].GetBinError(bin),
                                         thedict["statDown"].GetBinError(bin)]) / thedict["statUp"].GetBinContent(bin))
        hincmax.SetBinError(bin,  0.)
        hincsyst.SetBinError(bin, 0.)
        hincstat.SetBinError(bin, 0.)

        #print "bin:", bin, "incmax:", hincmax.GetBinContent(bin)*100, "incsyst:", hincsyst.GetBinContent(bin)*100, #listWithHistos[0].GetBinError(bin), listWithHistos[1].GetBinError(bin)

    #sys.exit()

    # Set maximum of the y axis
    if yaxismax == "auto":
        if (maxinctot >= 0.9):
            if maxinctot >= 5:
                hincmax.GetYaxis().SetRangeUser(0, 5)
            else:
                hincmax.GetYaxis().SetRangeUser(0, maxinctot + 0.1)

        else:
            hincmax.GetYaxis().SetRangeUser(0, 0.9)
    else:
        hincmax.GetYaxis().SetRangeUser(0, yaxismax)


    # Begin drawing; first the total, systematic and statistical lines
    thePlot.addHisto(hincmax,  'hist',      'Total',      'L')
    if not vl.onlyTotal:
        thePlot.addHisto(hincsyst, 'hist,same', 'Systematic',  'L')
        thePlot.addHisto(hincstat, 'hist,same', 'Statistical', 'L')

        # And now, the rest of the lines are up to vl.nuncs of the most dominant systematic sources
        iS = 0
        plottedsysts = 0
        maxnumofuncs = vl.nuncs

        while plottedsysts < maxnumofuncs and plottedsysts < len(uncList):
            if "lumi" in uncList[iS][0].lower():
                uncList[iS][1].SetLineColor(r.kBlack)
                uncList[iS][1].SetLineStyle( 4 )
            else:
                #uncList[iS][1].SetLineColor( vl.ColorMapList[iS] )
                uncList[iS][1].SetLineColor( vl.UncGroupsColourMap[uncList[iS][0].lower()] )
                uncList[iS][1].SetLineWidth( 2 )

            uncList[iS][1].SetFillColorAlpha(r.kBlue, 0.)
            thePlot.addHisto(uncList[iS][1], 'H,same', vl.SysNameTranslator[uncList[iS][0].lower()], 'L')
            #thePlot.addHisto(uncList[iS][1], 'H,same', uncList[iS][0], 'L')
            plottedsysts += 1

            #print iS, uncList[iS][0]
            iS += 1


    return (uncList, hincstat, hincsyst, hincmax)
