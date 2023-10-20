import ROOT as r
import warnings as wr
import sys, os, math, argparse
from copy import deepcopy
from multiprocessing import Pool
from array import array

sys.path.append('{cmsswpath}/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/differential/'.format(cmsswpath = os.environ['CMSSW_BASE']))
import errorPropagator as ep
import beautifulUnfoldingPlots as bp
import varList as vl
#import getLaTeXtable as tex

r.gROOT.SetBatch(True)

#===================================
def getXsecForSys(syst, thedict):
    data = deepcopy(thedict["data_obs"][""].Clone('data_%s'%syst))

    ttbar     = thedict['ttbar'  ][""]
    dy        = thedict['dy'     ][""]
    nonworz   = thedict['nonworz'][""]
    vvttbarv  = thedict['vvttv'  ][""]
    
    if syst in thedict['ttbar']:
        ttbar    = thedict['ttbar'][syst]

    if syst in thedict['dy']:
        dy       = thedict['dy'][syst]

    if syst in thedict['nonworz']:
        nonworz  = thedict['nonworz'][syst]

    if syst in thedict['vvttv']:
        vvttbarv = thedict['vvttv'][syst]

    data.Add(ttbar   , -1 )
    data.Add(dy      , -1 )
    data.Add(nonworz , -1 )
    data.Add(vvttbarv, -1 )
    
    return data


def SubtractBackgroundToData(pathtofile):
    tfile   = r.TFile.Open(pathtofile + '/forExtr.root', "read")

    histodict = {}; unclist = [ "" ]
    for key in tfile.GetListOfKeys():
        tmpnam = key.GetName()
        # In Asimov the stat. unc. is not included in the data_obs histogram
        if "data_obs" in key.GetName():
            ep.SetTheStatsUncs(tfile.Get(tmpnam))
        if "Up" not in tmpnam and "Down" not in tmpnam: # It is the nominal value of a process
            histodict[tmpnam.replace("x_", "")] = {}
            histodict[tmpnam.replace("x_", "")][""] = deepcopy(tfile.Get(tmpnam).Clone(tmpnam.replace("x_", "") + "_"))
        else:
            tmpproc = tmpnam.replace("x_", "").split("_")[0]
            tmpunc  = ("_".join(tmpnam.replace("x_", "").split("_")[1:]))
            #print(tmpnam, tmpproc, tmpunc, tmpvar)

            if tmpunc not in unclist:
                unclist.append(tmpunc)

            if tmpproc not in histodict:
                histodict[tmpproc] = {}

            histodict[tmpproc][tmpunc]= deepcopy(tfile.Get(tmpnam).Clone(tmpproc + "_" + tmpunc))

    tfile.Close(); del tfile

    out = r.TFile.Open(pathtofile + '/detectorsignal_bs.root', 'recreate')
    variations = {}
    for syst in unclist:
        variations[syst]  = deepcopy(getXsecForSys(syst, histodict))
        variations[syst].Write()
    out.Close()

    return variations


def ExtractSignalAndPlotIt(tsk):
    inpath, iY, iV = tsk

    dictofvars = SubtractBackgroundToData(inpath + "/" + iY + "/" + iV)
    #print dictofvars

    if not os.path.isdir(inpath + "/" + iY + "/" + iV + "/detectorplots/"):
        os.system("mkdir -p " + inpath + "/" + iY + "/detectorplots/")

    PlotDetectorLevelResults(inpath, iY, iV, dictofvars)
    return

def PlotDetectorLevelResults(inpath, iY, iV, thedict):
    #scaleval = 1
    thelumi = vl.TotalLumi if iY == "run3" else vl.LumiDict[iY]
    scaleval = 1/thelumi/1000 if vl.doxsec else 1

    for key in thedict:
        thedict[key].Scale(scaleval)

    #### 1) Plot result plot
    nominal_withErrors  = ep.propagateHisto(thedict, doSym = vl.doSym) # Up and down variations with all the systematics

    #for iB in range(1, nominal_withErrors[0].GetNbinsX() + 1):
        #print nominal_withErrors[0].GetBinContent(iB), nominal_withErrors[0].GetBinError(iB)
        #print nominal_withErrors[1].GetBinContent(iB), nominal_withErrors[1].GetBinError(iB)
    
    statOnlyList = [deepcopy(thedict[""]),  deepcopy(thedict[""])]
    #for iB in range(1, thedict[""].GetNbinsX() + 1):
    #    thedict[""].SetBinError(iB, 1e-5)

    #sys.exit()


    plot                = bp.beautifulUnfPlot('{var}_detector'.format(var = iV), iV)
    plot.doRatio        = True
    plot.doFit          = False
    plot.plotspath      = inpath + "/" + iY + "/detectorplots/"
    plot.displayedLumi  = vl.TotalLumi if iY == "run3" else vl.LumiDict[iY]
    plot.doPreliminary  = vl.doPre

    thedict[""].SetMarkerStyle(r.kFullCircle)
    thedict[""].GetXaxis().SetNdivisions(505, True)
    thedict[""].SetLineWidth(0)

    nominal_withErrors[0].SetFillColorAlpha(r.kOrange, 1.0)
    nominal_withErrors[0].SetLineColor(0)
    nominal_withErrors[0].SetLineWidth(0)
    nominal_withErrors[0].SetFillStyle(1001)

    statOnlyList[0].SetFillColorAlpha(r.kGray + 1, 1.0)
    statOnlyList[0].SetLineColor(0)
    statOnlyList[0].SetFillStyle(1001)
    #if iV != "Fiducial":
    #    tex.saveLaTeXfromhisto(thedict[""], iV, path = inpath + "/" + iY + "/tables", errhisto = nominal_withErrors[0], ty = "detector")

    if "legpos_detector" in vl.varList[iV]: legloc = vl.varList[iV]["legpos_detector"]
    else:                                   legloc = "TR"

    print('> Adding generated distribution with used software and others.')
    if not os.path.isfile(inpath + "/" + iY + "/" + iV + "/detector.root"):
        raise RuntimeError('The rootfile with the detector level information does not exist')
    
    tmptfile = r.TFile.Open(inpath + "/" + iY + "/" + iV + "/detector.root")
    tru                         = vl.giveMeOneComparison(tmptfile, "tw", scaleval, iV)
    twds                        = vl.giveMeOneComparison(tmptfile, "twds", scaleval, iV)
    #twttbarherwig               = vl.giveMeOneComparison(tmptfile, "twttbarherwig", scaleval, iV)
    twamc_dr               = vl.giveMeOneComparison(tmptfile, "twamcatnlo_dr", scaleval, iV)
    twamc_dr2              = vl.giveMeOneComparison(tmptfile, "twamcatnlo_dr2", scaleval, iV)
    twamc_ds               = vl.giveMeOneComparison(tmptfile, "twamcatnlo_ds", scaleval, iV)
    twamc_ds_runningBW     = vl.giveMeOneComparison(tmptfile, "twamcatnlo_ds_runningBW", scaleval, iV)
    twamc_ds_is            = vl.giveMeOneComparison(tmptfile, "twamcatnlo_ds_is", scaleval, iV)
    twamc_ds_is_runningBW  = vl.giveMeOneComparison(tmptfile, "twamcatnlo_ds_is_runningBW", scaleval, iV)
    """
    tru                         = vl.giveMeOneComparison(tmptfile, "bb4l", scaleval, iV)
    twttbardr                   = vl.giveMeOneComparison(tmptfile, "twttbardr", scaleval, iV)
    twttbards                   = vl.giveMeOneComparison(tmptfile, "twttbards", scaleval, iV)
    twttbarherwig               = vl.giveMeOneComparison(tmptfile, "twttbarherwig", scaleval, iV)
    twttbaramc_dr               = vl.giveMeOneComparison(tmptfile, "twttbaramc_dr", scaleval, iV)
    twttbaramc_dr2              = vl.giveMeOneComparison(tmptfile, "twttbaramc_dr2", scaleval, iV)
    twttbaramc_ds               = vl.giveMeOneComparison(tmptfile, "twttbaramc_ds", scaleval, iV)
    twttbaramc_ds_runningBW     = vl.giveMeOneComparison(tmptfile, "twttbaramc_ds_runningBW", scaleval, iV)
    twttbaramc_ds_is            = vl.giveMeOneComparison(tmptfile, "twttbaramc_ds_is", scaleval, iV)
    twttbaramc_ds_is_runningBW  = vl.giveMeOneComparison(tmptfile, "twttbaramc_ds_is_runningBW", scaleval, iV)
    """
    tmptfile.Close()


    themaxs = []
#    for el in [tru, twttbardr, twttbards, twttbarherwig, twttbaramc_dr, twttbaramc_dr2, twttbaramc_ds, twttbaramc_ds_runningBW, 
#               thedict[""], nominal_withErrors[0], nominal_withErrors[1]]:
    # Add here the additional samples
    for el in [tru, twds, twamc_dr, twamc_dr2, twamc_ds, twamc_ds_runningBW, twamc_ds_is, twamc_ds_is_runningBW,
               thedict[""], nominal_withErrors[0], nominal_withErrors[1]]:
        themaxs.append(vl.getAConservativeMaximum(el))
    tmpval = max(themaxs)

    nominal_withErrors[0].SetMaximum(tmpval)
    nominal_withErrors[1].SetMaximum(tmpval)
    thedict[""].SetMaximum(tmpval)

    if "yaxismax_detector" in vl.varList[iV]:
        plot.yaxisuplimit = vl.varList[iV]["yaxismax_detector"]

    plot.addHisto(nominal_withErrors,      'A2',     'Total unc.',                     'F', "total")
    plot.addHisto(statOnlyList,            '2,same', 'Stat unc.',                      'F', "stat")
    plot.addHisto(tru,                     'P,same', 'tW PH + P8',                     'P', 'mc')
    plot.addHisto(twds,                    'P,same', 'tW DS PH + P8',       'P', 'mc')
    #plot.addHisto(twttbarherwig,           'P,same', 'tW DR + t#bar{t} PH + H7',       'P', 'mc')
    plot.addHisto(twamc_dr,                'P,same', 'tW DR aMC + P8',      'P', 'mc')
    plot.addHisto(twamc_dr2,               'P,same', 'tW DR2 aMC + P8',     'P', 'mc')
    plot.addHisto(twamc_ds,                'P,same', 'tW DS aMC + P8',      'P', 'mc')
    plot.addHisto(twamc_ds_runningBW,      'P,same', 'tW DS dyn. aMC + P8', 'P', 'mc')
    plot.addHisto(twamc_ds_is,                 'P,same', 'tW DS aMC + P8',      'P', 'mc')
    plot.addHisto(twamc_ds_is_runningBW,   'P,same', 'tW DS dyn. aMC + P8', 'P', 'mc')
    
    """
    plot.addHisto(tru,                     'P,same', 'b#bar{b}l^{+}#nu l^{-}#nu PH + P8','P', 'mc')
    plot.addHisto(twttbardr,               'P,same', 'tW DR + t#bar{t} PH + P8',       'P', 'mc')
    plot.addHisto(twttbards,               'P,same', 'tW DS + t#bar{t} PH + P8',       'P', 'mc')
    plot.addHisto(twttbarherwig,           'P,same', 'tW DR + t#bar{t} PH + H7',       'P', 'mc')
    plot.addHisto(twttbaramc_dr,           'P,same', 'tW DR + t#bar{t} aMC + P8',      'P', 'mc')
    plot.addHisto(twttbaramc_dr2,          'P,same', 'tW DR2 + t#bar{t} aMC + P8',     'P', 'mc')
    plot.addHisto(twttbaramc_ds,           'P,same', 'tW DS + t#bar{t} aMC + P8',      'P', 'mc')
    plot.addHisto(twttbaramc_ds_runningBW, 'P,same', 'tW DS dyn. + t#bar{t} aMC + P8', 'P', 'mc')
    """
    
    plot.addHisto(thedict[""], 'P,E,same{s}'.format(s = ",X0" if "equalbinsunf" in vl.varList[iV] else ""), vl.labellegend, 'PE', 'data', redrawaxis = True)
    #plot.saveCanvas(legloc)
    plot.saveCanvasv2(legloc)
    del plot

    #for iB in range(1, nominal_withErrors[0].GetNbinsX() + 1):
    #    print(thedict[""].GetBinContent(iB), tru.GetBinContent(iB))
    
    #### 2) Plot relative uncertainty plot
    plot2               = bp.beautifulUnfPlot('{var}uncs_detector'.format(var = iV), iV)
    plot2.doFit         = False
    plot2.doPreliminary = vl.doPre
    plot2.plotspath     = inpath + "/" + iY + "/detectorplots/"
    plot2.displayedLumi = vl.TotalLumi if iY == "run3" else vl.LumiDict[iY]

    yaxismax_detectorunc = 1
    if "yaxismax_detectorunc" in vl.varList[iV]:
        yaxismax_detectorunc = vl.varList[iV]["yaxismax_detectorunc"]

    uncListorig, hincstat, hincsyst, hincmax = ep.drawTheRelUncPlot(nominal_withErrors, thedict, plot2, yaxismax_detectorunc, doSym = vl.doSym)
    #uncListorig, hincstat, hincsyst, hincmax = ep.drawTheRelUncPlotv2(nominal_withErrors, thedict, plot, yaxismax_detectorunc)

    if "legpos_detectorunc" in vl.varList[iV]: unclegpos = vl.varList[iV]["legpos_detectorunc"]
    else:                                      unclegpos = "TR"

    plot2.saveCanvas(unclegpos) # Cuidado porque el png no sale bien pero el canvas rootfile si

    ######out2 = r.TFile.Open(inpath + "/" + iY + "/" + iV + "/detectorsignal_bs.root", 'update')
    ######nom0 = deepcopy(nominal_withErrors[0].Clone("nom0"))
    ######nom1 = deepcopy(nominal_withErrors[1].Clone("nom1"))
    ######nom0.Write()
    ######nom1.Write()
    ######hincmax.Write()
    ######hincsyst.Write()
    ######out2.Close(); del out2
    return



if __name__=="__main__":
    parser = argparse.ArgumentParser(usage = "python nanoAOD_checker.py [options]", description = "Checker tool for the outputs of nanoAOD production (NOT postprocessing)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inpath',    '-i', metavar = 'inpath',     dest = "inpath",   required = False, default = "./temp/differential")
    parser.add_argument('--year',      '-y', metavar = 'year',       dest = "year",     required = False, default = "all")
    parser.add_argument('--variable',  '-v', metavar = 'variable',   dest = "variable", required = False, default = "all")
    parser.add_argument('--extraArgs', '-e', metavar = 'extra',      dest = "extra",    required = False, default = "")
    parser.add_argument('--nthreads',  '-j', metavar = 'nthreads',   dest = "nthreads", required = False, default = 0, type = int)
    parser.add_argument('--pretend',   '-p', action  = "store_true", dest = "pretend",  required = False, default = False)


    args     = parser.parse_args()
    year     = args.year
    extra    = args.extra
    nthreads = args.nthreads
    pretend  = args.pretend
    inpath   = args.inpath
    variable = args.variable

    #print("\n===== Unfolding procedures: Response matrices & ROOT files production =====")
    #print("> Setting binning, paths, and other details...")

    #print "\n> Drawing matrices and writing ROOT file (old one will be overwritten!)."


    #### First, find the tasks
    tasks = []
    theyears = []
    presentyears = next(os.walk(inpath))[1]

    if "2022" in presentyears:
        theyears.append("2022")
    if "2022PostEE" in presentyears:
        theyears.append("2022PostEE")
    if "run3" in presentyears:
        theyears.append("run3")

    if year.lower() != "all" and year in presentyears:
        theyears = [ year ]
    elif year.lower() != "all":
        raise RuntimeError("FATAL: the year requested is not in the provided input folder.")

    for iY in theyears:
        thevars = next(os.walk(inpath + "/" + iY))[1]

        if variable.lower() != "all" and variable in thevars:
            thevars = [ variable ]
        elif variable.lower() != "all":
            raise RuntimeError("FATAL: the variable requested is not in the provided input folder.")

        for iV in thevars:
            if any( [el in iV for el in vl.vetolist] ): continue
            tasks.append( (inpath, iY, iV) )


    if nthreads > 1:
        pool = Pool(nthreads)
        pool.map(ExtractSignalAndPlotIt, tasks)
        pool.close()
        pool.join()
    else:
        for tsk in tasks:
            ExtractSignalAndPlotIt(tsk)
