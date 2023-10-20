import ROOT  as r
import os, sys, argparse, math, json
import warnings as wr
from array import array
from copy  import deepcopy
from multiprocessing import Pool
import CombineHarvester.CombineTools.plotting as plot
import CombineHarvester.CombineTools.combine.rounding as rounding

sys.path.append('{cmsswpath}/src/CMGTools/TTHAnalysis/python/plotter/tw-run2/differential/'.format(cmsswpath = os.environ['CMSSW_BASE']))
#import varList as vl

r.PyConfig.IgnoreCommandLineOptions = True
r.TH1.AddDirectory(0)

comm1 = "combineTool.py -M Impacts -d {incard} --doInitialFit --robustFit 1 {ncores} {asimov} {extra} -m 1 -n {prefix} --out {outdir} --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000 --cminDefaultMinimizerStrategy 0" # --cminDefaultMinimizerType Minuit
comm2 = "combineTool.py -M Impacts -d {incard} --robustFit 1 --doFits {ncores} {asimov} {extra} -m 1 -n {prefix} --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000 --cminDefaultMinimizerStrategy 0" # --cminDefaultMinimizerType Minuit
comm3 = "combineTool.py -M Impacts -d {incard} -o impacts{prefix}.json {ncores} {asimov} {extra} -m 1 -n {prefix} --robustFit 1 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000 --cminDefaultMinimizerStrategy 0" # --cminDefaultMinimizerType Minuit
comm4 = "plotImpacts.py -i impacts{prefix}.json -o impacts{prefix} {bl}"


def confirm(message = "Do you wish to continue?"):
    """
    Ask user to enter y(es) or n(o) (case-insensitive).
    :return: True if the answer is Y.
    :rtype: bool
    """
    answer = ""
    while answer not in ["y", "n", "yes", "no"]:
        answer = raw_input(message + " [Y/N]\n").lower()
    return answer[0] == "y"


def checkNfits(tsk):
    thepath, thename = tsk
    failed = False

    thef = r.TFile(thepath + "/" + thename, "READ")
    if thef.limit.GetEntries() < 3:
        failed = True
    thef.Close(); del thef
    return (failed, thename)


def reviewThePreviousStepsFiles(thefolder, nth, verbose):
    thelist  = os.listdir(thefolder)
    thetasks = []
    theres   = []
    for el in thelist:
        if "higgsCombine_" in el:
            thetasks.append( (thefolder, el) )

    if nth > 1:
        pool = Pool(nth)
        theres = pool.map(checkNfits, thetasks)
        pool.close()
        pool.join()
    else:
        for tsk in thetasks:
            theres.append(checkNfits(tsk))

    if any(theres[:][0]):
        print("\n#######################################\nERROR!!\n#######################################\n\nOne or more fits to calculate the impacts have failed. In particular, the following parameters (either nuisances or POI) have one or more failed fits.")
        for el in theres:
            if el[0]: print("\t- " + "_".join(el[1].replace("higgsCombine_", "").replace("paramFit_", "").replace("initialFit_", "").replace(".MultiDimFit.mH1.root", "").split("_")[2:]))

        print("")
        if not confirm("Do you want to, in any case, create the impacts? Please, take into account that the nuisance or POI that has one or more failed fits will NOT be present in the impacts' plot."):
            sys.exit()
    elif verbose:
        print("\t- No failed fits!")

    return




def makeImpacts(task):
    inpath, year, region, physModel, ncores, pretend, verbose, extra, doobs, doBlind = task

    print('\n> Creating impacts for region(s)', region, 'from year', year, '\n')
    
    impactsoutpath = inpath + "/" + iY + "/" + ("Obs" * doObs) + "Impacts_" + region

    if "," not in region:
        physicsModel = 'text2workspace.py -m 125 {infile} -o {outfile}'.format(infile  = inpath + "/" + iY + "/{r}/cuts-tw-{r}.txt".format(r = region),
                                                                               outfile = inpath + "/" + iY + "/{r}/cuts-tw-{r}_ws.root".format(r = region),)
        if verbose:
            print("Text2Workspace command:", physicsModel, "\n")

        if not pretend:
            if os.path.isfile(inpath + "/" + iY + "/{r}/cuts-tw-{r}_ws.root".format(r = region)):
                if verbose:
                    print("    - Erasing old workspace...")
                os.system("rm " + inpath + "/" + iY + "/{r}/cuts-tw-{r}_ws.root".format(r = region))
            outstat = os.system(physicsModel)
            if outstat:
                raise RuntimeError("FATAL: text2workspace.py failed to execute for year {y} and regions {r}.".format(y = year, r = region))

    firstcomm = comm1.format(y      = year,
                             ncores = ("--parallel " + str(ncores)) if ncores else "",
                             asimov = ("" if doobs else "-t -1 --setParameters r=1") if physModel == "1POI" else ("--freezeParameters ttbar_norm" if doobs else "-t -1 --setParameters r_tw=1,r_ttbar=1 --freezeParameters ttbar_norm"),
                             incard = "../combcard_{r}.root".format(r = region.replace(",", "")) if "," in region else "../{r}/cuts-tw-{r}_ws.root".format(r = region),
                             outdir = "./",
                             extra  = extra,
                             prefix = year + "_" + region.replace(",", ""),
    )

    if verbose:
        print("\nFirst command:", firstcomm, "\n")

    if not pretend:
        outstat = os.system("cd " + impactsoutpath + "; " + firstcomm + "; cd -")
        if outstat:
            raise RuntimeError("FATAL: first command failed to execute for region {r} of year {y}.".format(v = region,
                                                                                                           y = year))

    secondcomm = comm2.format(y      = year,
                             ncores = ("--parallel " + str(ncores)) if ncores else "",
                             asimov = ("" if doobs else "-t -1 --setParameters r=1") if physModel == "1POI" else ("--freezeParameters ttbar_norm" if doobs else "-t -1 --setParameters r_tw=1,r_ttbar=1 --freezeParameters ttbar_norm"),
                             incard = "../combcard_{r}.root".format(r = region.replace(",", "")) if "," in region else "../{r}/cuts-tw-{r}_ws.root".format(r = region),
                             outdir = "./",
                             extra  = extra,
                             prefix = year + "_" + region.replace(",", ""),
    )

    if verbose:
        print("\nSecond command:", secondcomm, "\n")
    
    if not pretend:
        outstat = os.system("cd " + impactsoutpath + "; " + secondcomm + "; cd -")
        if outstat:
            raise RuntimeError("FATAL: second command failed to execute for region {r} of year {y}.".format(v = region,
                                                                                                            y = year))
    
    thirdcomm = comm3.format(y      = year,
                             ncores = ("--parallel " + str(ncores)) if ncores else "",
                             asimov = ("" if doobs else "-t -1 --setParameters r=1") if physModel == "1POI" else ("--freezeParameters ttbar_norm" if doobs else "-t -1 --setParameters r_tw=1,r_ttbar=1 --freezeParameters ttbar_norm"),
                             incard = "../combcard_{r}.root".format(r = region.replace(",", "")) if "," in region else "../{r}/cuts-tw-{r}_ws.root".format(r = region),
                             outdir = "./",
                             extra  = extra,
                             prefix = year + "_" + region.replace(",", ""),
    )

    if verbose and not pretend:
        print("\nChecking whether any fit might have failed.")

    if not pretend:
        reviewThePreviousStepsFiles(impactsoutpath, ncores if ncores else 1, verbose)


    if verbose:
        print("\nThird command:", thirdcomm, "\n")
    
    if not pretend:
        outstat = os.system("cd " + impactsoutpath + "; " + thirdcomm + "; cd -")
        if outstat:
            raise RuntimeError("FATAL: third command failed to execute for region {r} of year {y}.".format(v = region,
                                                                                                           y = year))

    fourthcomm = comm4.format(y      = year,
                              ncores = ("--parallel " + str(ncores)) if ncores else "",
                              asimov = ("" if doobs else "-t -1 --setParameters r=1") if physModel == "1POI" else ("--freezeParameters ttbar_norm" if doobs else "-t -1 --setParameters r_tw=1,r_ttbar=1 --freezeParameters ttbar_norm"),
                              incard = "../combcard_{r}.root".format(r = region.replace(",", "")) if "," in region else "../{r}/cuts-tw-{r}_ws.root".format(r = region),
                              outdir = "./",
                              bl     = "--blind" if doBlind else "",
                              prefix = year + "_" + region.replace(",", ""),
    )

    if verbose:
        print("Fourth command:", fourthcomm, "\n")

    if not pretend:
        outstat = os.system("cd " + impactsoutpath + "; " + fourthcomm + "; cd -")
        if outstat:
            raise RuntimeError("FATAL: fourth command failed to execute for region {r} of year {y}.".format(v = region,
                                                                                                            y = year))
                           
    
    print('\n> Region(s)', region, "' impacts produced.\n")
    return



if __name__ == '__main__':
#    vl.SetUpWarnings()
    r.gROOT.SetBatch(True)
    print("===== Fitting procedure with some uncertainty profiling\n")
    parser = argparse.ArgumentParser(usage = "python nanoAOD_checker.py [options]", description = "Checker tool for the outputs of nanoAOD production (NOT postprocessing)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inpath',     '-i', metavar = 'inpath',     dest = "inpath",   required = False, default = "./temp/differential/")
    parser.add_argument('--year',       '-y', metavar = 'year',       dest = "year",     required = False, default = "all")
    parser.add_argument('--region',    '-r', metavar = 'region',        dest = "region",  required = False, default = "1j1t")
    parser.add_argument('--physModel',    '-pM', metavar = 'physModel',     dest = "physModel",   required = False, default = "1POI")
    parser.add_argument('--extraArgs',  '-e', metavar = 'extra',      dest = "extra",    required = False, default = "")
    parser.add_argument('--nthreads',   '-j', metavar = 'nthreads',   dest = "nthreads", required = False, default = 0, type = int)
    parser.add_argument('--pretend',    '-p', action  = "store_true", dest = "pretend",  required = False, default = False)
    parser.add_argument('--verbose',    '-V', action  = "store_true", dest = "verbose",  required = False, default = False)
    parser.add_argument('--doObserved', '-O', action  = "store_true", dest = "doobs",    required = False, default = False)
    parser.add_argument('--blindSignalStrength','-b',action="store_true",dest="blindmu", required = False, default = False)


    args     = parser.parse_args()
    year     = args.year
    nthreads = args.nthreads
    pretend  = args.pretend
    inpath   = args.inpath
    verbose  = args.verbose
    region   = args.region
    physModel = args.physModel
    extra    = args.extra
    doObs    = args.doobs
    doBlind  = args.blindmu


    tasks = []
    theyears = []
    presentyears = next(os.walk(inpath))[1]

    if "2016" in presentyears:
        theyears.append("2016")
    if "2017" in presentyears:
        theyears.append("2017")
    if "2018" in presentyears:
        theyears.append("2018")
    if "run2" in presentyears:
        theyears.append("run2")

    if year.lower() != "all" and year in presentyears:
        theyears = [ year ]
    elif year.lower() != "all":
        raise RuntimeError("FATAL: the year requested is not in the provided input folder.")

    if region.lower() == "all":
        raise RuntimeError("FATAL: you must provide a region, or various of them!")
        
    for iY in theyears:
        if not os.path.isdir(inpath + "/" + iY + "/" + "/" + ("Obs" * doObs) + "Impacts_" + region):
            os.system("mkdir -p " + inpath + "/" + iY + "/" + ("Obs" * doObs) + "Impacts_" + region)
        tasks.append( (inpath, iY, region, physModel, nthreads, pretend, verbose, extra, doObs, doBlind) )

    #print tasks
    for task in tasks:
        print("\nProcessing " + str(task) + "\n")
        makeImpacts(task)
