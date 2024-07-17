from copy import deepcopy
import os, sys, argparse

sys.path.append('{cmsswpath}/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/differential/'.format(cmsswpath = os.environ['CMSSW_BASE']))
import varList as vl
from multiprocessing import Pool

#basecommand = 'combineTool.py -M MultiDimFit {algosettings} --setParameters {setpars} --rMin 0 --rMax 2 --floatOtherPOIs=1 -m 125 --split-points 1 --saveInactivePOI 1 {parallel} {queue} {extra} --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000'
#basecommand = 'combineTool.py -M MultiDimFit {algosettings} --setParameters {setpars}  --split-points 1 --saveInactivePOI 1 {parallel} {queue} {extra} --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000'
#basecommand = 'combineTool.py -M MultiDimFit {algosettings} --setParameters {setpars} --split-points 1 --floatOtherPOIs=1 --saveInactivePOI 1 {parallel} {queue} {extra} --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000 --setParameterRanges {minmaxlist} --squareDistPoiStep --autoRange 2 -m 125'

# basecommand = 'combine -M FitDiagnostics -m 125 --setParameters {setpars} --saveWorkspace --saveShapes --saveWithUncertainties --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000 --setParameterRanges {minmaxlist} {extra}'
#basecommand = 'combine -M FitDiagnostics -m 125 --setParameters {setpars} --saveWorkspace --saveShapes --saveWithUncertainties --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000 {extra}'
basecommand = 'combine -M FitDiagnostics -m 125 --setParameters {setpars} --saveWorkspace --saveShapes --saveWithUncertainties --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000 --robustHesse 1 {extra}'

userstorage = "/pool/phedex/userstorage/asoto/Proyectos/twRun3/"

#individual_scaff = {
#    #### Statistical
#    'mc_stat' : [],
#
#    #### Systematic
#    # Experimental
#    'jes': [
#        "jes_HF",
#        "jes_HF_2016",
#        "jes_HF_2017",
#        "jes_HF_2018",
#        "jes_BBEC1",
#        "jes_BBEC1_2016",
#        "jes_BBEC1_2017",
#        "jes_BBEC1_2018",
#        "jes_FlavorQCD",
#        "jes_RelativeSample_2016",
#        "jes_RelativeSample_2017",
#        "jes_RelativeSample_2018",
#        "jes_EC2",
#        "jes_EC2_2016",
#        "jes_EC2_2017",
#        "jes_EC2_2018",
#        "jes_RelativeBal",
#        "jes_Absolute",
#        "jes_Absolute_2016",
#        "jes_Absolute_2017",
#        "jes_Absolute_2018",
#        "unclenergy",
#    ],
#    'jer': [
#        "jer_2016",
#        "jer_2017",
#        "jer_2018",
#    ],
#    'trigger': [
#        "triggereff_2016",
#        "triggereff_2017",
#        "triggereff_2018",
#    ],
#    'pileup': [
#        "pileup",
#    ],
#    'elec': [
#        "elecidsf",
#        "elecrecosf",
#    ],
#    'muon': [
#        #"muonen_2016",
#        "muonen_2017",
#        "muonen_2018",
#        "muonidsf_stat_2016",
#        "muonidsf_stat_2017",
#        "muonidsf_stat_2018",
#        "muonidsf_syst",
#        "muonisosf_stat_2016",
#        "muonisosf_stat_2017",
#        "muonisosf_stat_2018",
#        "muonisosf_syst",
#    ],
#    'btag': [
#        "btagging_2016",
#        "btagging_2017",
#        "btagging_2018",
#        "btagging_corr",
#    ],
#    'mistag': [
#        "mistagging_2016",
#        "mistagging_2017",
#        "mistagging_2018",
#        "mistagging_corr",
#    ],
#    'lumi': [
#        "lumi",
#    ],
#
#    # Normalisation
#    'ttbar_norm' : [
#        "ttbar_norm",
#    ],
#    'nonworz_norm' : [
#        "nonworz_norm",
#    ],
#    'dy_norm' : [
#        "dy_norm",
#    ],
#    'vvttv_norm' : [
#        "vvttv_norm",
#    ],
#
#    # Modelling
#    'pdf' : [
#        "pdfhessian",
#    ],
#    'matching' : [
#        "ttbar_matching",
#    ],
#    'ttbar_scales' : [
#        "ttbar_scales",
#    ],
#    'tw_scales' : [
#        "tw_scales",
#    ],
#    'isr' : [
#        "isr_ttbar",
#        "isr_tw",
#    ],
#    'fsr' : [
#        "fsr",
#    ],
#    'colour' : [
#        "colour_rec_erdon",
#        "colour_rec_cr1",
#        "colour_rec_cr2",
#    ],
#    'ue' : [
#        "ue",
#    ],
#    'toppt' : [
#        "topptrew",
#    ],
#    'mtop' : [
#        "mtop",
#    ],
#    'ds' : [
#        "ds",
#    ],
#}

individual_scaff = {
    #### Statistical
    'mc_stat' : [],

    #### Systematic
    # Experimental
    'experimental': [
       "jes_HF",
       "jes_HF_2022",
       "jes_HF_2022PostEE",
       "jes_BBEC1",
       "jes_BBEC1_2022",
       "jes_BBEC1_2022PostEE",
       "jes_FlavorQCD",
       "jes_FlavorQCD_bottom",
       "jes_FlavorQCD_charm",
       "jes_FlavorQCD_light",
       "jes_RelativeSample_2022",
       "jes_RelativeSample_2022PostEE",
       "jes_EC2",
       "jes_EC2_2022",
       "jes_EC2_2022PostEE",
       "jes_RelativeBal",
       "jes_Absolute",
       "jes_Absolute_2022",
       "jes_Absolute_2022PostEE",
       "jer_2022",
       "jer_2022PostEE",
       "unclenergy_2022",
       "unclenergy_2022PostEE",
       "triggereff_2022",
       "triggereff_2022PostEE",
       "pileup",
       "elecidsf",
       "elecrecosf",
       "elecscale_2022",
       "elecscale_2022PostEE",
       "elecsmear_2022",
       "elecsmear_2022PostEE",
       #"muonidsf_stat_2022",
       "muonidsf_stat_2022PostEE",
       "muonidsf_syst",
       "muonisosf_stat_2022",
       "muonisosf_stat_2022PostEE",
       "muonisosf_syst",
       "btagging",
       "btagging_2022",
       "btagging_2022PostEE",
       "btagging_corr",
       "mistagging_2022",
       "mistagging_2022PostEE",
       "mistagging_corr",
       "mistagging",
        "lumi",
    ],

    # Normalisation
    'normalisation' : [
        "ttbar_norm",
        "nonworz_norm",
        "dy_norm",
        "vvttv_norm",
    ],

    # Modelling
    'modelling' : [
       "pdfhessian",
       "alphaS",
       "ttbar_matching",
       "ttbar_scales",
       "ttbar_scales_muR",
       "ttbar_scales_muF",
       "tw_scales",
       "tw_scales_muR",
       "tw_scales_muF",
       "isr_ttbar",
       "isr_tw",
       "fsr",
       "colour_rec_erdon",
       "colour_rec_cr1",
       "colour_rec_cr2",
       "ue",
       "topptrew",
       "mtop",
       "mtop3",
       "ds"
    ],
}


global_scaff = {
   'systematics': [
       "jes_HF",
       "jes_HF_2022",
       "jes_HF_2022PostEE",
       "jes_BBEC1",
       "jes_BBEC1_2022",
       "jes_BBEC1_2022PostEE",
       "jes_FlavorQCD",
       "jes_FlavorQCD_bottom",
       "jes_FlavorQCD_charm",
       "jes_FlavorQCD_light",
       "jes_RelativeSample_2022",
       "jes_RelativeSample_2022PostEE",
       "jes_EC2",
       "jes_EC2_2022",
       "jes_EC2_2022PostEE",
       "jes_RelativeBal",
       "jes_Absolute",
       "jes_Absolute_2022",
       "jes_Absolute_2022PostEE",
       "jer_2022",
       "jer_2022PostEE",
       "unclenergy_2022",
       "unclenergy_2022PostEE",
       "triggereff_2022",
       "triggereff_2022PostEE",
       "pileup",
       "elecidsf",
       "elecrecosf",
       "elecscale_2022",
       "elecscale_2022PostEE",
       "elecsmear_2022",
       "elecsmear_2022PostEE",
       #"muonidsf_stat_2022",
       "muonidsf_stat_2022PostEE",
       "muonidsf_syst",
       "muonisosf_stat_2022",
       "muonisosf_stat_2022PostEE",
       "muonisosf_syst",
       "btagging",
       "btagging_2022",
       "btagging_2022PostEE",
       "btagging_corr",
       "mistagging_2022",
       "mistagging_2022PostEE",
       "mistagging_corr",
       "mistagging",
        "lumi",
        "ttbar_norm",
        "nonworz_norm",
        "dy_norm",
        "vvttv_norm",
       "pdfhessian",
       "alphaS",
       "ttbar_matching",
       "ttbar_scales",
       "ttbar_scales_muR",
       "ttbar_scales_muF",
       "tw_scales",
       "tw_scales_muR",
       "tw_scales_muF",
       "isr_ttbar",
       "isr_tw",
       "fsr",
       "colour_rec_erdon",
       "colour_rec_cr1",
       "colour_rec_cr2",
       "ue",
       "topptrew",
       "mtop",
       "mtop3",
       "ds"    
   ],
}


def getBestFitParsValue(thefile):
    outdict = {}
    with open(thefile, "r") as f:
        for line in f:
            splits = line.strip().replace(" ", "").split(":")
            outdict[splits[0]] = splits[1:]
    return outdict


def doMeThisFitPlease(task):
    outpath, thecard, tmpfreeze, thegoodvalsdict, pois, gr, redo, pretend = task
    #print tmpfreeze
    if not os.path.isfile(outpath + "/higgsCombine{g}.FitDiagnostics.mH125.root".format(g = gr)) or redo:
        thecomm = basecommand.format(setpars     = ",".join([el + "=1" for el in pois] + [el + "=" + str(thegoodvalsdict[el][0]) for el in tmpfreeze if el in thegoodvalsdict]),
                                    extra        = '-n {g} {card} --freezeParameters {c} --out {o}'.format(card = thecard,
                                                                                                            g    = gr,
                                                                                                            c    = ",".join(tmpfreeze),
                                                                                                            o    = outpath),
                                    minmaxlist   = ":".join([el + "=0,3" for el in pois]),
                                    )
        print("\nCommand:", thecomm)
        if not pretend:
            os.system(thecomm)
            os.system("mv higgsCombine{u}.FitDiagnostics.mH125.root ".format(u = gr) + outpath + "/")

    return


def calculateRelativeUncertainties(task):
    inpath, iY, iV, nps, nth, redo = task
    
    bins_detector = vl.varList[iV]['bins_detector']
    ndetectorbins = len(bins_detector) - 1
    bins_particle = vl.varList[iV]['bins_particle']
    nparticlebins = len(bins_particle) - 1
    
    reluncpath = inpath + "/"+ iY + "/" + iV + "/sigextr_fit_combine/relativeuncs_fd"
    for iP in ["", "/individual", "/global"]:
        if not os.path.isdir(reluncpath + iP):
            os.system("mkdir -p " + reluncpath + iP)

    doBatch = False
    userstdirrw = userstorage.replace("phedex", "phedexrw").replace("ciencias", "cienciasrw") + "/" + inpath + iY + "/" + iV + "/relativeuncs_fd"
    userstdir   = userstdirrw.replace("phedexrw", "phedex").replace("cienciasrw", "ciencias")

    #print userstdirrw

    if not os.path.isdir(userstdirrw):
        os.system("mkdir -p " + userstdirrw)
        os.system("cp " + inpath + "/"+ iY + "/" + iV + "/sigextr_fit_combine/fit_output.root " + userstdirrw + "/")
        thecard = userstdirrw + "/fit_output.root"
    elif os.path.isfile(userstdir + "/fit_output.root"):
        thecard = userstdir + "/fit_output.root"
    else:
        thecard = userstdirrw + "/fit_output.root"

    thecard = inpath + "/"+ iY + "/" + iV + "/sigextr_fit_combine/fit_output.root"
    
    #print thecard
    #sys.exit()

    #thecard = "../../fit_output.root"
    POIs    = ["r_tW_{iB}".format(iB = i) for i in range(nparticlebins)]
    indGroup = deepcopy(individual_scaff)
    gloGroup = deepcopy(global_scaff)

    for iB in range(ndetectorbins):
        indGroup["mc_stat"].append(    "prop_bincard_" + iV + "_bin" + str(iB))
        gloGroup["systematics"].append("prop_bincard_" + iV + "_bin" + str(iB))


    thegoodvalsdict = getBestFitParsValue(inpath + "/"+ iY + "/" + iV + "/sigextr_fit_combine/fitdiagnostics/fitpars.txt")


    subtasks = []
    indfreeze = []
    for gr in vl.individual_list:
        indfreeze += indGroup[gr]
        subtasks.append( (reluncpath + "/individual/", thecard, deepcopy(indfreeze), thegoodvalsdict, POIs, gr, redo, pretend) )



    # And finally, the globals
    glofreeze = []
    for gr in vl.global_list:
        glofreeze += gloGroup[gr]
        subtasks.append( (reluncpath + "/global/", thecard, deepcopy(glofreeze), thegoodvalsdict, POIs, gr, redo, pretend) )

    #print subtasks

    if nth > 1:
        pool = Pool(nth)
        pool.map(doMeThisFitPlease, subtasks)
        pool.close()
        pool.join()
    else:
        for tsk in subtasks:
            doMeThisFitPlease(tsk)


    ### First, we obtain the individual-separated files.
    #cumulative = []
    #for gr in vl.individual_list:
        #cumulative += indGroup[gr]
        #if not os.path.isfile(reluncpath + "/individual/higgsCombine{g}.FitDiagnostics.mH125.root".format(g = gr, p = poi)) or redo:
            #thecomm = basecommand.format(setpars     = ",".join([el + "=1" for el in POIs] + [el + "=" + str(thegoodvalsdict[el][0]) for el in cumulative if el in thegoodvalsdict]),
                                        #extra        = '-n {g} {card} --freezeParameters {c} --out {o}'.format(card = thecard,
                                                                                                               #g    = gr,
                                                                                                               #c    = ",".join(cumulative),
                                                                                                               #o    = reluncpath + "/individual/"),
                                        #minmaxlist   = ":".join([el + "=0,3" for el in POIs]),
            #print "\nCommand:", thecomm
            #if not pretend:
                #os.system(thecomm)
                #os.system("mv higgsCombine{u}.FitDiagnostics.mH125.root ".format(u = gr) + reluncpath + "/individual/")



    ## And finally, the globals
    #cumulative = []
    #for gr in vl.global_list:
        #cumulative += gloGroup[gr]
        #thecomm = basecommand.format(setpars      = ",".join([el + "=1" for el in POIs] + [el + "=" + str(thegoodvalsdict[el][0]) for el in cumulative if el in thegoodvalsdict]),
                                     #extra        = '-n {g} {card} --freezeParameters {c} --out {o}'.format(card = thecard,
                                                                                                            #g    = gr,
                                                                                                            #c    = ",".join(cumulative),
                                                                                                            #o    = reluncpath + "/global/"),
                                     ##extra        = '-n {g} {card} --freezeParameters {c} --snapshotName FitDiagnostics'.format(card = inpath + "/"+ iY + "/" + iV + "/sigextr_fit_combine/fitdiagnostics/higgsCombine" + iY + "_" + iV + ".FitDiagnostics.mH120.root", g = gr, c = ",".join(cumulative)),
                                     #minmaxlist   = ":".join([el + "=0,3" for el in POIs]),
        #)

        #print "\nCommand:", thecomm
        #if not pretend:
            #os.system(thecomm)
            #os.system("mv higgsCombine{u}.FitDiagnostics.mH125.root ".format(u = gr) + reluncpath + "/global/")
    return




if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage = "python nanoAOD_checker.py [options]", description = "Checker tool for the outputs of nanoAOD production (NOT postprocessing)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inpath',      '-i',  metavar = 'inpath',     dest = "inpath",      required = False, default = "./temp/differential")
    parser.add_argument('--pretend',     '-p',  action  = "store_true", dest = "pretend",     required = False, default = False)
    parser.add_argument('--nPoints',     '-P',  metavar = 'npoints',    dest = "npoints",     required = False, default = 128, type = int)
    parser.add_argument('--queue',       '-q',  action  = "store_true", dest = "queue",       required = False, default = False)
    parser.add_argument('--ncores',      '-j',  metavar = 'ncores',     dest = "ncores",      required = False, default = -1, type = int)
    parser.add_argument('--redoAll',     '-r',  action  = "store_true", dest = "redo",        required = False, default = False)
    parser.add_argument('--year',        '-y',  metavar = 'year',       dest = "year",        required = False, default = "all")
    parser.add_argument('--variable',    '-v',  metavar = 'variable',   dest = "variable",    required = False, default = "all")

    args     = parser.parse_args()
    pretend  = args.pretend
    inpath   = args.inpath
    npoints  = args.npoints
    doBatch  = args.queue
    nthreads = args.ncores
    variable = args.variable
    year     = args.year
    redo     = args.redo

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

        theyears = []
        presentyears = next(os.walk(inpath))[1]
        if "2022" in presentyears:
            theyears.append("2022")
        if "2022PostEE" in presentyears:
            theyears.append("2022PostEE")
        if "run3" in presentyears:
            theyears.append("run3")

        for iV in thevars:
            if any( [el in iV for el in vl.vetolist] ): continue
            tasks.append( (inpath, iY, iV, npoints, nthreads, redo) )


    """
    if nthreads > 1:
        pool = Pool(nthreads)
        pool.map(calculateRelativeUncertainties, tasks)
        pool.close()
        pool.join()
    else:
        for tsk in tasks:
            calculateRelativeUncertainties(tsk)
    """
    print(tasks)

    for tsk in tasks:
        calculateRelativeUncertainties(tsk)
