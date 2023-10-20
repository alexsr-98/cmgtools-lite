import os, sys, enum, argparse
from multiprocessing import Pool
import warnings as wr
import ROOT as r

sys.path.append('{cmsswpath}/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/differential/'.format(cmsswpath = os.environ['CMSSW_BASE']))
import varList as vl

r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(True)

friendspath   = "/lustrefs/hdd_pool_dir/nanoAODv11/tw-run3/productions"
logpath       = friendspath + "/{p}/{y}/logs/cards_differential"

friendsscaff = "--Fs {P}/0_jecs --Fs {P}/1_lepsuncsAndParticle --Fs {P}/2_cleaning --Fs {P}/3_varstrigger --FMCs {P}/4_scalefactors"

slurmscaff    = "sbatch -c {nth} -p {queue} -J {jobname} -e {logpath}/log.%j.%x.err -o {logpath}/log.%j.%x.out --wrap '{command}'"

commandscaff  = '''python3 makeShapeCards_TopRun2.py --tree NanoAOD {mcafile} {cutsfile} "{variable}" "{bins}" {samplespaths} {friends} --od {outpath} -l {lumi} {nth} -f -L tw-run3/functions_tw.cc --neg {weights} --year {year} {asimovornot} {uncs} {extra} {name} --AP --storeAll --notVarsChanges --threshold 0.001'''
#### IMPORTANT NOTE: it is crucial to keep the xsec and lumi weights, although afterwards we
# will have to divide by the second one. The first one is relevant becasue of the dividing of the files.
# The second one, because of the proportions between MC simulations.

nomweight     = '''-W "MuonIDSF * MuonISOSF * ElecIDSF * ElecRECOSF * TrigSF * bTagWeight * puWeight"'''
genweight     = ""



def PythonListToString(theL):
    ret = "["
    for el in theL:
        ret += str(el)
        if el != theL[-1]:
            ret += ", "
    ret += "]"
    return ret


def ExecuteOrSubmitTask(tsk):
    prod, year, variable, asimov, nthreads, outpath, region, noUnc, useFibre, extra, pretend, queue = tsk
    if queue == "":
        thecomm = CardsCommand(prod, year, variable, asimov, nthreads, outpath, region, noUnc, useFibre, extra)
        print("Command: " + thecomm)

        if not pretend:
            os.system(thecomm)

    else:
        if not os.path.isdir(logpath.format(p = prod, y = yr)):
            os.system("mkdir -p " + logpath.format(p = prod, y = yr))

        thecomm = slurmscaff.format(nth     = nthreads,
                                    queue   = queue,
                                    jobname = "CMGTcardsforunfolding_" + variable + "_" + region,
                                    logpath = logpath.format(p = prod, y = yr),
                                    command = CardsCommand(prod, year, variable, asimov, nthreads, outpath, region, noUnc, useFibre, extra))

        print("Command: " + thecomm)

        if not pretend:
            os.system(thecomm)
    return


def CardsCommand(prod, year, var, isAsimov, nthreads, outpath, region, noUnc, useFibre, extra):
    mcafile_   = "tw-run3/differential/mca-differential/mca-tw-diff.txt" if "forExtr" not in region and "control" not in region and region != "nonfiducialAll" else "tw-run3/mca-tw.txt"
    cutsfile_  = "tw-run3/differential/cuts-differential/cuts-{reg}-{njnt}.txt".format(reg = region.replace("Response", "").replace("All", "") if ("forExtr" not in region and "but" not in region and "control" not in region) else "detectorparticle" if "but" in region else "detector",
                                                                                                  njnt = "1j1t" if "control" not in region else vl.diffControlReg )

    samplespaths_ = "-P " + friendspath + "/" + prod + ("/" + year) * (year != "run3")
    if useFibre: samplespaths_ = samplespaths_.replace("phedexrw", "phedex").replace("cienciasrw", "ciencias")

    nth_       = "" if nthreads == 0 else ("--split-factor=-1 -j " + str(nthreads))
    friends_   = friendsscaff
    outpath_   = outpath + "/" + year + "/" + (var if "control" not in region else ("controlReg" + vl.diffControlReg))

    if not os.path.isdir(outpath_):
        os.system("mkdir -p " + outpath_)

    thebins = (vl.varList[var]["bins_detector"] if (region == "detector" or region == "detectorparticlebutdetector" or "nonfiducial" in region or "forExtr" in region or "control" in region) else
               vl.varList[var]["bins_particle"] if (region == "particle" or region == "detectorparticle") else
              (vl.varList[var]["bins_particle"], vl.varList[var]["bins_detector"]) )

    bins_      = ""
    if type(thebins) != tuple:
        bins_ = PythonListToString(thebins)
    else:
        bins_ = PythonListToString(thebins[0]) + "*" + PythonListToString(thebins[1])

    variable_  = (vl.varList[var]["var_detector"] if (region == "detector" or region == "detectorparticlebutdetector" or "nonfiducial" in region or "forExtr" in region or "control" in region) else
                  vl.varList[var]["var_particle"] if (region == "particle" or region == "detectorparticle") else
                  #vl.varList[var]["var_particle"] + ":" + vl.varList[var]["var_detector"])
                  vl.varList[var]["var_detector"] + ":" + vl.varList[var]["var_particle"])
    name_      = "--binname " + region
    weights_   = (nomweight if (region == "detector" or "forExtr" in region or region == "detectorparticlebutdetector" or "control" in region) else
                  nomweight if region == "detectorparticleResponse" or "nonfiducial" in region else
                  #nomweight if region == "detectorparticleResponse" else
                  genweight)

    extra_ = extra
    if region not in ["detector", "particle"]:
        extra_ += " --xp twds,twherwig,twamcatnlo_dr,twamcatnlo_dr2,twamcatnlo_ds,twamcatnlo_ds_is,twamcatnlo_ds_runningBW,twamcatnlo_ds_is_runningBW"

    comm = commandscaff.format(outpath      = outpath_,
                               friends      = friends_,
                               samplespaths = samplespaths_,
                               lumi      = (vl.LumiDict[year] if year != "run3" else str(vl.LumiDict["2022"]) + "," + str(vl.LumiDict["2022PostEE"])) if "forExtr" not in region else (vl.LumiDict[year] if year != "run3" else str(vl.LumiDict["2022"]) + "," + str(vl.LumiDict["2022PostEE"])),
                               variable  = variable_,
                               bins      = bins_,
                               nth       = nth_,
                               year      = year if year != "run3" else "2022,2022PostEE",
                               asimovornot = "--asimov s+b" if isAsimov or region == "nonfiducialAll" else "--asimov s+b" if region != "forExtr" else "",
                               mcafile   = mcafile_,
                               cutsfile  = cutsfile_,
                               #uncs      = "" if region == "particle" else "--unc tw-run3/uncs-tw.txt --amc" if not noUnc and "Response" in region else "--unc tw-run3/uncs-tw_1j1tdiff_" + var.lower().replace("_", "") + ".txt --amc" if not noUnc else "--amc",
                               uncs      = "" if region == "particle" else "--unc tw-run3/uncs-tw.txt --amc" if not noUnc and "Response" in region else "--unc tw-run3/uncs-tw.txt --amc" if not noUnc else "--amc",
                               name      = name_,
                               weights   = weights_,
                               extra     = extra_)

    return comm


if __name__=="__main__":
    parser = argparse.ArgumentParser(usage = "python3 nanoAOD_checker.py [options]", description = "Checker tool for the outputs of nanoAOD production (NOT postprocessing)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--production','-P',  metavar = "prod",       dest = "prod",     required = True)
    parser.add_argument('--year',      '-y',  metavar = 'year',       dest = "year",     required = False, default = "2022")
    parser.add_argument('--queue',     '-q',  metavar = 'queue',      dest = "queue",    required = False, default = "")
    parser.add_argument('--extraArgs', '-e',  metavar = 'extra',      dest = "extra",    required = False, default = "")
    parser.add_argument('--nthreads',  '-j',  metavar = 'nthreads',   dest = "nthreads", required = False, default = 0, type = int)
    parser.add_argument('--pretend',   '-p',  action  = "store_true", dest = "pretend",  required = False, default = False)
    parser.add_argument('--outpath',   '-o',  metavar = 'outpath',    dest = "outpath",  required = False, default = "./temp/differential/cards")
    parser.add_argument('--region',    '-r',  metavar = 'region',     dest = "region",   required = False, default = "all")
    parser.add_argument('--nounc',     '-nu', action  = "store_true", dest = "nounc",    required = False, default = False)
    parser.add_argument('--variable',  '-v',  metavar = 'variable',   dest = "variable", required = False, default = "Lep1Lep2_DPhi")
    parser.add_argument('--asimov',    '-a',  action  = "store_true", dest = "asimov",   required = False, default = False)
    parser.add_argument('--useFibre',  '-f',  action  = "store_true", dest = "useFibre", required = False, default = False)
    parser.add_argument('--pureReg',   '-pr', action = "store_true",  dest = "doPureReg",required = False, default = False)


    args      = parser.parse_args()
    prod      = args.prod
    year      = args.year
    queue     = args.queue
    extra     = args.extra
    nthreads  = args.nthreads
    pretend   = args.pretend
    outpath   = args.outpath
    region    = args.region
    noUnc     = args.nounc
    variable  = args.variable
    asimov    = args.asimov
    useFibre  = args.useFibre
    doPureReg = args.doPureReg

    theregs  = ["detector", "particle", "detectorparticleResponse", "detectorparticlebutdetector",
    ##theregs  = ["detector", "detectorparticleResponse", "detectorparticlebutdetector",
                #"detectorparticle", "nonfiducial"]#, "forExtr", "controlReg"]
                "detectorparticle", "nonfiducial", "forExtr"]
                ##"detectorparticle", "nonfiducial", "forExtr", "controlReg"]


    #theregs  = ["detectorparticleResponse"]
    #theregs  = ["detector", "particle", "detectorparticlebutdetector",
                #"detectorparticle", "nonfiducial"]

    thevars  = vl.varList["Names"]["Variables"]
    theyears = ["2022", "2022PostEE", "run3"]
    tasks    = []

    if variable.lower() != "all":
        if "," in variable:
            thevars = variable.split(",")
        else:
            thevars = [ variable ]

    if region.lower() != "all":
        if "," in region:
            theregs = region.split(",")
        else:
            theregs = [ region ]

    if year.lower() != "all":
        if "," in year:
            theyears = year.split(",")
        else:
            theyears = [ year ]


    for reg in theregs:
        for yr in theyears:
            if reg == "controlReg":
                var = vl.ControlRegionVariableDict[vl.diffControlReg]
                tasks.append( (prod, yr, var, asimov, nthreads, outpath, reg, noUnc, useFibre, extra, pretend, queue) )
            else:
                for var in thevars:
                    tasks.append( (prod, yr, var, asimov, nthreads, outpath, reg, noUnc, useFibre, extra, pretend, queue) )

    #print tasks
    calculate = True
    for task in tasks:
        print("\nProcessing " + str(task) + "\n")

        #if str(task) == "('2020-09-20', 'run2', 'Lep1Lep2Jet1MET_Mt', True, 87, 'temp_2020_10_29_pruebasdiff', 'forExtr', False, True, '', False, '')":
            #calculate = True

        if calculate:
            ExecuteOrSubmitTask(task)
            #sys.exit()
