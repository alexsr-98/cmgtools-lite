import os, sys, enum, argparse
from multiprocessing import Pool
import warnings as wr
import ROOT as r

sys.path.append('{cmsswpath}/src/CMGTools/TTHAnalysis/python/plotter/twttbar-run2UL/differential/'.format(cmsswpath = os.environ['CMSSW_BASE']))
import varList as vl

r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(True)

friendspath   = "/lustrefs/hdd_pool_dir/nanoAODv9/bb4l/productions"
logpath       = friendspath + "/{p}/{y}/logs/cards_differential"

friendsscaff = "--FMCs {P}/0_jecs --Fs {P}/1_lepsuncsAndParticle --Fs {P}/2_cleaning --Fs {P}/3_varstrigger --FMCs {P}/4_scalefactors"

slurmscaff    = "sbatch -c {nth} -p {queue} -J {jobname} -e {logpath}/log.%j.%x.err -o {logpath}/log.%j.%x.out --wrap '{command}'"

commandscaff  = '''python3 makeShapeCards_TopRun2.py --tree NanoAOD {mcafile} {cutsfile} "{variable}" "{bins}" {samplespaths} {friends} --od {outpath} -l {lumi} {nth} -f -L twttbar-run2UL/functions_twttbar.cc --neg {weights} --year {year} {asimovornot} {uncs} {extra} {name} --AP --storeAll --notVarsChanges --threshold 0.001'''
#commandscaff  = '''python makeShapeCardsNew.py --tree NanoAOD {mcafile} {cutsfile} "{variable}" "{bins}" {samplespaths} {friends} --od {outpath} -l {lumi} {nth} -f -L twttbar-run2UL/functions_tw.cc --neg --threshold 0.01 {weights} --year {year} {asimovornot} {uncs} {extra} {name} --AP --storeAll'''

#### IMPORTANT NOTE: it is crucial to keep the xsec and lumi weights, although afterwards we
# will have to divide by the second one. The first one is relevant becasue of the dividing of the files.
# The second one, because of the proportions between MC simulations.

nomweight     = '''-W "MuonIDSF * MuonISOSF * ElecIDSF * ElecRECOSF * TrigSF * puWeight * bTagWeight * L1PreFiringWeight_Nom * jetPUidWeight"'''
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
    prod, year, variable, asimov, nthreads, outpath, region, noUnc, useFibre, extra, pretend, queue, njnt = tsk
    if queue == "":
        thecomm = CardsCommand(prod, year, variable, asimov, nthreads, outpath, region, noUnc, useFibre, extra, njnt)
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
                                    command = CardsCommand(prod, year, variable, asimov, nthreads, outpath, region, noUnc, useFibre, extra, njnt))

        print("Command: " + thecomm)

        if not pretend:
            os.system(thecomm)
    return


def CardsCommand(prod, year, var, isAsimov, nthreads, outpath, region, noUnc, useFibre, extra, njnt):
    mcafile_   = "twttbar-run2UL/differential/mca/mca-twttbar-diff.txt" if "forExtr" not in region and "control" not in region and region != "nonfiducialAll" else "twttbar-run2UL/differential/mca/mca-twttbar-diff_bs.txt"
    cutsfile_  = "twttbar-run2UL/differential/cuts/cuts-{reg}-{njnt}.txt".format(reg = region.replace("Response", "").replace("Square", "").replace("All", "") if ("forExtr" not in region and "but" not in region and "control" not in region) else "detectorparticle" if "but" in region else "detector",
                                                                                         njnt = njnt)

    samplespaths_ = "-P " + friendspath + "/" + prod + ("/" + year) * (year != "run2")
    if useFibre: samplespaths_ = samplespaths_.replace("phedexrw", "phedex").replace("cienciasrw", "ciencias")

    nth_       = "" if nthreads == 0 else ("--split-factor=-1 -j " + str(nthreads))
    friends_   = friendsscaff
    outpath_   = outpath + "/" + year + "/" + (var if "control" not in region else ("controlReg" + vl.diffControlReg))

    if not os.path.isdir(outpath_):
        os.system("mkdir -p " + outpath_)

    thebins = (vl.varList[var]["bins_detector"] if (region == "detector" or region == "detectorparticlebutdetector" or "nonfiducial" in region or "forExtr" in region or "control" in region) else
               vl.varList[var]["bins_particle"] if (region == "particle" or region == "detectorparticle") else
              (vl.varList[var]["bins_particle"], vl.varList[var]["bins_particle"]) if region == "detectorparticleResponseSquare" else
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
                  nomweight if region == "detectorparticleResponse" or "nonfiducial" in region or region == "detectorparticleResponseSquare" else
                  #nomweight if region == "detectorparticleResponse" else
                  genweight)

    extra_ = extra
    if region not in ["detector", "particle"]:
        extra_ += " --xp twttbardr,twttbards,twttbarherwig,twttbaramc_dr,twttbaramc_dr2,twttbaramc_ds,twttbaramc_ds_is,twttbaramc_ds_runningBW,twttbaramc_ds_is_runningBW,ttbar,ttbarherwig,ttbaramc,twdr,twds,twherwig,twamc_dr,twamc_ds,twamc_dr2,twamc_ds_is,twamc_ds_is_runningBW,twamc_ds_runningBW,bb4lv1"

    comm = commandscaff.format(outpath      = outpath_,
                               friends      = friends_,
                               samplespaths = samplespaths_,
                               lumi      = (vl.LumiDict[year] if year != "run2" else str(vl.LumiDict["2016apv"]) + "," + str(vl.LumiDict["2016"]) + "," + str(vl.LumiDict["2017"]) + "," + str(vl.LumiDict["2018"])) if "forExtr" not in region else (vl.LumiDict[year] if year != "run2" else str(vl.LumiDict["2016apv"]) + "," + str(vl.LumiDict["2016"]) + "," + str(vl.LumiDict["2017"]) + "," + str(vl.LumiDict["2018"])),
                               variable  = variable_,
                               bins      = bins_,
                               nth       = nth_,
                               year      = year if year != "run2" else "2016apv,2016,2017,2018",
                               asimovornot = "--asimov s+b" if isAsimov or region == "nonfiducialAll" else "--asimov s+b" if region != "forExtr" else "",
                               mcafile   = mcafile_,
                               cutsfile  = cutsfile_,
                               uncs      = "" if (region == "particle" or region == "detectorparticleResponseSquare") else "--unc twttbar-run2UL/differential/uncs/uncs-twttbar.txt --amc" if not noUnc and "Response" in region else "--unc twttbar-run2UL/differential/uncs/uncs-twttbar.txt --amc" if not noUnc else "--amc",
                               name      = name_,
                               weights   = weights_,
                               extra     = extra_)

    return comm


if __name__=="__main__":
    parser = argparse.ArgumentParser(usage = "python nanoAOD_checker.py [options]", description = "Checker tool for the outputs of nanoAOD production (NOT postprocessing)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--production','-P',  metavar = "prod",       dest = "prod",     required = True)
    parser.add_argument('--year',      '-y',  metavar = 'year',       dest = "year",     required = False, default = "2016")
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
    parser.add_argument('--njetsnbjets',      '-njnt',  metavar = 'njnt',       dest = "njnt",     required = False, default = "2j2t")



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
    njnt      = args.njnt

    theregs  = ["detector", "particle", "detectorparticleResponse", "detectorparticlebutdetector",
                "detectorparticle", "nonfiducial", "forExtr", "detectorparticleResponseSquare"]

    thevars  = vl.varList["Names"]["Variables"]
    theyears = ["2016", "2016apv", "2017", "2018", "run2"]
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
                tasks.append( (prod, yr, var, asimov, nthreads, outpath, reg, noUnc, useFibre, extra, pretend, queue, njnt) )
            else:
                for var in thevars:
                    tasks.append( (prod, yr, var, asimov, nthreads, outpath, reg, noUnc, useFibre, extra, pretend, queue, njnt) )

    #print tasks
    calculate = True
    for task in tasks:
        print("\nProcessing " + str(task) + "\n")

        #if str(task) == "('2020-09-20', 'run2', 'Lep1Lep2Jet1MET_Mt', True, 87, 'temp_2020_10_29_pruebasdiff', 'forExtr', False, True, '', False, '')":
            #calculate = True

        if calculate:
            ExecuteOrSubmitTask(task)
            #sys.exit()
