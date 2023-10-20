# -*- coding: UTF-8 -*-.
#!/usr/bin/env python
import os, sys, enum, argparse
from multiprocessing import Pool
import warnings as wr
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(True)

#### Settings
friendspath  = "/lustrefs/hdd_pool_dir/nanoAODv11/tw-run3/productions"

logpath      = friendspath + "/{p}/{y}/logs/cards_inclusive"

lumidict     = {"2022"       : 7.78,
                "2022PostEE" : 20.67}

friendsscaff = "--Fs {P}/0_jecs --Fs {P}/1_lepsuncsAndParticle --Fs {P}/2_cleaning --Fs {P}/3_varstrigger --FMCs {P}/4_scalefactors --Fs {P}/6_mvasRF" # --Fs {P}/5_mvas

commandscaff = '''python3 makeShapeCards_TopRun2.py --tree NanoAOD {mcafile} {cutsfile} "{variable}" "{bins}" {samplespaths} {friends} --od {outpath} -l {lumi} {nth} -f -L tw-run3/functions_tw.cc --neg --threshold 0.01 -W "MuonIDSF * MuonISOSF * ElecIDSF * ElecRECOSF * TrigSF * bTagWeight * puWeight" --year {year} {asimovornot} {pseudodataorNot} {uncs} {extra} --AP --storeAll'''

slurmscaff   = "sbatch -c {nth} -p {queue} -J {jobname} -e {logpath}/log.%j.%x.err -o {logpath}/log.%j.%x.out --wrap '{command}'"

listofforcedshape = "btagging_2022PostEE,btagging_2022,btagging_corr,elecidsf,elecrecosf,fsr,isr_ttbar,isr_tw,jer,jes_Absolute,jes_Absolute_2022PostEE,jes_Absolute_2022,jes_BBEC1,jes_BBEC1_2022PostEE,jes_EC2,jes_EC2_2022PostEE,jes_EC2_2022,jes_FlavorQCD,jes_HF,jes_HF_2022PostEE,jes_HF_2022,jes_RelativeBal,jes_RelativeSample_2022PostEE,jes_RelativeSample_2022,lumi,mistagging_2022PostEE,mistagging_2022,mistagging_corr,mtop,muonidsf,pdfhessian,pileup,topptrew,triggereff,ttbar_scales,tw_scales,ds,colour_rec_erdon,colour_rec_cr1,colour_rec_cr2"




def GeneralExecutioner(task):
    prod, year, variable, bines, asimov, pseudoData, nthreads, outpath, region, queue, extra, pretend, useFibre, noUnc = task

    if not os.path.isdir(outpath + "/" + str(year)) and not pretend:
        os.system("mkdir -p " + outpath + "/" + str(year))

    if queue != "":
        if not os.path.isdir(logpath.format(y = year, p = prod)) and not pretend:
            os.system("mkdir -p " + logpath.format(y = year, p = prod))

        jobname_   = "Card_{r}_{y}".format(y = year, r = region)
        submitcomm = slurmscaff.format(nth  = nthreads,
                                    queue   = queue,
                                    jobname = jobname_,
                                    logpath = logpath.format(y = year, p = prod),
                                    command = CardsCommand(prod, year, variable, bines, asimov, pseudoData, nthreads, outpath, region, noUnc, useFibre, extra))
        print("Command:", submitcomm)
        if not pretend: os.system(submitcomm)
    else:
        execcomm = CardsCommand(prod, year, variable, bines, asimov, pseudoData, nthreads, outpath, region, noUnc, useFibre, extra)
        print("Command:", execcomm)
        if not pretend: os.system(execcomm)


    return






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




def CardsCommand(prod, year, var, bines, isAsimov, isPseudoData, nthreads, outpath, region, noUnc, useFibre, extra):
    mcafile_   = "tw-run3/mca-tw.txt"
    cutsfile_  = "tw-run3/cuts-tw-{reg}.txt".format(reg = region)

    samplespaths_ = "-P " + friendspath + "/" + prod + ("/" + year) * (year != "run3")
    if useFibre: samplespaths_ = samplespaths_.replace("phedexrw", "phedex").replace("cienciasrw", "ciencias")

#    nth_       = "" if nthreads == 0 else ("--split-factor=-1 -j " + str(nthreads))
    nth_       = "" if nthreads == 0 else ("--split-factor=0 -j " + str(nthreads))
    friends_   = friendsscaff
    outpath_   = outpath + "/" + year + "/" + region

    extra_ = extra + " --forceShape {l}".format(l = listofforcedshape)

    comm = commandscaff.format(outpath      = outpath_,
                               friends      = friends_,
                               samplespaths = samplespaths_,
                               lumi      = lumidict[year] if year != "run3" else str(lumidict["2022"]) + "," + str(lumidict["2022PostEE"]),
                               variable  = var,
                               bins      = bines,
                               nth       = nth_,
                               year      = year if year != "run3" else "2022,2022PostEE",
                               asimovornot = "--asimov s+b" if isAsimov else "",
                               pseudodataorNot = "--pseudoData s+b" if isPseudoData else "",
                               mcafile   = mcafile_,
                               cutsfile  = cutsfile_,
                               #uncs      = "--unc tw-run2/uncs-tw_{r}mvaNEW.txt --amc".format(r = region) if not noUnc else "--amc",
                               uncs      = "--unc tw-run3/uncs-tw.txt --amc" if not noUnc else "--amc",
                               extra     = extra_)

    return comm


def ExecuteOrSubmitTask(tsk):
    prod, year, variable, bines, asimov, pseudoData, nthreads, outpath, region, noUnc, useFibre, extra, pretend, queue = tsk

    if not os.path.isdir(outpath + "/" + year + "/" + region):
        os.system("mkdir -p " + outpath + "/" + year + "/" + region)

    if queue == "":
        thecomm = CardsCommand(prod, year, variable, bines, asimov, pseudoData, nthreads, outpath, region, noUnc, useFibre, extra)
        print("Command: " + thecomm)

        if not pretend:
            os.system(thecomm)

    else:
        if not os.path.isdir(logpath.format(p = prod, y = yr)):
            os.system("mkdir -p " + logpath.format(p = prod, y = yr))

        thecomm = slurmscaff.format(nth     = nthreads,
                                    queue   = queue,
                                    jobname = "CMGTcards_" + year + "_" + (variable.replace("(", "").replace(")", "") if not "min" in variable else "Jet2_Pt") + "_" + region,
                                    logpath = logpath.format(p = prod, y = yr),
                                    command = CardsCommand(prod, year, variable, bines, asimov, pseudoData, nthreads, outpath, region, noUnc, useFibre, extra))

        print("Command: " + thecomm)

        if not pretend:
            os.system(thecomm)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage = "python3 nanoAOD_checker.py [options]", description = "Checker tool for the outputs of nanoAOD production (NOT postprocessing)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--production','-P', metavar = "prod",       dest = "prod",     required = True)
    parser.add_argument('--year',      '-y', metavar = 'year',       dest = "year",     required = False, default = "2022")
    parser.add_argument('--queue',     '-q', metavar = 'queue',      dest = "queue",    required = False, default = "")
    parser.add_argument('--extraArgs', '-e', metavar = 'extra',      dest = "extra",    required = False, default = "")
    parser.add_argument('--nthreads',  '-j', metavar = 'nthreads',   dest = "nthreads", required = False, default = 0, type = int)
    parser.add_argument('--pretend',   '-p', action  = "store_true", dest = "pretend",  required = False, default = False)
    parser.add_argument('--outpath',   '-o', metavar = 'outpath',    dest = "outpath",  required = False, default = "./temp/cards")
    parser.add_argument('--region',    '-r', metavar = 'region',     dest = "region",   required = False, default = "1j1t")
    parser.add_argument('--nounc',     '-nu',action  = "store_true", dest = "nounc",    required = False, default = False)
    parser.add_argument('--variable',  '-v', metavar = 'variable',   dest = "variable", required = False, default = "getBDtW(tmvaBDT_1j1b)")
    parser.add_argument('--bines',     '-b', metavar = 'bines',      dest = "bines",    required = False, default = "")
    parser.add_argument('--asimov',    '-a', action  = "store_true", dest = "asimov",   required = False, default = False)
    parser.add_argument('--pseudoData','-pd', action = "store_true", dest = "pseudoData", required = False, default = False)
    parser.add_argument('--useFibre',  '-f', action  = "store_true", dest = "useFibre", required = False, default = False)


    args     = parser.parse_args()
    prod     = args.prod
    year     = args.year
    queue    = args.queue
    extra    = args.extra
    nthreads = args.nthreads
    pretend  = args.pretend
    outpath  = args.outpath
    region   = args.region
    noUnc    = args.nounc
    variable = args.variable
    bines    = args.bines
    asimov   = args.asimov
    pseudoData = args.pseudoData
    useFibre = args.useFibre

    # Check that asimov and pseudoData are not both true
    if asimov and pseudoData:
        print("ERROR: asimov and pseudoData cannot be both true.")
        sys.exit()

    #print CardsCommand(prod, year, variable, bines, asimov, nthreads, outpath, region, noUnc, useFibre, extra)

    theregs  = ["1j1t", "2j1t", "2j2t"]
    #thevars  = ["getBDtW(tmvaBDT_1j1b)", "getBDtWOther(tmvaBDT_2j1b)", "min(max(Jet2_Pt, 30.), 189.)"] #Actual
    thevars  = ["getRFtW(mvaRF_1j1b)", "getRFtWOther(mvaRF_2j1b)", "min(max(Jet2_Pt, 30.), 189.)"] #With RF
    thebins  = ["[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5]",
                "[0.5,1.5,2.5,3.5,4.5,5.5,6.5]",
                "[30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.]"]              
    theyears = ["2022", "2022PostEE"]
    tasks    = []
    
    ##### Attempt without MVAs
    #thevars  = ["min(max(LepGood_pt_corrAll[0], 20.), 139.)", "min(max(LepGood_pt_corrAll[0], 20.), 139.)", "min(max(Jet2_Pt, 30.), 189.)"]
    #thebins  = ["[20,30,40,50,60,70,80,90,100,110,120,130,140]",
    #            "[20,30,40,50,60,70,80,90,100,110,120,130,140]",
    #            "[30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.]"]         

    ##### Attempt with SF channels also
    theregs  = ["1j1t", "2j1t", "2j2t", "1j1t-mm", "1j1t-ee"]
    thevars  = ["getRFtW(mvaRF_1j1b)", "getRFtWOther(mvaRF_2j1b)", "min(max(Jet2_Pt, 30.), 189.)", "getRFtW_mm(mvaRF_1j1b_mm)", "getRFtW_ee(mvaRF_1j1b_ee)"] #With RF
    thebins  = ["[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5]",
                "[0.5,1.5,2.5,3.5,4.5,5.5,6.5]",
                "[30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.]",
                "[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5]",
                "[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5]"]              
    theyears = ["2022", "2022PostEE"]
    tasks    = []

    ##### Attempt with extra emu regions
    theregs  = ["1j1t", "2j1t", "2j2t", "1j0t","3j1t","3j2t"]
    thevars  = ["getRFtW(mvaRF_1j1b)", "getRFtWOther(mvaRF_2j1b)", "min(max(Jet2_Pt, 30.), 189.)", "getRFtW1j0b(mvaRF_1j1b)", "min(max(Jet_pt_nom[iJetSel30_Recl[1]], 30.), 149.)", "min(max(Jet_pt_nom[iJetSel30_Recl[2]], 30.), 149.)"] #With RF
    thebins  = ["[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5]",
                "[0.5,1.5,2.5,3.5,4.5,5.5,6.5]",
                "[30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.]",
                "[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5]",
                "[30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.]",
                "[30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.]"]              
    theyears = ["2022", "2022PostEE"]
    tasks    = []

    ##### Attempt with extra emu regions
    theregs  = ["1j1t", "2j1t", "2j2t","3j2t"]
    thevars  = ["getRF1j1b_opt_20bins(mvaRF_1j1b)", "getRF2j1b_opt_12bins(mvaRF_2j1b)", "min(max(minimax, 0.), 449.)", "min(max(Jet_pt_nom[iJetSel30_Recl[1]], 30.), 149.)"] #With RF
    thebins  = ["[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5]",
                "[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5]",
                "[0.0, 32.142857142857146, 64.28571428571429, 96.42857142857144, 128.57142857142858, 160.71428571428572, 192.8571428571429, 225.00000000000003, 257.14285714285717, 289.28571428571433, 321.42857142857144, 353.5714285714286, 385.7142857142858, 417.8571428571429, 450.0]",
                "[30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.]"]              
    theyears = ["2022", "2022PostEE"]
    tasks    = []


    ##### Attempt with SF channels also  (only SF)
    ##theregs  = ["1j1t-mm", "1j1t-ee"]
    ##thevars  = ["getRFtW_mm(mvaRF_1j1b_mm)", "getRFtW_ee(mvaRF_1j1b_ee)"] #With RF
    ##thebins  = ["[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5]",
    ##            "[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5]"]              
    ##theyears = ["2022", "2022PostEE"]
    ##tasks    = []
    
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

    if bines != "":
        thebins = [bines]
    
    for yr in theyears:
        for i in range(len(theregs)):
            if "{year}" in thevars[i] and yr != "run3":
                tasks.append( (prod, yr, thevars[i].format(year = yr), thebins[i], asimov, pseudoData, nthreads, outpath, theregs[i], noUnc, useFibre, extra, pretend, queue) )
            elif "{year}" in thevars[i] and yr == "run3":
                tasks.append( (prod, yr, thevars[i].format(year = ""), thebins[i], asimov, pseudoData, nthreads, outpath, theregs[i], noUnc, useFibre, extra, pretend, queue) )
            else:
                tasks.append( (prod, yr, thevars[i], thebins[i], asimov, pseudoData, nthreads, outpath, theregs[i], noUnc, useFibre, extra, pretend, queue) )

    #print tasks
    calculate = True
    for task in tasks:
        print("\nProcessing " + str(task) + "\n")

        if calculate:
            ExecuteOrSubmitTask(task)
