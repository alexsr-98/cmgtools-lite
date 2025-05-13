# -*- coding: UTF-8 -*-.
#!/usr/bin/env python
import os, sys, enum, argparse
from multiprocessing import Pool
import warnings as wr
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(True)

#### Settings
friendspath  = "/lustrefs/hdd_pool_dir/nanoAODv12/tw-2022Run3/productions"
#friendspath = "/pool/phedexrw/userstorage/asoto/tw-run3/productions"

logpath      = friendspath + "/{p}/{y}/logs/plots"

lumidict     = {"2022"       : 8.0,
                "2022PostEE" : 26.67,}
                #"2022PostEE" : 20.67}
lumidict_byEra = {"2022C": 5.0707,
                  "2022D": 3.0063,
                  "2022E": 5.8783,
                  "2022F": 18.0070,
                  "2022G": 3.1219,}

friendsscaff = "--Fs {P}/0_jecs --Fs {P}/1_lepsuncsAndParticle --Fs {P}/2_cleaning --Fs {P}/3_varstrigger --FMCs {P}/4_scalefactors_full2022_comb --FMCs {P}/8_topPtReweight13to13p6" #--FMCs {P}/7_caloWeight ## CentralCaloWeight  /// puWeight //// CentralPVWeight

#-W 'MuonIDSF * MuonISOSF * ElecIDSF * ElecRECOSF * TrigSF * puWeight * bTagWeight * TopPtWeight'

commandscaff = "python3 mcPlots.py --tree NanoAOD --pdir {outpath} {friends} {samplespaths} -f -l {lumi} {nth} --year {year} {ratio} -W 'MuonIDSF * MuonISOSF * ElecIDSF * ElecRECOSF * TrigSF * puWeight * bTagWeight * TopPtWeightNNLO * TopPtWeight13to13p6' --showRatio --fixRatioRange --legendColumns 1 --legendWidth 0.15 --legendFontSize 0.05616 --noCms --topSpamSize 1.1 --addSpamSize 1.2 --lspam '#splitline{{#scale[1.2]{{#bf{{CMS}}}}}}{{#scale[1.0]{{#it{{{cmsText}}}}}}}' --neg --showMCError -L tw-run3/functions_tw.cc {selplot} {mcafile} {cutsfile} {plotsfile} {extra} --AP --noStatTotLegendOnRatio --addspam '{nameregion}' --lspamPosition 0.21 .835 .35 .875 --TotalUncRatioStyle 3444 0 --noStatUncOnRatio --YTitleOffset 1.8 1.8 --CanvasSize 600 450 --TotalUncRatioColor 922 922 --addspamPosition .365 .845 .55 .885 --transparentLegend --PrincipalPadDimensions 0.00 0.25 1.00 1.00 --RatioPadDimensions 0.00 0.00 1.00 0.25 --LeftRightMargins 0.16 0.03 --ratioYLabel 'Data / Pred.     ' --labelsSize 24 --axisTitleSize 26 --labelsFont 43 --BottomMarginRatio 0.42 --XTitleOffsetRatio 1.0 --noXErrData --printBin 'bin' --printBinUnity --noExpoShift --no-elist --maxDigitsYaxis 4 --UncLegendName Uncertainty --ratioLinePersonalization 2 1 1 --lspamShiftLegend --noAddspamShift --ratioYNDiv 503 --histoLineWidth 0 --removeMarkerUncertainty --centerRatioYTitle --perBin" # --NotDrawRatioLine --lspam '#scale[1.1]{{#bf{{CMS}}}}'

# For blind plots
commandblindscaff = "python3 mcPlots.py --tree NanoAOD --pdir {outpath} {friends} {samplespaths} -f -l {lumi} {nth} --year {year} {ratio} -W 'MuonIDSF * ElecIDSF * TrigSF' --showRatio --fixRatioRange --legendColumns 1 --legendWidth 0.07 --legendFontSize 0.032 --noCms --topSpamSize 1.1 --lspam '#splitline{{#scale[1.1]{{#bf{{CMS}}}}}}{{#scale[0.9]{{#it{{{cmsText}}}}}}}' --neg --showMCError -L tw-run3/functions_tw.cc {selplot} {mcafile} {cutsfile} {plotsfile} {extra} --AP --noStatTotLegendOnRatio --addspam '{nameregion}' --lspamPosition 0.21 .845 .35 .885 --TotalUncRatioStyle 3444 0 --noStatUncOnRatio --YTitleOffset 2.1 2.1 --TotalUncRatioColor 922 922 --addspamPosition .41 .855 .6 .895 --transparentLegend --LeftRightMargins 0.16 0.03  --labelsSize 22 --labelsFont 43 --BottomMarginRatio 0.42 --XTitleOffsetRatio 4.8 --noXErrData --printBin 'bin' --printBinUnity --noExpoShift --no-elist --maxDigitsYaxis 4 --UncLegendName Uncertainty --ratioLinePersonalization 2 1 1 --lspamShiftLegend --noAddspamShift --ratioYNDiv 503 --histoLineWidth 0 --removeMarkerUncertainty --centerRatioYTitle --perBin" #--ratioYNDiv 210 --NotDrawRatioLine

# For trigger efficiency plots
commandtriggscaff = "python3 mcPlots.py --tree NanoAOD --pdir {outpath} {friends} {samplespaths} -f -l {lumi} {nth} --year {year} {ratio} -W 'MuonIDSF * MuonISOSF * ElecIDSF * ElecRECOSF * puWeight' --showRatio --fixRatioRange --legendColumns 1 --legendWidth 0.15 --legendFontSize 0.05616 --noCms --topSpamSize 1.1 --addSpamSize 1.2 --lspam '#splitline{{#scale[1.1]{{#bf{{CMS}}}}}}{{#scale[0.9]{{#it{{{cmsText}}}}}}}' --neg -L tw-run3/functions_tw.cc {selplot} {mcafile} {cutsfile} {plotsfile} {extra} --AP --noStatTotLegendOnRatio --addspam '{nameregion}' --lspamPosition 0.21 .845 .35 .885 --TotalUncRatioStyle 3444 0 --noStatUncOnRatio --YTitleOffset 1.8 1.8 --CanvasSize 600 450 --TotalUncRatioColor 922 922 --addspamPosition .365 .845 .55 .885 --transparentLegend --PrincipalPadDimensions 0.00 0.25 1.00 1.00 --RatioPadDimensions 0.00 0.00 1.00 0.25 --LeftRightMargins 0.16 0.03 --ratioYLabel 'Data / Pred.     ' --labelsSize 22 --axisTitleSize 26 --labelsFont 43 --BottomMarginRatio 0.42 --XTitleOffsetRatio 1.0 --noXErrData --printBin 'bin' --printBinUnity --noExpoShift --no-elist --maxDigitsYaxis 4 --UncLegendName Uncertainty --ratioLinePersonalization 2 1 1 --lspamShiftLegend --noAddspamShift --ratioYNDiv 503 --histoLineWidth 0 --removeMarkerUncertainty --centerRatioYTitle --perBin" # --NotDrawRatioLine --lspam '#scale[1.1]{{#bf{{CMS}}}}'


slurmscaff   = 'sbatch {extraS} -c {nth} -p {queue} -J {jobname} -e {logpath}/log.%j.%x.err -o {logpath}/log.%j.%x.out --wrap "{command}"'



def GeneralExecutioner(task):
    prod, year, nthreads, outpath, selplot, region, ratiorange, queue, extra, pretend, useFibre, doUncs, doBlind, extraslurm, era, CMSText = task

    if era == None: # if era is not specified, use the year
        era = year
    
    if not os.path.isdir(outpath + "/" + str(era)) and not pretend:
        os.system("mkdir -p " + outpath + "/" + str(era))

    if queue != "":
        if not os.path.isdir(logpath.format(y = year, p = prod)) and not pretend:
            os.system("mkdir -p " + logpath.format(y = year, p = prod))

        for reg in region.split(","):
            jobname_   = "CMGTplotter_{y}_{r}_{s}".format(y = era, r = reg, s = "all" if not len(selplot) else ".".join(selplot))
            submitcomm = slurmscaff.format(extraS = extraslurm,
                                          nth     = nthreads,
                                          queue   = queue,
                                          jobname = jobname_,
                                          logpath = logpath.format(y = year, p = prod),
                                          command = PlottingCommand(prod, year, nthreads, outpath, selplot, reg, ratiorange, extra, useFibre, doUncs, doBlind, era, CMSText))
            print("Command:", submitcomm)
            if not pretend: os.system(submitcomm)
    else:
        for reg in region.split(","):
            execcomm = PlottingCommand(prod, year, nthreads, outpath, selplot, reg, ratiorange, extra, useFibre, doUncs, doBlind, era, CMSText)
            print("Command:", execcomm)
            if not pretend: os.system(execcomm)


    return



def PlottingCommand(prod, year, nthreads, outpath, selplot, region, ratio, extra, useFibre, doUncs, doBlind, era, CMSText):
    thecomm = commandscaff if not doBlind else commandblindscaff
    uncFile = "tw-run3/uncs-tw.txt"
    mcafile_    = "tw-run3/mca-tw.txt"
    if era != str(year):
        mcafile_ = "tw-run3/mca-tw-{e}.txt".format(e = era)
    cutsfile_   = "tw-run3/cuts-tw-{reg}.txt".format(reg = region if ("_" not in region) else region.split("_")[0] if ("differential" not in region and "cut" not in region) else region)
    plotsfile_  = "tw-run3/plots-tw/plots-tw-{reg}.txt".format(reg = region.replace("SF", "").replace("-ee", "").replace("-mm", "").replace("1j0t","1j1t").replace("2j0t","2j1t").replace("3j2t","2j2t").replace("3j1t","2j2t"))

    extraBlankSpaces = ""
    if CMSText == "Supplementary":
        extraBlankSpaces = "    "

    if "trigger" in region:
        mcafile_    = "tw-run3/triggerSFs/mca-trig.txt"
        uncFile = "tw-run3/triggerSFs/uncs-trig.txt"
        cutsfile_   = "tw-run3/triggerSFs/cuts-tw-{reg}.txt".format(reg = region)
        plotsfile_ = "tw-run3/triggerSFs/plots-trigger/plots-tw-{reg}.txt".format(reg = region if ("_" not in region) else region.split("_")[0])
        thecomm = commandtriggscaff
        
    
    samplespaths_ = "-P " + friendspath + "/" + prod + ("/" + year) * (year != "run3")
    if True: samplespaths_ = samplespaths_.replace("phedexrw", "phedex").replace("cienciasrw", "ciencias")

    nth_        = "" if nthreads == 0 else ("--split-factor=-1 -j " + str(nthreads))

    #friends_    = friendsscaff + (" --Fs {P}/5_mvas --Fs {P}/6_mvas_new" * ("MVA" in region))
    friends_    = friendsscaff + (" --Fs {P}/6_mvasRF" * ("MVA" in region))
    outpath_    = outpath + "/" + era + "/" + (region if "_" not in region else (region.split("_")[0] + "/" + region.split("_")[1]))
    selplot_    = " ".join( [ "--sP {p}".format(p = sp) for sp in selplot ] ) if len(selplot) else ""
    ratio_      = "--maxRatioRange " + ratio
    # Set the name of the region
    if "-ee" in region:
        nameregion_ = extraBlankSpaces + "    e^{#pm}e^{#mp}"
    elif "-mm" in region:
        nameregion_ = extraBlankSpaces + "    #mu^{#pm}#mu^{#mp}"
    else:
        nameregion_ = extraBlankSpaces + "    e^{#pm}#mu^{#mp}"
    if "_" not in region and "nojets" not in region:
        nameregion_ += " (" + region.replace("t", "b").replace("plus", "+") + ")"
    elif "differential" in region and "nojets" not in region:
        #nameregion_ += " (" + region.split("_")[0].replace("t", "b").replace("plus", "+") + "+0j_{#scale[1.2]{#font[12]{l}}}" + ")"
        nameregion_ += " (" + region.split("_")[0].replace("t", "b").replace("plus", "+") + "+0j_{#scale[1.2]{l}}" + ")"
    elif "MVA" in region and "nojets" not in region:
        nameregion_ += " (" + region.split("_")[0].replace("t", "b").replace("plus", "+") + ")"
    
    
    # Select luminosity based on era and year
    if era != str(year):
        luminosities = lumidict_byEra[era]
    else:
        luminosities = lumidict[year] if year != "run3" else str(lumidict["2022"]) + "," + str(lumidict["2022PostEE"])
    
    
    comm = thecomm.format(outpath      = outpath_,
                          friends      = friends_,
                          samplespaths = samplespaths_,
                          lumi         = luminosities,
                          nth          = nth_,
                          year         = year if year != "run3" else "2022,2022PostEE",
                          selplot      = selplot_,
                          mcafile      = mcafile_,
                          cutsfile     = cutsfile_,
                          plotsfile    = plotsfile_,
                          ratio        = ratio_,
                          extra        = ("--unc {uncFile} ").format(uncFile = uncFile) * doUncs + extra,
                          nameregion   = nameregion_,
                          cmsText      = CMSText)

    return comm


def confirm(message = "Do you wish to continue?"):
    """
    Ask user to enter y(es) or n(o) (case-insensitive).
    :return: True if the answer is Y.
    :rtype: bool
    """
    answer = "y"
    while answer not in ["y", "n", "yes", "no"]:
        answer = input(message + " [Y/N]\n").lower()
    return answer[0] == "y"



if __name__=="__main__":
    parser = argparse.ArgumentParser(usage = "python3 plotterHelper.py [options]", description = "Helper for plotting.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--production','-P', metavar = "prod",      dest = "prod",    required = True)
    parser.add_argument('--year',      '-y', metavar = 'year',      dest = "year",    required = False, default = "2022")
    parser.add_argument('--queue',     '-q', metavar = 'queue',     dest = "queue",   required = False, default = "")
    parser.add_argument('--extraArgs', '-e', metavar = 'extra',     dest = "extra",   required = False, default = "")
    parser.add_argument('--extraSlurmArgs','-eS',metavar='extraslurm',dest="extraslurm",required=False, default = "")
    parser.add_argument('--nthreads',  '-j', metavar = 'nthreads',  dest = "nthreads",required = False, default = 0, type = int)
    parser.add_argument('--select-plot','--sP',action = "append",   dest = "selplot", required = False, default = [])
    parser.add_argument('--pretend',   '-p', action = "store_true", dest = "pretend", required = False, default = False)
    parser.add_argument('--outpath',   '-o', metavar = 'outpath',   dest = "outpath", required = False, default = "./temp/varplots")
    parser.add_argument('--region',    '-r', metavar = 'region',    dest = "region",  required = False, default = "1j1t")
    parser.add_argument('--maxRatioRange', "-R", metavar = 'ratiorange', dest = "ratiorange", required = False, default = "0.8 1.2")
    parser.add_argument('--useFibre',  "-f", action = "store_true", dest = "useFibre", required = False, default = False)
    parser.add_argument('--uncertainties', "-u", action = "store_true", dest = "doUncs", required = False, default = False)
    parser.add_argument('--blind', action = "store_true", dest = "blindornot", required = False, default = False)
    parser.add_argument('--splitEras', action = "store_true", dest = "splitEras", required = False, default = False)
    parser.add_argument('--normaliseToData', action = "store_true", dest = "normToData", required = False, default = False)
    parser.add_argument('--cmsText',      '-cT', metavar = 'cmsText',      dest = "cmsText",    required = False, default = "Preliminary")



    args     = parser.parse_args()
    prod     = args.prod
    year     = args.year
    queue    = args.queue
    extra    = args.extra
    extraslurm=args.extraslurm
    nthreads = args.nthreads
    selplot  = args.selplot
    pretend  = args.pretend
    outpath  = args.outpath
    region   = args.region
    useFibre = args.useFibre
    ratiorange = args.ratiorange
    doUncs = args.doUncs
    normToData = args.normToData
    cmsText = args.cmsText
    

    if args.blindornot: extra += " --xp data"
    doBlind = args.blindornot

    Eras = [] # If era is empty, it will plot all the eras
    if args.splitEras:
        for era in lumidict_byEra:
            Eras.append(era)

    # Add the string for normalising to data
    if normToData:
        extra += " --scaleBkgSigToData ttbar --scaleBkgSigToData dy --scaleBkgSigToData nonworz --scaleBkgSigToData vvttv --scaleBkgSigToData tw"
    
    if queue != "":
        print("> Plotting jobs will be sent to the cluster.")
        if year == "all":
            print("   - All years and the combination will be plotted.")
            cont = False
            if   pretend:
                cont = True
            elif confirm("Three jobs per requested region will be sent to queue {q} with {j} requested threads to plot in each year and in the combination {pls}. Do you want to continue?".format(q = queue, j = nthreads, pls = "all the plots" if not len(selplot) else " and ".join(selplot))):
                cont = True

            if cont:
                for y in ["2022", "2022PostEE", "run3"]:
                    GeneralExecutioner( (prod, y, nthreads, outpath, selplot, region, ratiorange, queue, extra, pretend, useFibre, doUncs, doBlind, extraslurm, cmsText) )
        else:
            cont = False
            if   pretend:
                cont = True
            elif confirm("One job per requested region will be sent to queue {q} with {j} requested threads to plot in year {y} {pls}. Do you want to continue?".format(q = queue, j = nthreads, y = year, pls = "all the plots" if not len(selplot) else " and ".join(selplot))):
                cont = True

            if cont:
                if len(Eras) > 0:
                    for era in Eras:
                        GeneralExecutioner( (prod, year, nthreads, outpath, selplot, region, ratiorange, queue, extra, pretend, useFibre, doUncs, doBlind, extraslurm, era, cmsText) )
                else:
                    GeneralExecutioner( (prod, year, nthreads, outpath, selplot, region, ratiorange, queue, extra, pretend, useFibre, doUncs, doBlind, extraslurm, None, cmsText) )
    else:
        print("> Local execution chosen.")
        if year == "all":
            print("   - All years and the combination will be plotted.")
            for y in ["2022", "2022PostEE", "run3"]:
                GeneralExecutioner( (prod, y, nthreads, outpath, selplot, region, ratiorange, queue, extra, pretend, useFibre, doUncs, doBlind, extraslurm, cmsText) )
        else:
            if len(Eras) > 0:
                for era in Eras:
                    GeneralExecutioner( (prod, year, nthreads, outpath, selplot, region, ratiorange, queue, extra, pretend, useFibre, doUncs, doBlind, extraslurm, era, cmsText) )
            else:
                GeneralExecutioner( (prod, year, nthreads, outpath, selplot, region, ratiorange, queue, extra, pretend, useFibre, doUncs, doBlind, extraslurm, None, cmsText) )
