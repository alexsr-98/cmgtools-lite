import ROOT as r
import warnings as wr
import sys, os, argparse, math
from array import array
from copy import deepcopy
from multiprocessing import Pool

sys.path.append('{cmsswpath}/src/CMGTools/TTHAnalysis/python/plotter/twttbar-run2UL/differential/'.format(cmsswpath = os.environ['CMSSW_BASE']))
import varList as vl

r.gROOT.SetBatch(True)
vl.SetUpWarnings()

friendspath   = "/lustrefs/hdd_pool_dir/nanoAODv9/bb4l/productions"
logpath       = friendspath + "/{p}/{y}/logs/cards_differential"
friendsscaff = "--FMCs {P}/0_jecs --Fs {P}/1_lepsuncsAndParticle --Fs {P}/2_cleaning --Fs {P}/3_varstrigger --FMCs {P}/4_scalefactors"

slurmscaff    = "sbatch -c {nth} -p {queue} -J {jobname} -e {logpath}/log.%j.%x.err -o {logpath}/log.%j.%x.out --wrap '{command}'"

commandscaff  = '''python3 makeShapeCards_TopRun2.py --tree NanoAOD {mcafile} {cutsfile} "{variable}" "{bins}" {samplespaths} {friends} --od {outpath} -l {lumi} {nth} -f -L {func} --neg -W "{w}" --year {year} {asimovornot} {uncs} {extra} {name} --AP --storeAll --notVarsChanges --threshold 0.001'''

theweights    = 'MuonIDSF * MuonISOSF * ElecIDSF * ElecRECOSF * TrigSF * puWeight * bTagWeight * L1PreFiringWeight_Nom * jetPUidWeight'
minchunkbytes = 1000
redofiles     = True
#redofiles     = False

#### IMPORTANT NOTE: it is crucial to keep the xsec and lumi weights, although afterwards we
# will have to divide by the second one. The first one is relevant because of the dividing of the files.
# The second one, because of the proportions between MC simulations.


def PythonListToString(theL):
    ret = "["
    for el in theL:
        ret += str(el)
        if el != theL[-1]:
            ret += ", "
    ret += "]"
    return ret


def ExecuteOrSubmitTask(tsk):
    prod, year, variable, asimov, nthreads, outpath, noUnc, useFibre, extra, pretend, queue = tsk
    if queue == "":
        thecomm = CardsCommand(prod, year, variable, asimov, nthreads, outpath, noUnc, useFibre, extra)
        print("Command: " + thecomm)

        if not pretend:
            os.system(thecomm)

    else:
        if not os.path.isdir(logpath.format(p = prod, y = yr)):
            os.system("mkdir -p " + logpath.format(p = prod, y = yr))

        thecomm = slurmscaff.format(nth     = nthreads,
                                    queue   = queue,
                                    jobname = "CMGTcardsForDiff_" + variable + "_" + year,
                                    logpath = logpath.format(p = prod, y = yr),
                                    command = CardsCommand(prod, year, variable, asimov, nthreads, outpath, noUnc, useFibre, extra))

        print("Command: " + thecomm)
        if not pretend:
            os.system(thecomm)

    return


def produceAdHocMCAandUncsfiles(thepath, they, thev):
    outmain = ""
    with open("./twttbar-run2UL/differential/mca/mca_scaffold.txt", "r") as scaff:
        for line in scaff:
            if "TOREPLACE" in line:
                outmain += line.replace("TOREPLACEAPV2016", thepath + "/mca-bb4l-2016apv.txt").replace("TOREPLACE2016", thepath + "/mca-bb4l-2016.txt").replace("TOREPLACE2017", thepath + "/mca-bb4l-2017.txt").replace("TOREPLACE2018", thepath + "/mca-bb4l-2018.txt").replace("TOREPLACESYS2016", thepath + "/mca-bb4l-2016-uncs.txt").replace("TOREPLACESYS2017", thepath + "/mca-bb4l-2017-uncs.txt").replace("TOREPLACESYS2018", thepath + "/mca-bb4l-2018-uncs.txt")
            else:
                outmain += line

    foutmain = open(thepath + "/" + "mca-total.txt", "w")
    foutmain.write(outmain)
    foutmain.close(); del foutmain

    cutstring = "(Dresschannel == 1) && (nDressSelLep == 2) && (DressLep1_Pt > 25) && (DressLep2_Pt > 20) && (DressisSS == 0) && (DressminMllAFAS > 20) && (nDressSelJet == 2) && (nDressBSelJet == 2)"
    pathtooriginalmca = "./twttbar-run2UL/differential/mca/"
    if they == "run2":
        for iY in ["2016apv", "2016", "2017", "2018"]:
            outyear = ""
            with open(pathtooriginalmca + "/mca-{y}-bb4l.txt".format(y = iY), "r") as scaff:
                for line in scaff:
                    tmpline = line.strip()
                    if tmpline == "":       continue
                    elif tmpline[0] == "#": continue

                    for iB in range(len(vl.varList[thev]["bins_particle"]) - 1):
                        tmpcut = "(" + str(vl.varList[thev]["bins_particle"][iB]) + "<=" + vl.varList[thev]["var_particle"] + ") && (" + vl.varList[thev]["var_particle"] + "<" + str(vl.varList[thev]["bins_particle"][iB + 1]) + ")"
                        if "Origchannel" not in tmpline:
                            outyear += tmpline.replace(";", ": " + cutstring + " && " + tmpcut + " ;").replace("bb4l+", "bb4l_partbin{b}+".format(b = iB)).replace("FillColor=798", "FillColor=" + str(vl.coloursForDiffWithCombine[iB + 1]) + ",Label='bb4l bin {b}'".format(b = iB + 1)) + "\n"
                        else:
                            outyear += tmpline.replace(";", "&& " + cutstring + " && " + tmpcut + " ;").replace("bb4l+", "bb4l_partbin{b}+".format(b = iB)).replace("FillColor=798", "FillColor=" + str(vl.coloursForDiffWithCombine[iB + 1]) + ",Label='bb4l bin {b}'".format(b = iB + 1)) + "\n"

                    outyear += tmpline.replace(";", ": (!(" + cutstring + ")) ;").replace("bb4l+", "bb4l_bkg").replace("FillColor=798", "FillColor=" + str(vl.coloursForDiffWithCombine[0]) + ",Label='bb4l bkg.'") + "\n"

            foutyear = open(thepath + "/" + "mca-bb4l-{y}.txt".format(y = iY), "w")
            foutyear.write(outyear)
            foutyear.close(); del foutyear

            # now systs
            outyear = ""
            #with open(pathtooriginalmca + "/mca-{y}-bb4l-uncs.txt".format(y = iY), "r") as scaff:
                #for line in scaff:
                    #tmpline = line.strip()
                    #if tmpline == "":       continue
                    #elif tmpline[0] == "#": continue

                    #for iB in range(len(vl.varList[thev]["bins_particle"]) - 1):
                        #tmpcut = "(" + str(vl.varList[thev]["bins_particle"][iB]) + "<=" + vl.varList[thev]["var_particle"] + ") && (" + vl.varList[thev]["var_particle"] + "<" + str(vl.varList[thev]["bins_particle"][iB + 1]) + ")"
                        #outyear += tmpline.replace(";", ": " + cutstring + " && " + tmpcut + " ;").replace(tmpline.split(":")[0].replace(" ", ""), tmpline.split(":")[0].replace(" ", "") + "_partbin{b}".format(b = iB)) + "\n"

                    #outyear += tmpline.replace(";", ": (!(" + cutstring + ")) ;").replace(tmpline.split(":")[0].replace(" ", ""), tmpline.split(":")[0].replace(" ", "") + "_bkg") + "\n"

            #foutyear = open(thepath + "/" + "mca-bb4l-{y}-uncs.txt".format(y = iY), "w")
            #foutyear.write(outyear)
            #foutyear.close(); del foutyear
    else:
        outyear = ""
        with open(pathtooriginalmca + "/mca-{y}-bb4l.txt".format(y = they), "r") as scaff:
            for line in scaff:
                tmpline = line.strip()
                if tmpline == "":       continue
                elif tmpline[0] == "#": continue

                for iB in range(len(vl.varList[thev]["bins_particle"]) - 1):
                    tmpcut = "(" + str(vl.varList[thev]["bins_particle"][iB]) + "<=" + vl.varList[thev]["var_particle"] + ") && (" + vl.varList[thev]["var_particle"] + "<" + str(vl.varList[thev]["bins_particle"][iB + 1]) + ")"
                    if "Origchannel" not in tmpline:
                        outyear += tmpline.replace(";", ": " + cutstring + " && " + tmpcut + " ;").replace("bb4l+", "bb4l_partbin{b}+".format(b = iB)).replace("FillColor=798", "FillColor=" + str(vl.coloursForDiffWithCombine[iB + 1]) + ",Label='bb4l bin {b}'".format(b = iB + 1)) + "\n"
                    else:
                        outyear += tmpline.replace(";", "&& " + cutstring + " && " + tmpcut + " ;").replace("bb4l+", "bb4l_partbin{b}+".format(b = iB)).replace("FillColor=798", "FillColor=" + str(vl.coloursForDiffWithCombine[iB + 1]) + ",Label='bb4l bin {b}'".format(b = iB + 1)) + "\n"

                outyear += tmpline.replace(";", ": (!(" + cutstring + ")) ;").replace("bb4l+", "bb4l_bkg").replace("FillColor=798", "FillColor=" + str(vl.coloursForDiffWithCombine[0]) + ",Label='bb4l bkg.'") + "\n"

            foutyear = open(thepath + "/" + "mca-bb4l-{y}.txt".format(y = they), "w")
            foutyear.write(outyear)
            foutyear.close(); del foutyear
        ## now systs
        #outyear = ""
        #with open(pathtooriginalmca + "/mca-{y}-bb4l-uncs.txt".format(y = they), "r") as scaff:
            #for line in scaff:
                #tmpline = line.strip()
                #if tmpline == "":       continue
                #elif tmpline[0] == "#": continue

                #for iB in range(len(vl.varList[thev]["bins_particle"]) - 1):
                    #tmpcut = "(" + str(vl.varList[thev]["bins_particle"][iB]) + "<=" + vl.varList[thev]["var_particle"] + ") && (" + vl.varList[thev]["var_particle"] + "<" + str(vl.varList[thev]["bins_particle"][iB + 1]) + ")"
                    #outyear += tmpline.replace(";", ": " + cutstring + " && " + tmpcut + " ;").replace(tmpline.split(":")[0].replace(" ", ""), tmpline.split(":")[0].replace(" ", "") + "_partbin{b}".format(b = iB)) + "\n"

                #outyear += tmpline.replace(";", ": (!(" + cutstring + ")) ;").replace(tmpline.split(":")[0].replace(" ", ""), tmpline.split(":")[0].replace(" ", "") + "_bkg") + "\n"

        #foutyear = open(thepath + "/" + "mca-bb4l-{y}-uncs.txt".format(y = they), "w")
        #foutyear.write(outyear)
        #foutyear.close(); del foutyear
    
    
    ################################# INCERTIDUMBRES
    outuncs = ""
    with open("./twttbar-run2UL/differential/uncs/uncs-twttbar.txt", "r") as scaff:
        for line in scaff:
            if "alias" in line:
                outuncs += line.replace("bb4l", "|".join(["bb4l_bkg"] + ["bb4l_partbin{b}".format(b = iB) for iB in range(len(vl.varList[thev]["bins_particle"]) - 1)]))
            else:
                tmpline = line.strip()
                if tmpline.split(":")[0].replace(" ", "") in vl.VariatedSamplesUncertaintySources and "$bb4lProc" in tmpline.split(":")[1]:
                    stringtoreplace = None
                    if "EstimateFromXbins" in tmpline:
                        stringtoreplace = "EstimateFromXbins=1"
                    elif "fitOrder=" in tmpline and "symmAfterFit" not in tmpline:
                        stringtoreplace = "fitOrder=" + tmpline.split("fitOrder=")[1][0]
                    elif "fitOrder=" in tmpline:
                        stringtoreplace = ["fitOrder=" + tmpline.split("fitOrder=")[1][0], "symmAfterFit=1"]
                    

                    for iB in range(len(vl.varList[thev]["bins_particle"]) - 1):
                        twotimestmpline = tmpline.replace("$bb4lProc", "bb4l_partbin{b}".format(b = iB)) + "\n"
                        for substr in tmpline.split(";")[0].split(":"):
                            if "syst" in substr:
                                tmpsubstr = substr.replace("'", "").replace(" ", "").replace("[", "").replace("]", "")
                                if "\," in tmpsubstr:
                                    for secsubstr in tmpsubstr.split("\,"):
                                        twotimestmpline = twotimestmpline.replace(secsubstr, secsubstr + "_partbin{b}".format(b = iB))
                                else:
                                    twotimestmpline = twotimestmpline.replace(tmpsubstr, tmpsubstr + "_partbin{b}".format(b = iB))
                        if (("EstimateFromXbins" in tmpline) or ("fitOrder" in tmpline)) and ("particle_smoothing" in vl.varList[thev]):
                            theunc = tmpline.split(":")[0].replace(" ", "")
                            if theunc in vl.varList[thev]["particle_smoothing"]:
                                if isinstance(vl.varList[thev]["particle_smoothing"][theunc], dict): ### Fine tuned depending on bin
                                    if iB in vl.varList[thev]["particle_smoothing"][theunc]:
                                        if "NOM" != vl.varList[thev]["particle_smoothing"][theunc][iB]:
                                            if isinstance(stringtoreplace, list):
                                                twotimestmpline = twotimestmpline.replace(stringtoreplace[0], vl.varList[thev]["particle_smoothing"][theunc][iB])
                                                for k in range(1, len(stringtoreplace)): twotimestmpline = twotimestmpline.replace(stringtoreplace[k], "")
                                            else:
                                                twotimestmpline = twotimestmpline.replace(stringtoreplace, vl.varList[thev]["particle_smoothing"][theunc][iB])
                                    elif -1 in vl.varList[thev]["particle_smoothing"][theunc]:
                                        if isinstance(stringtoreplace, list):
                                            twotimestmpline = twotimestmpline.replace(stringtoreplace[0], vl.varList[thev]["particle_smoothing"][theunc][-1])
                                            for k in range(1, len(stringtoreplace)): twotimestmpline = twotimestmpline.replace(stringtoreplace[k], "")
                                        else:
                                            twotimestmpline = twotimestmpline.replace(stringtoreplace, vl.varList[thev]["particle_smoothing"][theunc][-1])
                                else: ### The same for all non-bkg
                                    if isinstance(stringtoreplace, list):
                                        twotimestmpline = twotimestmpline.replace(stringtoreplace[0], vl.varList[thev]["particle_smoothing"][theunc])
                                        for k in range(1, len(stringtoreplace)): twotimestmpline = twotimestmpline.replace(stringtoreplace[k], "")
                                    else:
                                        twotimestmpline = twotimestmpline.replace(stringtoreplace, vl.varList[thev]["particle_smoothing"][theunc])
                                        

                        outuncs += twotimestmpline
                    
                    # Now, we add the background
                    twotimestmpline = tmpline.replace("$bb4lProc", "bb4l_bkg") + "\n"
                    for substr in tmpline.split(";")[0].split(":"):
                        if "syst" in substr:
                            tmpsubstr = substr.replace("'", "").replace(" ", "").replace("[", "").replace("]", "")
                            if "\," in tmpsubstr:
                                for secsubstr in tmpsubstr.split("\,"):
                                    twotimestmpline = twotimestmpline.replace(secsubstr, secsubstr + "_bkg")
                            else:
                                twotimestmpline = twotimestmpline.replace(tmpsubstr, tmpsubstr + "_bkg")
                    
                    outuncs += twotimestmpline
                    


                elif "pdf" in tmpline.split(":")[0] and "syst" in tmpline.split(":")[1] and "altSamplePDF" not in tmpline.split(":")[3]:
                    for iB in range(len(vl.varList[thev]["bins_particle"]) - 1):
                        outuncs += tmpline.replace("syst_tw_pdf_2016", "syst_tw_pdf_2016_partbin{b}".format(b = iB)).replace("syst_tw_pdf_2017", "syst_tw_pdf_2017_partbin{b}".format(b = iB)).replace("syst_tw_pdf_2018", "syst_tw_pdf_2018_partbin{b}".format(b = iB)) + "\n"
                    outuncs += tmpline.replace("syst_tw_pdf_2016", "syst_tw_pdf_2016_bkg").replace("syst_tw_pdf_2017", "syst_tw_pdf_2017_bkg").replace("syst_tw_pdf_2018", "syst_tw_pdf_2018_bkg") + "\n"
                    
                elif "pdf" in tmpline.split(":")[0] and "altSamplePDF" in tmpline.split(":")[3]:
                    for iB in range(len(vl.varList[thev]["bins_particle"]) - 1):
                        outuncs += tmpline.replace("$bb4lProc", "bb4l_partbin{b}".format(b = iB)).replace("syst_tw_pdf_2016", "syst_tw_pdf_2016_partbin{b}".format(b = iB)).replace("syst_tw_pdf_2017", "syst_tw_pdf_2017_partbin{b}".format(b = iB)).replace("syst_tw_pdf_2018", "syst_tw_pdf_2018_partbin{b}".format(b = iB)) + "\n"
                    outuncs += tmpline.replace("$bb4lProc", "bb4l_bkg").replace("syst_tw_pdf_2016", "syst_tw_pdf_2016_bkg").replace("syst_tw_pdf_2017", "syst_tw_pdf_2017_bkg").replace("syst_tw_pdf_2018", "syst_tw_pdf_2018_bkg") + "\n"
                    
                else:
                    tmpline = line.strip()
                    stringtoreplace = None
                    if "EstimateFromXbins" in tmpline:
                        stringtoreplace = "EstimateFromXbins=1"
                    elif "fitOrder=" in tmpline and "symmAfterFit" not in tmpline:
                        stringtoreplace = "fitOrder=" + tmpline.split("fitOrder=")[1][0]
                    elif "fitOrder=" in tmpline:
                        stringtoreplace = ["fitOrder=" + tmpline.split("fitOrder=")[1][0], "symmAfterFit=1"]

                    if "$bb4lProc" in tmpline and stringtoreplace and ("particle_smoothing" in vl.varList[thev]):
                        theunc = tmpline.split(":")[0].replace(" ", "")
                        if theunc in vl.varList[thev]["particle_smoothing"]:
                            for iB in range(len(vl.varList[thev]["bins_particle"]) - 1):
                                twotimestmpline = tmpline.replace("$bb4lProc", "bb4l_partbin{b}".format(b = iB)) + "\n"
                                if isinstance(vl.varList[thev]["particle_smoothing"][theunc], dict): ### Fine tuned depending on bin
                                    if iB in vl.varList[thev]["particle_smoothing"][theunc]:
                                        if "NOM" != vl.varList[thev]["particle_smoothing"][theunc][iB]:
                                            if isinstance(stringtoreplace, list):
                                                twotimestmpline = twotimestmpline.replace(stringtoreplace[0], vl.varList[thev]["particle_smoothing"][theunc][iB])
                                                for k in range(1, len(stringtoreplace)): twotimestmpline = twotimestmpline.replace(stringtoreplace[k], "")
                                            else:
                                                twotimestmpline = twotimestmpline.replace(stringtoreplace, vl.varList[thev]["particle_smoothing"][theunc][iB])
                                    elif -1 in vl.varList[thev]["particle_smoothing"][theunc]:
                                        if isinstance(stringtoreplace, list):
                                            twotimestmpline = twotimestmpline.replace(stringtoreplace[0], vl.varList[thev]["particle_smoothing"][theunc][-1])
                                            for k in range(1, len(stringtoreplace)): twotimestmpline = twotimestmpline.replace(stringtoreplace[k], "")
                                        else:
                                            twotimestmpline = twotimestmpline.replace(stringtoreplace, vl.varList[thev]["particle_smoothing"][theunc][-1])
                                else: ### The same for all non-bkg
                                    if isinstance(stringtoreplace, list):
                                        twotimestmpline = twotimestmpline.replace(stringtoreplace[0], vl.varList[thev]["particle_smoothing"][theunc])
                                        for k in range(1, len(stringtoreplace)): twotimestmpline = twotimestmpline.replace(stringtoreplace[k], "")
                                    else:
                                        twotimestmpline = twotimestmpline.replace(stringtoreplace, vl.varList[thev]["particle_smoothing"][theunc])
                                
                                outuncs += twotimestmpline
                            
                            twotimestmpline = tmpline.replace("$bb4lProc", "bb4l_bkg") + "\n"
                            outuncs += twotimestmpline
                        else:
                            outuncs += tmpline + "\n"
                    else:
                        outuncs += tmpline + "\n"
                    
    foutuncs = open(thepath + "/" + "twttbar-uncs-adhoc.txt", "w")
    foutuncs.write(outuncs)
    foutuncs.close(); del foutuncs
    return


def CardsCommand(prod, year, var, isAsimov, nthreads, outpath, noUnc, useFibre, extra_):
    cutsfile_  = "twttbar-run2UL/differential/cuts/cuts-detector-2j2t.txt"

    samplespaths_ = "-P " + friendspath + "/" + prod + ("/" + year) * (year != "run2")
    if useFibre: samplespaths_ = samplespaths_.replace("phedexrw", "phedex").replace("cienciasrw", "ciencias")

    nth_       = "" if nthreads == 0 else ("--split-factor=-1 -j " + str(nthreads))
    friends_   = friendsscaff
    outpath_   = outpath + "/" + year + "/" + var + "/sigextr_fit_combine/bincards"

    if not os.path.isdir(outpath_):
        os.system("mkdir -p " + outpath_)

    mcafile_   = outpath_ + "/" + "mca-total.txt"
    if not os.path.isfile(mcafile_) or redofiles:
        produceAdHocMCAandUncsfiles(outpath_, year, var)

    variable_  = vl.varList[var]["var_detector"]
    bins_      = PythonListToString(vl.varList[var]["bins_detector"])
    name_      = "--binname card_" + var

    theweights_ = theweights

    comm = commandscaff.format(outpath      = outpath_,
                               friends      = friends_,
                               samplespaths = samplespaths_,
                               lumi         = (vl.LumiDict[year] if year != "run2" else str(vl.LumiDict["2016apv"]) + "," + str(vl.LumiDict["2016"]) + "," + str(vl.LumiDict["2017"]) + "," + str(vl.LumiDict["2018"])),
                               variable     = variable_,
                               bins         = bins_,
                               nth          = nth_,
                               year         = year if year != "run2" else "2016apv,2016,2017,2018",
                               asimovornot  = "--asimov s+b" if isAsimov else "",
                               mcafile      = mcafile_,
                               cutsfile     = cutsfile_,
                               uncs         = "--unc {adhocu} --amc".format(adhocu = outpath_ + "/twttbar-uncs-adhoc.txt") if (not noUnc) else "--amc",
                               name         = name_,
                               func         = "twttbar-run2UL/functions_twttbar.cc",
                               extra        = extra_,
                               w            = theweights_,)
    return comm


def CheckProducedCardsByTask(task):
    prod, year, variable, asimov, nthreads, outpath, noUnc, useFibre, extra, pretend, queue, thebin = task

    ch = "card_bin" + str(thebin) + ".root"

    chkpath = outpath + "/" + year + "/" + variable + "/sigextr_fit_combine/bincards/" + ch

    #print chkpath

    if not os.path.isfile(chkpath):                                 #### 1st: existance
        print("# Card {chk} has not been found.".format(chk = ch))
        return False
    elif os.path.getsize(chkpath) <= minchunkbytes:                 #### 2nd: size
        print("# Card {chk} has less size than the minimum.".format(chk = ch))
        return False
    else:
        fch = r.TFile.Open(chkpath, "READ")
        if not fch:                                                 #### 3rd: ROOT access (corruption)
            print("# Card {chk} cannot be accessed: it is corrupted.".format(chk = ch))
            return False

        fch.Close(); del fch
    return True


if __name__=="__main__":
    parser = argparse.ArgumentParser(usage = "python nanoAOD_checker.py [options]", description = "Checker tool for the outputs of nanoAOD production (NOT postprocessing)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--year',      '-y', metavar = 'year',       dest = "year",     required = False, default = "all")
    parser.add_argument('--variable',  '-v', metavar = 'variable',   dest = "variable", required = False, default = "all")
    parser.add_argument('--step',      '-s', metavar = 'step',       dest = "step",     required = False,  default = 0, type = int)
    parser.add_argument('--extraArgs', '-e', metavar = 'extra',      dest = "extra",    required = False, default = "")
    parser.add_argument('--nthreads',  '-j', metavar = 'nthreads',   dest = "nthreads", required = False, default = 0, type = int)
    parser.add_argument('--pretend',   '-p', action  = "store_true", dest = "pretend",  required = False, default = False)
    parser.add_argument('--production','-P', metavar = "prod",       dest = "prod",     required = True)
    parser.add_argument('--queue',     '-q', metavar = 'queue',      dest = "queue",    required = False, default = "")
    parser.add_argument('--outpath',   '-o', metavar = 'outpath',    dest = "outpath",  required = False, default = "./temp/differential")
    parser.add_argument('--nounc',     '-nu',action  = "store_true", dest = "nounc",    required = False, default = False)
    parser.add_argument('--useFibre',  '-f', action  = "store_true", dest = "useFibre", required = False, default = False)
    parser.add_argument('--check',     '-c', action  = "store_true", dest = "check",    required = False, default = False)
    parser.add_argument('--asimov',    '-a', action  = "store_true", dest = "asimov",   required = False, default = False)


    args     = parser.parse_args()
    prod     = args.prod
    step     = args.step
    queue    = args.queue
    year     = args.year
    check    = args.check
    extra    = args.extra
    nthreads = args.nthreads
    outpath  = args.outpath
    noUnc    = args.nounc
    asimov    = args.asimov
    pretend  = args.pretend
    variable = args.variable
    useFibre = args.useFibre
    thevars  = vl.varList["Names"]["Variables"]
    theyears = ["2016", "2017", "2018", "run2"]
    vetolist = ["plots", "Fiducial", "control", "tables"]

    if check:
        wr.warn("WARNING: check is only implemented for the first step.")

    if step == 0:
        if not check:
            print("> Preparing to submit the cards for the fit...")
        else:
            print("> Preparing to check all the produced cards...")

        tasks    = []

        if "," in variable:
            thevars = variable.split(",")
        elif variable.lower() != "all":
            thevars = [ variable ]

        if year.lower() != "all":
            if "," in year:
                theyears = year.split(",")
            else:
                theyears = [ year ]

        for yr in theyears:
            for var in thevars:
                if "Fiducial" in var: continue
                tasks.append( (prod, yr, var, asimov, nthreads, outpath, noUnc, useFibre, extra, pretend, queue) )

        #print tasks
        #for el in tasks: print el
        #sys.exit()

        if not check:
            if queue != "":
                print("> Executing...")
            else:
                print("> Launching jobs...")

            calculate = True
   #         calculate = False
            for task in tasks:
                print("    - Processing " + str(task))
  #              if str(task) == "('2021-06-09', 'run2', 'Lep1_Pt', True, 64, 'temp_2021_06_29_tolosplots/differential/', False, True, '', True, '', 0)":
                    #print "OJOCUIDAOOOOOOOOOOOOOO"
                    #continue
 #                   calculate = True

                if calculate:
                    ExecuteOrSubmitTask(task)
                    #sys.exit()
                    #calculate = False

        else:
            print("> Beginning the checks on produced cards...")
            taskstoredo = []
            for task in tasks:
                if not CheckProducedCardsByTask(task):
                    taskstoredo.append(task)

            if not len(taskstoredo):
                print("> Everything is ok!")
            else:
                print("> The following cards should be remade:")
                for el in taskstoredo: print(el)

                if vl.confirm("\n> Do you want to redo these cards?"):
                    for el in taskstoredo: ExecuteOrSubmitTask(el)
