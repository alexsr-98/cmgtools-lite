import ROOT as r
import warnings as wr
import sys, os, argparse, math
from array import array
from copy import deepcopy
from multiprocessing import Pool

sys.path.append('{cmsswpath}/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/differential/'.format(cmsswpath = os.environ['CMSSW_BASE']))
import varList as vl

r.gROOT.SetBatch(True)
vl.SetUpWarnings()

friendspath   = "/lustrefs/hdd_pool_dir/nanoAODv11/tw-run3/productions"
logpath       = friendspath + "/{p}/{y}/logs/cards_differential"

friendsscaff = "--Fs {P}/0_jecs --Fs {P}/1_lepsuncsAndParticle --Fs {P}/2_cleaning --Fs {P}/3_varstrigger --FMCs {P}/4_scalefactors"

slurmscaff    = "sbatch -c {nth} -p {queue} -J {jobname} -e {logpath}/log.%j.%x.err -o {logpath}/log.%j.%x.out --wrap '{command}'"

commandscaff  = '''python3 makeShapeCards_TopRun2.py --tree NanoAOD {mcafile} {cutsfile} "{variable}" "{bins}" {samplespaths} {friends} --od {outpath} -l {lumi} {nth} -f -L {func} --neg -W "{w}" --year {year} {asimovornot} {uncs} {extra} {name} --AP --storeAll --notVarsChanges --threshold 0.001'''


theweights    = "MuonIDSF * MuonISOSF * ElecIDSF * ElecRECOSF * TrigSF * bTagWeight * puWeight"
minchunkbytes = 1000
redofiles     = True
#redofiles     = False

#### IMPORTANT NOTE: it is crucial to keep the xsec and lumi weights, although afterwards we
# will have to divide by the second one. The first one is relevant becasue of the dividing of the files.
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
    with open("./tw-run3/differential/mca_scaffold.txt", "r") as scaff:
        for line in scaff:
            if "TOREPLACE" in line:
                outmain += line.replace("TOREPLACE2022PostEE", thepath + "/mca-tw-2022PostEE.txt").replace("TOREPLACE2022", thepath + "/mca-tw-2022.txt").replace("TOREPLACESYS2022PostEE", thepath + "/mca-tw-2022PostEE-uncs.txt").replace("TOREPLACESYS2022", thepath + "/mca-tw-2022-uncs.txt")
            else:
                outmain += line

    foutmain = open(thepath + "/" + "mca-total.txt", "w")
    foutmain.write(outmain)
    foutmain.close(); del foutmain

    cutstring = "(Dresschannel == 1) && (nDressSelLep >= 2) && (DressLep1_Pt > 25) && (DressLep2_Pt > 20) && (DressisSS == 0) && (DressminMllAFAS > 20) && (nDressSelJet == 1) && (nDressBSelJet == 1) && (nDressSelLooseJet == 0)"
    pathtooriginalmca = "./tw-run3/mca-tw-includes/"
    if they == "run3":
        for iY in ["2022", "2022PostEE"]:
            outyear = ""
            with open(pathtooriginalmca + "/mca-{y}-tw.txt".format(y = iY), "r") as scaff:
                for line in scaff:
                    tmpline = line.strip()
                    if tmpline == "":       continue
                    elif tmpline[0] == "#": continue

                    for iB in range(len(vl.varList[thev]["bins_particle"]) - 1):
                        tmpcut = "(" + str(vl.varList[thev]["bins_particle"][iB]) + "<=" + vl.varList[thev]["var_particle"] + ") && (" + vl.varList[thev]["var_particle"] + "<" + str(vl.varList[thev]["bins_particle"][iB + 1]) + ")"
                        outyear += tmpline.replace(";", ": " + cutstring + " && " + tmpcut + " ;").replace("tw+", "tw_partbin{b}+".format(b = iB)).replace("FillColor=798", "FillColor=" + str(vl.coloursForDiffWithCombine[iB + 1]) + ",Label='tW bin {b}'".format(b = iB + 1)) + "\n"

                    outyear += tmpline.replace(";", ": (!(" + cutstring + ")) ;").replace("tw+", "tw_bkg").replace("FillColor=798", "FillColor=" + str(vl.coloursForDiffWithCombine[0]) + ",Label='tW bkg.'") + "\n"

            foutyear = open(thepath + "/" + "mca-tw-{y}.txt".format(y = iY), "w")
            foutyear.write(outyear)
            foutyear.close(); del foutyear

            # now systs
            outyear = ""
            with open(pathtooriginalmca + "/mca-{y}-tw-uncs.txt".format(y = iY), "r") as scaff:
                for line in scaff:
                    tmpline = line.strip()
                    if tmpline == "":       continue
                    elif tmpline[0] == "#": continue

                    for iB in range(len(vl.varList[thev]["bins_particle"]) - 1):
                        tmpcut = "(" + str(vl.varList[thev]["bins_particle"][iB]) + "<=" + vl.varList[thev]["var_particle"] + ") && (" + vl.varList[thev]["var_particle"] + "<" + str(vl.varList[thev]["bins_particle"][iB + 1]) + ")"
                        outyear += tmpline.replace(";", ": " + cutstring + " && " + tmpcut + " ;").replace(tmpline.split(":")[0].replace(" ", ""), tmpline.split(":")[0].replace(" ", "") + "_partbin{b}".format(b = iB)) + "\n"

                    outyear += tmpline.replace(";", ": (!(" + cutstring + ")) ;").replace(tmpline.split(":")[0].replace(" ", ""), tmpline.split(":")[0].replace(" ", "") + "_bkg") + "\n"

            foutyear = open(thepath + "/" + "mca-tw-{y}-uncs.txt".format(y = iY), "w")
            foutyear.write(outyear)
            foutyear.close(); del foutyear
    else:
        outyear = ""
        with open(pathtooriginalmca + "/mca-{y}-tw.txt".format(y = they), "r") as scaff:
            for line in scaff:
                tmpline = line.strip()
                if tmpline == "":       continue
                elif tmpline[0] == "#": continue

                for iB in range(len(vl.varList[thev]["bins_particle"]) - 1):
                    tmpcut = "(" + str(vl.varList[thev]["bins_particle"][iB]) + "<=" + vl.varList[thev]["var_particle"] + ") && (" + vl.varList[thev]["var_particle"] + "<" + str(vl.varList[thev]["bins_particle"][iB + 1]) + ")"
                    outyear += tmpline.replace(";", ": " + cutstring + " && " + tmpcut + " ;").replace("tw+", "tw_partbin{b}+".format(b = iB)).replace("FillColor=798", "FillColor=" + str(vl.coloursForDiffWithCombine[iB + 1]) + ",Label='tW bin {b}'".format(b = iB + 1)) + "\n"

                outyear += tmpline.replace(";", ": (!(" + cutstring + ")) ;").replace("tw+", "tw_bkg").replace("FillColor=798", "FillColor=" + str(vl.coloursForDiffWithCombine[0]) + ",Label='tW bkg.'") + "\n"

            foutyear = open(thepath + "/" + "mca-tw-{y}.txt".format(y = they), "w")
            foutyear.write(outyear)
            foutyear.close(); del foutyear
        # now systs
        outyear = ""
        with open(pathtooriginalmca + "/mca-{y}-tw-uncs.txt".format(y = they), "r") as scaff:
            for line in scaff:
                tmpline = line.strip()
                if tmpline == "":       continue
                elif tmpline[0] == "#": continue

                for iB in range(len(vl.varList[thev]["bins_particle"]) - 1):
                    tmpcut = "(" + str(vl.varList[thev]["bins_particle"][iB]) + "<=" + vl.varList[thev]["var_particle"] + ") && (" + vl.varList[thev]["var_particle"] + "<" + str(vl.varList[thev]["bins_particle"][iB + 1]) + ")"
                    outyear += tmpline.replace(";", ": " + cutstring + " && " + tmpcut + " ;").replace(tmpline.split(":")[0].replace(" ", ""), tmpline.split(":")[0].replace(" ", "") + "_partbin{b}".format(b = iB)) + "\n"

                outyear += tmpline.replace(";", ": (!(" + cutstring + ")) ;").replace(tmpline.split(":")[0].replace(" ", ""), tmpline.split(":")[0].replace(" ", "") + "_bkg") + "\n"

        foutyear = open(thepath + "/" + "mca-tw-{y}-uncs.txt".format(y = they), "w")
        foutyear.write(outyear)
        foutyear.close(); del foutyear
    
    
    ################################# INCERTIDUMBRES
    outuncs = ""
    # Sin Suavizado
    #
    # Con Suavizado
    #with open("./tw-run3/uncs-tw_1j1tdiff_" + thev.lower().replace("_", "") + ".txt", "r") as scaff:
    with open("./tw-run3/uncs-tw.txt", "r") as scaff:
        for line in scaff:
            if "alias" in line:
                outuncs += line.replace("tw", "|".join(["tw_bkg"] + ["tw_partbin{b}".format(b = iB) for iB in range(len(vl.varList[thev]["bins_particle"]) - 1)]))
            else:
                tmpline = line.strip()
                if tmpline.split(":")[0].replace(" ", "") in vl.VariatedSamplesUncertaintySources and "$tWProc" in tmpline.split(":")[1]:
                    stringtoreplace = None
                    if "EstimateFromXbins" in tmpline:
                        stringtoreplace = "EstimateFromXbins=1"
                    elif "fitOrder=" in tmpline and "symmAfterFit" not in tmpline:
                        stringtoreplace = "fitOrder=" + tmpline.split("fitOrder=")[1][0]
                    elif "fitOrder=" in tmpline:
                        stringtoreplace = ["fitOrder=" + tmpline.split("fitOrder=")[1][0], "symmAfterFit=1"]
                    

                    for iB in range(len(vl.varList[thev]["bins_particle"]) - 1):
                        twotimestmpline = tmpline.replace("$tWProc", "tw_partbin{b}".format(b = iB)) + "\n"
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
                    twotimestmpline = tmpline.replace("$tWProc", "tw_bkg") + "\n"
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
                        outuncs += tmpline.replace("$tWProc", "tw_partbin{b}".format(b = iB)).replace("syst_tw_pdf_2016", "syst_tw_pdf_2016_partbin{b}".format(b = iB)).replace("syst_tw_pdf_2017", "syst_tw_pdf_2017_partbin{b}".format(b = iB)).replace("syst_tw_pdf_2018", "syst_tw_pdf_2018_partbin{b}".format(b = iB)) + "\n"
                    outuncs += tmpline.replace("$tWProc", "tw_bkg").replace("syst_tw_pdf_2016", "syst_tw_pdf_2016_bkg").replace("syst_tw_pdf_2017", "syst_tw_pdf_2017_bkg").replace("syst_tw_pdf_2018", "syst_tw_pdf_2018_bkg") + "\n"
                    
                else:
                    tmpline = line.strip()
                    stringtoreplace = None
                    if "EstimateFromXbins" in tmpline:
                        stringtoreplace = "EstimateFromXbins=1"
                    elif "fitOrder=" in tmpline and "symmAfterFit" not in tmpline:
                        stringtoreplace = "fitOrder=" + tmpline.split("fitOrder=")[1][0]
                    elif "fitOrder=" in tmpline:
                        stringtoreplace = ["fitOrder=" + tmpline.split("fitOrder=")[1][0], "symmAfterFit=1"]

                    if "$tWProc" in tmpline and stringtoreplace and ("particle_smoothing" in vl.varList[thev]):
                        theunc = tmpline.split(":")[0].replace(" ", "")
                        if theunc in vl.varList[thev]["particle_smoothing"]:
                            for iB in range(len(vl.varList[thev]["bins_particle"]) - 1):
                                twotimestmpline = tmpline.replace("$tWProc", "tw_partbin{b}".format(b = iB)) + "\n"
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
                            
                            twotimestmpline = tmpline.replace("$tWProc", "tw_bkg") + "\n"
                            outuncs += twotimestmpline
                        else:
                            outuncs += tmpline + "\n"
                    else:
                        outuncs += tmpline + "\n"
                    
    foutuncs = open(thepath + "/" + "tw-uncs-adhoc.txt", "w")
    foutuncs.write(outuncs)
    foutuncs.close(); del foutuncs
    return


def CardsCommand(prod, year, var, isAsimov, nthreads, outpath, noUnc, useFibre, extra_):
    cutsfile_  = "tw-run3/differential/cuts-differential/cuts-detector-1j1t.txt"

    samplespaths_ = "-P " + friendspath + "/" + prod + ("/" + year) * (year != "run3")
    if useFibre: samplespaths_ = samplespaths_.replace("phedexrw", "phedex").replace("cienciasrw", "ciencias")

    #nth_       = "" if nthreads == 0 else ("--split-factor=-1 -j " + str(nthreads))
    nth_       = "" if nthreads == 0 else ("--split-factor=0 -j " + str(nthreads))
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
                               lumi         = (vl.LumiDict[year] if year != "run3" else str(vl.LumiDict["2022"]) + "," + str(vl.LumiDict["2022PostEE"])),
                               variable     = variable_,
                               bins         = bins_,
                               nth          = nth_,
                               year         = year if year != "run3" else "2022,2022PostEE",
                               asimovornot  = "--asimov s+b" if isAsimov else "",
                               mcafile      = mcafile_,
                               cutsfile     = cutsfile_,
                               uncs         = "--unc {adhocu} --amc".format(adhocu = outpath_ + "/tw-uncs-adhoc.txt") if (not noUnc) else "--amc",
                               #uncs         = "--unc {adhocu} --amc".format(adhocu = "tw-run3/uncs-tw.txt") if (not noUnc) else "--amc",
                               name         = name_,
                               func         = "tw-run3/functions_tw.cc",
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



#def createCardsForEachSys(tsk):
#    inpath, iY, iV = tsk
#    print("    - Creating individual cards for variable " + iV + " of year " + iY)
#
#    path = inpath + "/" + iY + "/" + iV + "/sigextr_fit_combine/bincards"
#
#    #if not os.path.isdir(path + "/cards_{s}".format(s = syst)):
#        #os.system("mkdir " + path + "/cards_{s}".format(s = syst))
#
#    validfirstentries = {"*", "#", "bin", "observation", "shapes", "rate", "process"}
#
#    #listofuncsthataffecttheproc = []
#    #hasthissystitsowncard = False
#    #for s,el in vl.ModifiedSystMap.iteritems():
#        #if isinstance(el, dict):
#            #if el["tw" if not vl.unifttbar else "ttbar"]:
#                #listofuncsthataffecttheproc.append(s)
#        #elif el:
#            #listofuncsthataffecttheproc.append(s)
#
#    #if syst.replace("Up", "").replace("Down", "") in listofuncsthataffecttheproc:
#        #hasthissystitsowncard = True
#
#    txtfiletoopen  = path + "/card_bin{b}.txt"
#    rootfiletoopen = path + "/card_bin{b}.root"
#
#    for iB in range(len(vl.varList[iV]["bins_detector"]) - 1):
#        readF = open(txtfiletoopen.format(b = iB), 'r')
#
#        outtext = ""; theprocs = []
#        for line in readF:
#            entries = line.split()
#            if entries[0] in validfirstentries or entries[0][0] == "#":
#                tmpline = deepcopy(line)
#
#                if   "process" in entries[0] and not len(theprocs):
#                    for iP in range(1, len(entries)):
#                        theprocs.append(entries[iP])
#                elif "shapes" in entries[0]:
#                    tmpline = tmpline.replace("card_bin{b}.root".format(b = iB),
#                                              "finalcard_bin{b}_year{y}.root".format(b = iB, y = iY))
#
#                outtext += tmpline
#            elif len(entries) > 1:
#                if entries[0] == "jes" or entries[0] == "tw_matching": continue
#
#                if   ("shape" in entries[1]):
#                    goodline = ["shape"]
#                    ##Check to see whether there are zero evs. in any of the variations in any bin
#                    tmpRfile = r.TFile(rootfiletoopen.format(b = iB), "READ")
#                    for iP in range(2, len(entries)):
#                        toforget = False
#                        if "1" in entries[iP]:
#                            #print "x_" + theprocs[iP - 2] + "_" + entries[0].replace(" ", "")
#                            histoUp = tmpRfile.Get("x_" + theprocs[iP - 2] + "_" + entries[0].replace(" ", "") + "Up")
#                            histoDn = tmpRfile.Get("x_" + theprocs[iP - 2] + "_" + entries[0].replace(" ", "") + "Down")
#                            for jB in range(1, histoUp.GetNbinsX() + 1):
#                                if all([histoUp.GetBinContent(jB) < 1e-5, histoDn.GetBinContent(jB) < 1e-5]):
#                                    toforget = True
#                                    break
#                                if any([histoUp.GetBinContent(jB) < 1e-5, histoDn.GetBinContent(jB) < 1e-5]):
#                                    towarn = True
#                                    break
#
#                        if toforget:
#                            goodline.append("-")
#                        else:
#                            goodline.append(entries[iP])
#
#                    tmpRfile.Close(); del tmpRfile
#                    outtext += " ".join([entries[0]] + goodline) + "\n"
#                elif ("lnN"   in entries[1] and "norm" in entries[0]):
#                    outtext += line
#                else:
#                    goodline = ["shape"]
#                    for i in range(2, len(entries)):
#                        goodline.append("-" if "-" in entries[i] else "1")
#                    #print " ".join([entries[0]] + goodline)
#                    outtext += " ".join([entries[0]] + goodline) + "\n"
#
#        readF.close()
#        #print theprocs
#        outF  = open(path + "/finalcard_bin{b}_year{y}.txt".format(b = iB, y = iY), 'w')
#        outF.write(outtext)
#        outF.close()
#
#        print("    - Creating new rootfile card for bin", iB)
#        therootfile = r.TFile.Open(rootfiletoopen.format(b = iB), "READ")
#        tmpdictofthings = {}
#        copyname = ""; theproc = ""
#        for key in therootfile.GetListOfKeys():
#            copyname = key.GetName()
#            isaproc  = False
#            if len(copyname.split("_")) >= 3:
#                theproc  = deepcopy(copyname.split("_")[1] if ("bkg" not in copyname.split("_")[2] and "partbin" not in copyname.split("_")[2]) else "_".join(copyname.split("_")[1:3]))
#            else:
#                theproc  = deepcopy(copyname.split("_")[1])
#            #print theproc
#
#            if (copyname.split("_")[1] in vl.ProcessesNames and
#                len(copyname.split("_")) <= 2) or copyname == "x_data_obs":
#                isaproc = True
#            #else:
#                #theunc  = "_".join(copyname.split("_")[2:])
#                #pureunc = theunc.replace("Up", "").replace("Down", "")
#
#            if any([therootfile.Get(copyname).GetBinContent(jB) < 1e-5 for jB in range(1, therootfile.Get(copyname).GetNbinsX() + 1)]):
#                if isaproc:
#                    raise RuntimeError("FATAL: less than 1e-5 found in process " + theproc)
#                othervar = therootfile.Get(key.GetName().replace("Up" if "Up" in copyname else "Down", "Down" if "Up" in copyname else "Up"))
#                if any([(therootfile.Get(copyname).GetBinContent(jB) <= 0 and othervar.GetBinContent(jB) <= 0) for jB in range(1, therootfile.Get(copyname).GetNbinsX() + 1)]):
#                    continue
#                newhisto = deepcopy(therootfile.Get(copyname).Clone(copyname))
#                nomproc  = therootfile.Get("x_" + theproc)
#                for jB in range(1, newhisto.GetNbinsX() + 1):
#                    if newhisto.GetBinContent(jB) < 1e-5:
#                        if othervar.GetBinContent(jB) <= 1e-5:
#                            newhisto.SetBinContent(jB, 1e-5)
#                        else:
#                            thevar = othervar.GetBinContent(jB) - nomproc.GetBinContent(jB)
#                            if abs(thevar) >= nomproc.GetBinContent(jB):
#                                newhisto.SetBinContent(jB, 1e-5)
#                            else:
#                                newhisto.SetBinContent(jB, max(nomproc.GetBinContent(jB) - thevar, 1e-5))
#
#                tmpdictofthings[copyname] = deepcopy(newhisto)
#            else:
#                tmpdictofthings[copyname] = deepcopy(therootfile.Get(key.GetName()).Clone(copyname))
#
#        therootfile.Close(); del therootfile
#
#        theoutrootfile = r.TFile.Open(path + "/finalcard_bin{b}_year{y}.root".format(b = iB, y = iY), "RECREATE")
#        for key in tmpdictofthings:
#            tmpdictofthings[key].Write()
#        theoutrootfile.Close(); del theoutrootfile
#    towarndict = {}
#    return towarndict



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
    theyears = ["2022", "2022PostEE", "run3"]

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
#            listofuncsthataffecttheproc = []
#            for key,el in vl.ModifiedSystMap.iteritems():
#                if yr != "run3":
#                    skip = False
#                    for y in ["2022", "2022PostEE"]:
#                        if  y in key and yr != y:
#                            skip = True
#
#                    for y in [["22", "22PostEE"]]:
#                        if "".join(y) in key and not any([sub in yr for sub in y]):
#                            skip = True
#
#                    if skip:
#                        continue
#
#                listofuncsthataffecttheproc.append(key)
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

    #else: # preparing final cards
        #tasks = []
        #theyears = []
        #presentyears = next(os.walk(outpath))[1]

        #if "2016" in presentyears:
            #theyears.append("2016")
        #if "2017" in presentyears:
            #theyears.append("2017")
        #if "2018" in presentyears:
            #theyears.append("2018")
        #if "run2" in presentyears:
            #theyears.append("run2")


        #if year.lower() != "all" and year in presentyears:
            #theyears = [ year ]
        #elif year.lower() != "all":
            #raise RuntimeError("FATAL: the year requested is not in the provided input folder.")

        #for iY in theyears:
            #thevars = next(os.walk(outpath + "/" + iY))[1]

            #if variable.lower() != "all" and variable in thevars:
                #thevars = [ variable ]
            #elif variable.lower() != "all":
                #raise RuntimeError("FATAL: the variable requested is not in the provided input folder.")

            #for iV in thevars:
                #if "plots" in iV or "table" in iV: continue
                #if not os.path.isdir(outpath + "/" + iY + "/" + iV + "/sigextr_fit_combine"): continue


                #tasks.append( (outpath, iY, iV) )
        ##print tasks
        ##for el in tasks: print el
        ##sys.exit()
        #print "> Executing..."
        #if nthreads > 1:
            #pool = Pool(nthreads)
            #pool.map(createCardsForEachSys, tasks)
            #pool.close()
            #pool.join()
        #else:
            #for tsk in tasks:
                #createCardsForEachSys(tsk)
        #print "> Done!"
