import ROOT as r
#import tabulate as tb
import prettytable as pt
import warnings as wr
from copy import deepcopy
import sys, os, argparse
from collections import OrderedDict
from decimal import Decimal

sys.path.append('{cmsswpath}/src/CMGTools/TTHAnalysis/python/plotter/twttbar-run2UL/differential/'.format(cmsswpath = os.environ['CMSSW_BASE']))
import varList as vl
from tableFormats import latexTable
r.gROOT.SetBatch(True)


def saveLaTeXfromhisto(histo, varname, path = "./", outname = "", ty = "detector", errhisto = None):
    table   = [[], []]; headers = []
    usehisto = r.TH1D()

    if errhisto != None:
        usehisto = deepcopy(histo.Clone("usehisto"))
        for bn in range(1, usehisto.GetNbinsX() + 1): usehisto.SetBinError(bn, errhisto.GetBinError(bn))
    else: usehisto = deepcopy(histo.Clone("usehisto"))


    nbins   = usehisto.GetNbinsX()
    limits  = [int(round(usehisto.GetXaxis().GetBinLowEdge(bn), 0)) for bn in range(1, nbins + 1)]

    #### Construction of headers
    headers.append(vl.varList[varname]["printname"])
    if "DPhi" not in varname:
        for b in range(1, nbins + 1):
            if b == nbins:
                headers.append("[{inf}, $\\infty$)".format(inf = limits[b - 1]))
            else:
                headers.append("[{inf}, {sup})".format(inf = limits[b - 1], sup = limits[b]))
    else:
        for b in range(1, nbins + 1):
            if b == nbins:
                headers.append("[{inf}, 1]".format(inf = round(usehisto.GetXaxis().GetBinLowEdge(b), 2)))
            else:
                headers.append("[{inf}, {sup})".format(inf = round(usehisto.GetXaxis().GetBinLowEdge(b), 2), sup = round(usehisto.GetXaxis().GetBinLowEdge(b + 1), 2)))

    #### Construction of the table
    table[0].append(vl.TableDict[ty]["prefix"] + vl.varList[varname]["mathprintname"] + "\\right)$ " + vl.TableDict[ty]["dimensions"])
    table[1].append("Uncertainty " + vl.TableDict[ty]["dimensions"])
    for bn in range(1, nbins + 1):
        tmpval, tmpinc = getvalueandincstrings(usehisto, bn)
        table[0].append(tmpval); table[1].append(tmpinc)

    # To assure that we don't remove trailing zeros
    #frm = tuple([""] + [".{l}f".format(l = len(table[0][ind].split(".")[1] )) for ind in range(1, len(table[0])) ])
    frm = tuple([""] + [
    ".{l}f".format(l=len(table[0][ind].split(".")[1])) 
    if "." in table[0][ind] else ".0f" 
    for ind in range(1, len(table[0]))
    ])

    #finaltab = tb.tabulate(table, headers = headers, floatfmt = frm, numalign = "decimal", stralign = "center", colalign=("left",))
    #print(finaltab)

    #print "\n\n"

    #finallatextab = tb.tabulate(table, headers = headers, floatfmt = frm, tablefmt = "latex_raw")
    #print(finallatextab)
    print(headers)
    finallatextab = pt.PrettyTable(headers)

    finallatextab.add_row(table[0])
    finallatextab.add_row(table[1])

    #finallatextab.float_format = '0.2'
    finallatextab = latexTable(finallatextab)
    thetab = finallatextab.get_string()

    substr = "".join(["c"] + ["c" for x in range(nbins)])
    thetab = thetab.replace(substr, "".join(["l|"] + ["c|" for x in range(nbins - 1)] + ["c"] ))


    #print thetab
    #sys.exit()

    if not os.path.isdir(path): os.system("mkdir -p " + path)

    name = outname if outname != "" else (varname + "_" + ty)

    outfile = open(path + "/" + name + ".tex", "w")
    outfile.write(thetab)
    outfile.close(); del outfile

    return


def getvalueandincstrings(histo, b):
    valstring = ""; incstring = ""; # NOTE: Python by default uses up to 12 characters when
                                    #       transforming a float number into a string.
    valstring = str(Decimal(histo.GetBinContent(b)))
    incstring = str(Decimal(histo.GetBinError(b)))
    ndeczeros = 0
    numofdecs = 0

    if "e" in (valstring + incstring):
        raise RuntimeError("Scientific notation detected: case not supported. Histo name: " + histo.GetName() + ", value: " + valstring + ", unc.: " + incstring)

    if "." in incstring:
        incparts = incstring.split(".")

        if incparts[0].replace("0", "") == "": # There is no non-decimal part
            # Count figures that are zero
            if incparts[1][0] == "0":
                selfig = "0"
                while selfig == "0":
                    selfig = incparts[1][ndeczeros + 1]
                    ndeczeros += 1
                    if ndeczeros == len(incparts[1]) and selfig == "0":
                        raise RuntimeError("The uncertainty is completely zero!")

            numofdecs = ndeczeros + 1 # number of total decimals that we will take

            # Check if the first significant figure is a one or not.
            if incparts[1][ndeczeros] == "1":
                numofdecs += 1
        else:
            if   len(incparts[0]) == 1 and incparts[0][0] != 1: numofdecs = 0
            elif len(incparts[0]) == 1:                         numofdecs = 1
            else:
                wr.warn("WARNING: the result is a real number with only integer part and more than one digit. This should be saved with scientific notation, but it is not implemented. The full integer will be returned.")
                numofdecs = 0

        finalvalstring = str(round(histo.GetBinContent(b), numofdecs))
        finalincstring = str(round(histo.GetBinError(b), numofdecs))

        if   len(finalincstring) < len(finalvalstring):
            finalincstring += "0"
        elif len(finalincstring) > len(finalvalstring):
            finalvalstring += "0"
    else:
        finalvalstring = str(histo.GetBinContent(b))
        finalincstring = str(histo.GetBinError(b))

    return finalvalstring, finalincstring


#def getyieldsLaTeXtable(xuanpath, outname = "yields", path = "./results/tables"):
    #xuanchivo = open(xuanpath, "r")

    #recordthings = False
    #procs = []; table = [[]];

    #procs.append("Source")
    #table[0].append("Amount")

    #for line in xuanchivo.readlines():
        #if "====" in line: recordthings = True
        #if recordthings and "|" in line:
            #sublines = line.replace(" ", "").split("|")
            #procs.append(sublines[0].replace("Background", "") if not sublines[0] in vl.ProcessNameTranslator else vl.ProcessNameTranslator[sublines[0]])
            #table[0].append(round(float(sublines[1].split("+-")[0]), 0))
    #xuanchivo.close()

    #allbkg = 0
    #for ind in range(1, len(procs)):
        ##print(procs[ind].lower())
        #if ("signal"  not in procs[ind].lower() and
            #"total"   not in procs[ind].lower() and
            #"data"    not in procs[ind].lower()):
            #allbkg += table[0][ind]

    #for ind in range(len(procs)):
        #if "signal" in procs[ind].lower():
            #procs.insert(ind, "Total background")
            #table[0].insert(ind, round(float(allbkg), 0))
            #break

    #finallatextab = tb.tabulate(table, headers = procs, tablefmt = "latex_raw")

    #substr = "".join(["l"] + ["r" for x in range(len(table[0]) - 1)])
    #finallatextab = finallatextab.replace(substr, "".join(["l|"] + ["c|" for x in range(len(table[0]) - 2)] + ["c"] ))

    ##print finallatextab
    #if not os.path.isdir(path): os.system("mkdir -p " + path)

    #name = outname if outname != "" else (varname + "_" + ty)

    #outfile = open(path + "/" + outname + ".tex", "w")
    #outfile.write(finallatextab)
    #outfile.close(); del outfile
    #return


#def getcondnumLaTeXtable(path, variables, outname = "condnum", addasimov = False):
    #table   = []; headers = []

    #if len(variables) == 0: raise RuntimeError("No variables given.")

    ## First, obtain the info from the txts
    #listofsysts  = []
    #listoftuples = []
    #filledsysts  = False
    #for var in variables:
        #tmptupleofcontents = (vl.varList[var]["printname"], [])

        #thef = path + "/" + var + "/responseplots/condnums.txt"
        #if not os.path.isfile(thef):
            #raise RuntimeError("FATAL: file " + thef + " not found." )
        #tmpfile = open(thef, "r")

        #for line in tmpfile.readlines():
            #if ":" not in line: continue
            #sublines = line.strip().replace(" ", "").split(":")

            #if len(sublines) > 1 and not filledsysts:
                #thegoodname = sublines[0].replace("Up", "").replace("Down", "")

                #print thegoodname

                #listofsysts.append(vl.PrintSysNameTranslator[thegoodname] + " " + ("Up" if "Up" in sublines[0] else "Down"))

            #tmptupleofcontents[1].append(float(sublines[1]))

        #listoftuples.append(tmptupleofcontents)
        #filledsysts = True
        #tmpfile.close()

    #listoftuples.insert(0, ("Syst./Variable", listofsysts))
    #thedict = OrderedDict(listoftuples)

    ##sys.exit()

    #finallatextab = tb.tabulate(thedict, headers = "keys", floatfmt = ".2f", tablefmt = "latex_raw")

    #substr = "".join(["l"] + ["r" for x in range(len(variables))])
    #finallatextab = finallatextab.replace(substr, "".join(["l|"] + ["c|" for x in range(len(variables) - 1)] + ["c"] ))

    #print finallatextab

    #if not os.path.isdir(path + "/tables"):
        #os.system("mkdir -p " + path + "/tables")

    #outfile = open(path + "/tables/" + outname + ".tex", "w")
    #outfile.write(finallatextab)
    #outfile.close(); del outfile
    #return


def getcondnumLaTeXtable(path, variables, outname = "condnum", addasimov = False):
    table   = []; headers = []

    if len(variables) == 0: raise RuntimeError("No variables given.")

    # First, obtain the info from the txts
    listofsysts  = []
    listoftuples = []
    filledsysts  = False
    for var in variables:
        tmptupleofcontents = (vl.varList[var]["printname"], [])

        thef = path + "/" + var + "/responseplots/condnums.txt"
        if not os.path.isfile(thef):
            raise RuntimeError("FATAL: file " + thef + " not found." )
        tmpfile = open(thef, "r")

        for line in tmpfile.readlines():
            if ":" not in line: continue
            sublines = line.strip().replace(" ", "").split(":")

            if len(sublines) > 1 and not filledsysts:
                thegoodname = sublines[0].replace("Up", "").replace("Down", "")

                #print thegoodname
                if "tw_matching" in thegoodname or thegoodname == "jes": continue
                if not "lumi" in thegoodname:
                    listofsysts.append(vl.PrintSysNameTranslator[thegoodname] + " " + ("Up" if "Up" in sublines[0] else "Down" if not "nominal" in thegoodname else ""))
                else:
                    listofsysts.append(thegoodname.replace("_", "\\_") + " " + ("Up" if "Up" in sublines[0] else "Down" if "nominal" not in sublines[0].lower() else ""))

            tmptupleofcontents[1].append(float(sublines[1]))

        listoftuples.append(tmptupleofcontents)
        filledsysts = True
        tmpfile.close()

    listoftuples.insert(0, ("Syst./Variable", listofsysts))
    thedict = OrderedDict(listoftuples)

    #sys.exit()

    #finallatextab = tb.tabulate(thedict, headers = "keys", floatfmt = ".2f", tablefmt = "latex_raw")

    theheaders = [el[0] for el in listoftuples]

    finallatextab = pt.PrettyTable(theheaders)

    for iS in range(len(listofsysts)):
        finallatextab.add_row([el[1][iS] for el in listoftuples])


    finallatextab.float_format = '0.2'
    finallatextab = latexTable(finallatextab)
    thetab = finallatextab.get_string()

    substr = "".join(["c"] + ["c" for x in range(len(variables))])
    thetab = thetab.replace(substr, "".join(["l|"] + ["c|" for x in range(len(variables) - 1)] + ["c"] ))
    #print thetab
    #sys.exit()

    if not os.path.isdir(path + "/tables"):
        os.system("mkdir -p " + path + "/tables")

    outfile = open(path + "/tables/" + outname + ".tex", "w")
    outfile.write(thetab)
    outfile.close(); del outfile
    return


#def getgoftestsLaTeXtable(variables, path = "./results/tables", txttablespath = "./", outname = "", ty = "unfolding"):
    #table   = []; headers = []

    #if len(variables) == 0: raise RuntimeError("No variables given.")

    ## First, obtain the info from the txts
    #listofsamples = []
    #filledsamples = False
    #for var in variables:
        #tmplistofcontents = [ vl.varList[var]["printnamenodim"] ]
        #tmpfile = open(txttablespath + "/" + var + "_goftests_" + ty + ".txt", "r")
        #tmpsample = ""
        #for line in tmpfile.readlines(): # 1 p-value, 2 test statistic value
            #if "/" not in line: continue
            #sublines = line.replace(" ", "").split("/")
            #if "p-value" in sublines[1]:
                #tmplistofcontents.append( str( round(float(sublines[1].split(":")[1]), 3) ) )
                #if not filledsamples: listofsamples.append(vl.GOFTranslator[sublines[0]])
            #else:
                #tmplistofcontents[-1] += " | " + str( round( float(sublines[1].split(":")[1]), 3) )

        #table.append(tmplistofcontents)
        #filledsamples = True
        #tmpfile.close();


    #listofsamples.insert(0, "Variable/MC sample")
    ##sys.exit()

    #finallatextab = tb.tabulate(table, headers = listofsamples, tablefmt = "latex_raw")

    #substr = "".join(["l"] + ["l" for x in range(len(listofsamples) - 1)])
    #finallatextab = finallatextab.replace(substr, "".join(["l|"] + ["c|" for x in range(len(listofsamples) - 2)] + ["c"] ))


    #print finallatextab
    #if not os.path.isdir(path): os.system("mkdir -p " + path)

    #name = outname if outname != "" else "goftests"

    #outfile = open(path + "/" + name +"_" + ty + ".tex", "w")
    #outfile.write(finallatextab)
    #outfile.close(); del outfile
    #return


def getgoftestsLaTeXtable(variables, path = "./results/tables", txttablespath = "./", outname = "", ty = "particlefidbin"):
    table   = []; headers = []

    if len(variables) == 0: raise RuntimeError("No variables given.")

    # First, obtain the info from the txts
    listofsamples = []
    filledsamples = False
    for var in variables:
        tmplistofcontents = [ vl.varList[var]["printnamenodim"] ]
        tmpfile = open(txttablespath + "/" + var + "/goftestplots/" + var + "_goftests_" + ty + ".txt", "r")
        tmpsample = ""
        for line in tmpfile.readlines(): # 1 p-value, 2 test statistic value
            if "/" not in line: continue
            sublines = line.replace(" ", "").split("/")
            if "p-value" in sublines[1]:
                tmplistofcontents.append( str( round(float(sublines[1].split(":")[1]), 3) ) )
                if not filledsamples: listofsamples.append(vl.GOFTranslator[sublines[0]])
#            else:
#                tmplistofcontents[-1] += " | " + str( round( float(sublines[1].split(":")[1]), 3) )

        table.append(tmplistofcontents)
        filledsamples = True
        tmpfile.close();


    listofsamples.insert(0, "Variable/MC sample")
    #sys.exit()

    ### VEYO
    #finallatextab = tb.tabulate(table, headers = listofsamples, tablefmt = "latex_raw")

    #substr = "".join(["l"] + ["l" for x in range(len(listofsamples) - 1)])
    #finallatextab = finallatextab.replace(substr, "".join(["l|"] + ["c|" for x in range(len(listofsamples) - 2)] + ["c"] ))
    #print finallatextab

    ### NOVO
    finallatextab = pt.PrettyTable(listofsamples)

    for el in table:
        finallatextab.add_row(el)

    finallatextab.float_format = '0.2'
    finallatextab = latexTable(finallatextab)
    thetab = finallatextab.get_string()

    substr = "".join(["c"] + ["c" for x in range(len(variables))])
    thetab = thetab.replace(substr, "".join(["l|"] + ["c|" for x in range(len(variables) - 1)] + ["c|"] ))
    print(thetab)
    #sys.exit()


    if not os.path.isdir(path): os.system("mkdir -p " + path)

    name = outname if outname != "" else "goftests"

    outfile = open(path + "/" + name +"_" + ty + ".tex", "w")
    outfile.write(thetab)
    outfile.close(); del outfile
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage = "python nanoAOD_checker.py [options]", description = "Checker tool for the outputs of nanoAOD production (NOT postprocessing)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inpath',    '-i', metavar = 'inpath',     dest = "inpath",   required = True,  default = "./temp/differential")
    parser.add_argument('--type',      '-t', metavar = 'thetype',    dest = "thetype",  required = False, default = "condnumtable")
    parser.add_argument('--year',      '-y', metavar = 'year',       dest = "year",     required = False, default = "all")
    parser.add_argument('--variable',  '-v', metavar = 'variable',   dest = "variable", required = False, default = "all")


    args     = parser.parse_args()
    year     = args.year
    inpath   = args.inpath
    variable = args.variable
    thetype  = args.thetype

    #print("\n===== Unfolding procedures: Response matrices & ROOT files production =====")
    #print("> Setting binning, paths, and other details...")

    #print "\n> Drawing matrices and writing ROOT file (old one will be overwritten!)."


    #### First, find the tasks
    tasks = []
    if year == "all":
        if variable == "all":
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

            for iY in theyears:
                thevars = next(os.walk(inpath + "/" + iY))[1]
                actualvars = []
                for iV in thevars:
                    if any([el in iV for el in vl.vetolist]): continue
                    if iV == "Fiducial": continue
                    actualvars.append(iV)
                tasks.append( (inpath + "/" + iY, actualvars) )

    #tasks = [ (inpath, "2016", "Lep1Lep2Jet1MET_Mt") ]

    #print tasks

    #sys.exit()


    if thetype == "condnumtable":
        print("> Producing LaTeX tables for the condition numbers")
        thefunc = "getcondnumLaTeXtable"

        for tsk in tasks:
            eval(thefunc + "(tsk[0], tsk[1])")

print("> Done!")









