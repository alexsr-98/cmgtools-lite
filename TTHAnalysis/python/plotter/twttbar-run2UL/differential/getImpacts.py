import ROOT  as r
import os, sys, argparse, math, json
import warnings as wr
from array import array
from copy  import deepcopy
from multiprocessing import Pool
import CombineHarvester.CombineTools.plotting as plot
import CombineHarvester.CombineTools.combine.rounding as rounding

sys.path.append('{cmsswpath}/src/CMGTools/TTHAnalysis/python/plotter/twttbar-run2UL/differential/'.format(cmsswpath = os.environ['CMSSW_BASE']))
import varList as vl

r.PyConfig.IgnoreCommandLineOptions = True
r.TH1.AddDirectory(0)

#comm1 = "combineTool.py -M Impacts -d {incard} --doInitialFit --robustFit 1 {ncores} {asimov} {extra} -m 1 -n {prefix} --out {outdir} --redefineSignalPOIs {pois} --floatOtherPOIs 1 --robustHesse 1 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000"
#comm2 = "combineTool.py -M Impacts -d {incard} --robustFit 1 --doFits {ncores} {asimov} {extra} -m 1 -n {prefix} --redefineSignalPOIs {pois} --robustHesse 1 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000"
#comm3 = "combineTool.py -M Impacts -d {incard} -o impacts{prefix}.json {ncores} {asimov} {extra} -m 1 -n {prefix} --redefineSignalPOIs {pois} --robustHesse 1 --robustFit 1 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000"

comm1 = "combineTool.py -M Impacts -d {incard} --doInitialFit --robustFit 1 {ncores} {asimov} {extra} -m 1 -n {prefix} --out {outdir} --redefineSignalPOIs {pois} --floatOtherPOIs 1 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000 --cminDefaultMinimizerStrategy 0 --robustHesse 1"
comm2 = "combineTool.py -M Impacts -d {incard} --robustFit 1 --doFits {ncores} {asimov} {extra} -m 1 -n {prefix} --redefineSignalPOIs {pois} --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000 --cminDefaultMinimizerStrategy 0 --robustHesse 1"
comm3 = "combineTool.py -M Impacts -d {incard} -o impacts{prefix}.json {ncores} {asimov} {extra} -m 1 -n {prefix} --redefineSignalPOIs {pois} --robustFit 1  --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000 --cminDefaultMinimizerStrategy 0 --robustHesse 1"

def GetRounded(nom, e_hi, e_lo):
    if e_hi < 0.0:
        e_hi = 0.0
    if e_lo < 0.0:
        e_lo = 0.0
    rounded = rounding.PDGRoundAsym(nom, e_hi if e_hi != 0.0 else 1.0, e_lo if e_lo != 0.0 else 1.0)
    s_nom = rounding.downgradePrec(rounded[0],rounded[2])
    s_hi = rounding.downgradePrec(rounded[1][0][0],rounded[2]) if e_hi != 0.0 else '0'
    s_lo = rounding.downgradePrec(rounded[1][0][1],rounded[2]) if e_lo != 0.0 else '0'
    return (s_nom, s_hi, s_lo)


def plotImpacts(inputjson, outputname, outputpath, npois, doBlind):
    left_margin   = 0.2
    pullsprop     = 0.37
    label_size    = 0.021
    doEraseStatNuis = False
    canv_width    = 700
    #canv_width    = 900
    canv_height   = 600
    per_page      = 30
    max_pages     = None
    cms_label     = "Internal"
    colour_groups = None
    units         = None
    mainPOI       = None
    xaxistitleoffset = 1
    xaxistitlesize = 0.0265

    externalPullDef = False

    # Dictionary to translate parameter names
    #vl.SysNameTranslator

    # Load the json output of combineTool.py -M Impacts
    data = {}
    with open(inputjson) as jsonfile:
        data = json.load(jsonfile)

    # Set the global plotting style
    plot.ModTDRStyle(l      = left_margin, 
                     b      = 0.10, 
                     width  = canv_width, 
                     height = canv_height)

    # We will assume the first POI is the one to plot
    POIs = [ele['name'] for ele in data['POIs']]
    #POI  = POIs[0]
    #if mainPOI:
        #POI = mainPOI

    poiDict = {}
    for ele in data['POIs']:
        poiDict[ele['name'].encode("utf-8").decode("utf-8")] = ele
    
    #for el in data["params"]:
        #if el["name"] == "dy_norm": print el
    #sys.exit()

    # Sort parameters by the average of the absolute impact on all POIs
    data['params'].sort(key = lambda x: vl.mean([abs(x['impact_' + iP]) for iP in POIs]), reverse = True)
    
    # Delete stat. nuis.?
    if doEraseStatNuis:
        for p in data['params']: 
            print(p)
            if p['name'].startswith('prop_binch'): 
                print('deleting', p)
                data['params'].remove(p)
    
    # Set the number of parameters per page (show) and the number of pages (n)
    show = per_page
    n    = int(math.ceil(float(len(data['params'])) / float(show)))
    if max_pages is not None and max_pages > 0:
        n = max_pages

    colour_hists       = {}
    colour_group_hists = {}

    if colour_groups is not None:
        colour_groups = {
            x.split('=')[0]: int(x.split('=')[1]) for x in colour_groups.split(',')
        }

    seen_types = set()

    for name, col in vl.nuisanceColours.items():
        colour_hists[name] = r.TH1F()
        plot.Set(colour_hists[name], FillColor = col, Title = name)

    if colour_groups is not None:
        for name, col in colour_groups.items():
            colour_group_hists[name] = r.TH1F()
            plot.Set(colour_group_hists[name], FillColor = col, Title = name)

    for page in range(n):
        canv     = r.TCanvas(outputname, outputname)
        n_params = len(data['params'][show * page:show * (page + 1)])
        pdata    = data['params'][show * page:show * (page + 1)]
        print('\t- Doing page %i, have %i parameters' % (page, n_params))

        boxes = []
        for i in range(n_params):
            y1 = r.gStyle.GetPadBottomMargin()
            y2 = 1. - r.gStyle.GetPadTopMargin()
            h  = (y2 - y1) / float(n_params)
            y1 = y1 + float(i) * h
            y2 = y1 + h
            box = r.TPaveText(0, y1, 1, y2, 'NDC')
            plot.Set(box, TextSize=0.02, BorderSize=0, FillColor=0, TextAlign=12, Margin=0.005)
            if i % 2 == 0:
                box.SetFillColor(18)
#            if (n_params - i + page * show - 3) <= 0:
            if (n_params - i + page * show) <= 0:
                box.AddText('-')
            else:
#                box.AddText('%i' % (n_params - i + page * show-3))
                box.AddText('%i' % (n_params - i + page * show))
            box.Draw()
            boxes.append(box)

        # Crate and style the pads
        pads = plot.MultiRatioSplitColumns([pullsprop] + [(1. - pullsprop - 0.01)/npois] * npois, [0.] * (npois + 1), [0.] * (npois + 1))
        #pads = plot.MultiRatioSplitColumns([(1. - pullsprop - 0.01)/npois] * npois * 2, [0.] * npois * 2, [0.] * npois * 2)
        for iP in range(npois):
            pads[iP + 1].SetGrid(1, 0)
        pads[0].SetGrid(1, 0)
        pads[0].SetTickx(1)
        pads[1].SetTickx(1)

        h_pulls = r.TH2F("pulls", "pulls", 6, -2.9, 2.9, n_params, 0, n_params)
        g_pulls = r.TGraphAsymmErrors(n_params)
        g_impacts_hi_list = []; g_impacts_lo_list = []
        max_impact_list = [0.] * npois
        for iP in range(npois):
            g_impacts_hi_list.append(r.TGraphAsymmErrors(n_params))
            g_impacts_lo_list.append(r.TGraphAsymmErrors(n_params))
            
        g_check = r.TGraphAsymmErrors()
        g_check_i = 0
        
        text_entries = []
        redo_boxes   = []
        print("\t- Beginning parameter analysis")
        #### Pulls
        for p in range(n_params):
            i   = n_params - (p + 1)
            pre = pdata[p]['prefit']
            fit = pdata[p]['fit']
            tp  = pdata[p]['type']
            seen_types.add(tp)
            if pdata[p]['type'] != 'Unconstrained':
                pre_err_hi = (pre[2] - pre[1])
                pre_err_lo = (pre[1] - pre[0])

                pull = fit[1] - pre[1]
                pull = (pull/pre_err_hi) if pull >= 0 else (pull/pre_err_lo)
                pull_hi = fit[2] - pre[1]
                pull_hi = (pull_hi/pre_err_hi) if pull_hi >= 0 else (pull_hi/pre_err_lo)
                pull_hi = pull_hi - pull
                pull_lo = fit[0] - pre[1]
                pull_lo = (pull_lo/pre_err_hi) if pull_lo >= 0 else (pull_lo/pre_err_lo)
                pull_lo =  pull - pull_lo

                g_pulls.SetPoint(i, pull, float(i) + 0.5)
                g_pulls.SetPointError(
                    i, pull_lo, pull_hi, 0., 0.)
            else:
                # Hide this point
                g_pulls.SetPoint(i, 0., 9999.)
                y1 = r.gStyle.GetPadBottomMargin()
                y2 = 1. - r.gStyle.GetPadTopMargin()
                x1 = r.gStyle.GetPadLeftMargin()
                h = (y2 - y1) / float(n_params)
                y1 = y1 + ((float(i)+0.5) * h)
                x1 = x1 + (1 - pads[0].GetRightMargin() -x1)/2.
                s_nom, s_hi, s_lo = GetRounded(fit[1], fit[2] - fit[1], fit[1] - fit[0])
                text_entries.append((x1, y1, '%s^{#plus%s}_{#minus%s}' % (s_nom, s_hi, s_lo)))
                redo_boxes.append(i)
             
            for iP in range(npois):
                g_impacts_hi_list[iP].SetPoint(i, 0, float(i) + 0.5)
                g_impacts_lo_list[iP].SetPoint(i, 0, float(i) + 0.5)
            
            
            imp_list = []
            for nam in poiDict:
                imp_list.append(pdata[p][nam])
                
            for iP in range(npois):
                g_impacts_hi_list[iP].SetPointError(i, 0, imp_list[iP][2] - imp_list[iP][1], 0.5, 0.5)
                g_impacts_lo_list[iP].SetPointError(i, imp_list[iP][1] - imp_list[iP][0], 0, 0.5, 0.5)
                max_impact_list[iP] = max(
                    max_impact_list[iP], abs(imp_list[iP][1] - imp_list[iP][0]), abs(imp_list[iP][2] - imp_list[iP][1]))
                
            col = vl.nuisanceColours.get(tp, 2)
            if colour_groups is not None and len(pdata[p]['groups']) >= 1:
                for p_grp in pdata[p]['groups']:
                    if p_grp in colour_groups:
                        col = colour_groups[p_grp]
                        break
            h_pulls.GetYaxis().SetBinLabel(
                #i + 1, ('#color[%i]{%s}'% (col, vl.SysNameTranslator[pdata[p]['name'].encode('utf-8')])))
                i + 1, ('#color[%i]{%s}'% (col, pdata[p]['name'].encode('utf-8'))))

        # Style and draw the pulls histo
        plot.Set(h_pulls.GetXaxis(), TitleSize = xaxistitlesize, LabelSize = 0.015, TitleOffset=xaxistitleoffset-0.09, Title='(#hat{#theta}-#theta_{0})/#Delta#theta')
        plot.Set(h_pulls.GetYaxis(), LabelSize = label_size, TickLength=0.0)
        h_pulls.GetYaxis().LabelsOption('v')
        h_pulls.Draw()

        for i in redo_boxes:
            newbox = boxes[i].Clone()
            newbox.Clear()
            newbox.SetY1(newbox.GetY1()+0.005)
            newbox.SetY2(newbox.GetY2()-0.005)
            newbox.SetX1(r.gStyle.GetPadLeftMargin()+0.001)
            newbox.SetX2(0.7-0.001)
            newbox.Draw()
            boxes.append(newbox)
        latex = r.TLatex()
        latex.SetNDC()
        latex.SetTextFont(42)
        latex.SetTextSize(0.02)
        latex.SetTextAlign(22)
        for entry in text_entries:
            latex.DrawLatex(*entry)

        # Go to the other pad and draw the impacts histo
        h_impacts_list = []
        for iP in range(npois):
            pads[iP + 1].cd()
            if max_impact_list[iP] == 0.: max_impact_list[iP] = 1E-6  # otherwise the plotting gets screwed up
            h_impacts_list.append(r.TH2F("impacts_" + str(iP), 
                                         "impacts_" + str(iP), 
                                         6, 
                                         -max_impact_list[iP] * 1.1, 
                                         max_impact_list[iP] * 1.1, 
                                         n_params, 
                                         0, 
                                         n_params))
            plot.Set(h_impacts_list[-1].GetXaxis(), 
                     LabelSize  = 0.015, 
                     TitleSize  = xaxistitlesize, 
                     TitleOffset=xaxistitleoffset,
                     Ndivisions = 505, 
                     Title      = '#Delta#hat{#mu_{%s}}'%(iP))
            plot.Set(h_impacts_list[-1].GetYaxis(), 
                     LabelSize  = 0, 
                     TickLength = 0.0)
            h_impacts_list[-1].Draw()
        print(len(pads))

        # Back to the first pad to draw the pulls graph
        pads[0].cd()
        plot.Set(g_pulls, 
                 MarkerSize = 0.8, 
                 LineWidth  = 2)
        g_pulls.Draw('PSAME')

        # And back to the second pad to draw the impacts graphs
        alpha  = 0.7
        method = 'default'
        if 'method' in data and data['method'] in vl.lo_colour:
            method = data['method']
        
        for iP in range(npois):
            pads[iP + 1].cd()
            g_impacts_hi_list[iP].SetFillColor(plot.CreateTransparentColor(vl.hi_colour[method], alpha))
            g_impacts_hi_list[iP].Draw('2SAME')
            g_impacts_lo_list[iP].SetFillColor(plot.CreateTransparentColor(vl.lo_colour[method], alpha))
            g_impacts_lo_list[iP].Draw('2SAME')
            pads[iP + 1].RedrawAxis()
            

        legend = r.TLegend(0.02, 0.02, 0.30, 0.04, '', 'NBNDC')
        legend.SetNColumns(3)
        legend.AddEntry(g_pulls, 'Pull', 'LP')
        legend.AddEntry(g_impacts_hi_list[0], '+1#sigma Impact', 'F')
        legend.AddEntry(g_impacts_lo_list[0], '-1#sigma Impact', 'F')
        legend.Draw()

        leg_width = pads[0].GetLeftMargin() - 0.01
        if colour_groups is not None:
            legend2 = r.TLegend(0.01, 0.94, leg_width, 0.99, '', 'NBNDC')
            legend2.SetNColumns(2)
            for name, h in colour_group_hists.items():
#                legend2.AddEntry(h, vl.SysNameTranslator[name], 'F')
                legend2.AddEntry(h, name, 'F')
            legend2.Draw()
        elif len(seen_types) > 1:
            legend2 = r.TLegend(0.01, 0.94, leg_width, 0.99, '', 'NBNDC')
            legend2.SetNColumns(2)
            for name, h in colour_hists.items():
                if name == 'Unrecognised': continue
                legend2.AddEntry(h, name, 'F')
            legend2.Draw()

        plot.DrawCMSLogo(pads[0], 'CMS', cms_label, 0, 0.27, 0.0, 0.00, cmsTextSize = 0.5)
        s_nom_list = []; s_hi_list = []; s_lo_list = []
        
        for iP in range(npois):
            s_nom, s_hi, s_lo = GetRounded(poiDict["r_bb4l_" + str(iP)]["fit"][1], 
                                           poiDict["r_bb4l_" + str(iP)]["fit"][2] - poiDict["r_bb4l_" + str(iP)]["fit"][1], 
                                           poiDict["r_bb4l_" + str(iP)]["fit"][1] - poiDict["r_bb4l_" + str(iP)]["fit"][0])
            s_nom_list.append(s_nom); s_hi_list.append(s_hi); s_lo_list.append(s_lo);
            
        if not doBlind:
            for iP in range(npois):
                plot.DrawTitle(pads[iP + 1], '#hat{#mu_{%s}} = %s^{#plus%s}_{#minus%s}%s' % (
                    iP, s_nom_list[iP], s_hi_list[iP], s_lo_list[iP],
                    '' if units is None else ' ' + units), 3, 0.27, 0.3)

        extra = ''
        if page == 0:
            extra = '('
        if page == n - 1:
            extra = ')'
        canv.Print('{path}/{nam}.pdf{ex}'.format(path = outputpath, nam = outputname, ex = extra))
    return


def makeImpacts(task):
    inpath, year, varName, ncores, pretend, verbose, extra, doobs, doblind = task

    print('\n> Creating impacts for variable', varName, 'from year', year, '\n')
    bins_detector = vl.varList[varName]['bins_detector']
    ndetectorbins   = len(bins_detector) - 1
    bins_particle = vl.varList[varName]['bins_particle']
    nparticlebins   = len(bins_particle) - 1

    thepois = ",".join( ["r_bb4l_" + str(i) for i in range(nparticlebins)] )

    impactsoutpath = inpath + "/" + year + "/" + varName + "/sigextr_fit_combine/{o}Impacts/".format(o = "Obs" if doObs else "")

    asimov_ = "--setParameters "
    for idx in range(nparticlebins - 1):
        asimov_ += "r_bb4l_{iBp}=1,".format(iBp = idx)
    asimov_ += "r_bb4l_{iBp}=1".format(iBp = nparticlebins - 1)
    if not doobs:
        asimov_ +=  " -t -1 "


    firstcomm = comm1.format(y      = year,
                             ncores = ("--parallel " + str(ncores)) if ncores else "",
                             asimov = asimov_,
                             incard = "../fit_output.root",
                             outdir = "./",
                             var    = varName,
                             extra  = extra,
                             prefix = year + "_" + varName,
                             pois   = thepois,
    )

    if verbose:
        print("First command:", firstcomm, "\n")

    if not pretend:
        outstat = os.system("cd " + impactsoutpath + "; " + firstcomm + "; cd -")
        if outstat:
            raise RuntimeError("FATAL: first command failed to execute for variable {v} of year {y}.".format(v = varName,
                                                                                                             y = year))

    secondcomm = comm2.format(y     = year,
                             ncores = ("--parallel " + str(ncores)) if ncores else "",
                             asimov = asimov_,
                             incard = "../fit_output.root",
                             outdir = "./",
                             var    = varName,
                             extra  = extra,
                             prefix = year + "_" + varName,
                             pois   = thepois,
    )

    if verbose:
        print("Second command:", secondcomm, "\n")

    if not pretend:
        outstat = os.system("cd " + impactsoutpath + "; " + secondcomm + "; cd -")
        if outstat:
            raise RuntimeError("FATAL: second command failed to execute for variable {v} of year {y}.".format(v = varName,
                                                                                                              y = year))

    thirdcomm = comm3.format(y      = year,
                             ncores = ("--parallel " + str(ncores)) if ncores else "",
                             asimov = asimov_,
                             incard = "../fit_output.root",
                             outdir = "./",
                             var    = varName,
                             extra  = extra,
                             prefix = year + "_" + varName,
                             pois   = thepois,
    )

    if verbose:
        print("Third command:", thirdcomm, "\n")

    if not pretend:
        outstat = os.system("cd " + impactsoutpath + "; " + thirdcomm + "; cd -")
        if outstat:
            raise RuntimeError("FATAL: third command failed to execute for variable {v} of year {y}.".format(v = varName,
                                                                                                             y = year))

    
    plotImpacts(impactsoutpath + "/impacts{y}_{v}.json".format(y = year, v = varName), "impacts_{y}_{v}".format(y = year, v = varName), impactsoutpath, nparticlebins, doblind)
    
    print('\n> Variable', varName, "' impacts produced.\n")
    return



if __name__ == '__main__':
    vl.SetUpWarnings()
    r.gROOT.SetBatch(True)
    print("===== Fitting procedure with some uncertainty profiling\n")
    parser = argparse.ArgumentParser(usage = "python nanoAOD_checker.py [options]", description = "Checker tool for the outputs of nanoAOD production (NOT postprocessing)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inpath',     '-i', metavar = 'inpath',     dest = "inpath",   required = False, default = "./temp/differential/")
    parser.add_argument('--year',       '-y', metavar = 'year',       dest = "year",     required = False, default = "all")
    parser.add_argument('--variable',   '-v', metavar = 'variable',   dest = "variable", required = False, default = "all")
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
    variable = args.variable
    verbose  = args.verbose
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

    for iY in theyears:
        thevars = next(os.walk(inpath + "/" + iY))[1]

        if variable.lower() != "all" and variable in thevars and "," not in variable:
            thevars = [ variable ]
        elif "," in variable:
            if len(set(variable.split(",")) - set(thevars)):
                raise RuntimeError("FATAL: one of the variables requested is not in the provided input folder.")
            else:
                thevars = variable.split(",")
        elif variable.lower() != "all":
            raise RuntimeError("FATAL: the variable requested is not in the provided input folder.")

        for iV in thevars:
            if "plots" in iV or "Fiducial" in iV or "table" in iV or "response" in iV: continue
            if not os.path.isdir(inpath + "/" + iY + "/" + iV + "/sigextr_fit_combine/{o}Impacts".format(o = "Obs" if doObs else "")):
                os.system("mkdir -p " + inpath + "/" + iY + "/" + iV + "/sigextr_fit_combine/{o}Impacts".format(o = "Obs" if doObs else ""))
            tasks.append( (inpath, iY, iV, nthreads, pretend, verbose, extra, doObs, doBlind) )

    #print tasks
    for task in tasks:
        print("\nProcessing " + str(task) + "\n")
        makeImpacts(task)
