from copy import deepcopy
import ROOT as r
import warnings as wr
import os, sys, math, array, argparse
from multiprocessing import Pool
from array import array

sys.path.append('{cmsswpath}/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/differential/'.format(cmsswpath = os.environ['CMSSW_BASE']))
import beautifulUnfoldingPlots as bp
import errorPropagator as ep
import getLaTeXtable as tex
import varList as vl
import tdrstyle, CMS_lumi

vl.SetUpWarnings()
r.gROOT.SetBatch(True)
verbose = False
doArXiv     = True



def BibhuFunction(dv, df, dfinal, dcov, varName):
    for key in dv:
        #if key == "LeadingLepPt_asimov": print key
        
        for bin in range(1, dv[key].GetNbinsX() + 1):
            goodunc = r.TMath.Sqrt( dv[key].GetBinError(bin)**2 / df[key.replace(varName, "Fiducial")].GetBinContent(bin)**2 +
                                    dv[key].GetBinContent(bin)**2 * df[key.replace(varName, "Fiducial")].GetBinError(bin)**2 / df[key.replace(varName, "Fiducial")].GetBinContent(bin)**4 -
                                    2 * dv[key].GetBinContent(bin) / df[key.replace(varName, "Fiducial")].GetBinContent(bin)**3 * sum([dcov[key].GetBinContent(bin, j) for j in range(1, dv[key].GetNbinsX() + 1)]) )
            
            dfinal[key].SetBinError(bin, goodunc)
    return


def BibhuFunctionForCovs(cov, vals, var, outnam = "CovMat_fidbin"):
    partbins = array("d", vl.varList[var]["bins_particle"])
    outcov   = r.TH2D(outnam, "", len(vl.varList[var]["bins_particle"]) - 1, partbins,
                                  len(vl.varList[var]["bins_particle"]) - 1, partbins)


    fidval  = sum([vals.GetBinContent(i) for i in range(1,   vals.GetNbinsX() + 1)])
    fidunc2 = sum([cov.GetBinContent(i, j) for i in range(1, cov.GetNbinsX() + 1) for j in range(1, cov.GetNbinsX() + 1) ])

    for iB in range(1, outcov.GetNbinsX() + 1):
        for jB in range(iB, outcov.GetNbinsY() + 1):
            outcov.SetBinContent(iB, jB,
                                (cov.GetBinContent(iB, jB)/fidval**2 +
                                 vals.GetBinContent(iB) * vals.GetBinContent(jB) * fidunc2 / fidval**4 -
                                 vals.GetBinContent(iB) / (fidval**3) * sum([cov.GetBinContent(jB, k) for k in range(1, vals.GetNbinsX() + 1)]) -
                                 vals.GetBinContent(jB) / (fidval**3) * sum([cov.GetBinContent(iB, k) for k in range(1, vals.GetNbinsX() + 1)])
                                )
            )
            outcov.SetBinContent(jB, iB, outcov.GetBinContent(iB, jB))

    return outcov



def NormaliseAndPlot(tsk):
    inpath, iY, iV = tsk

    thedict, covmatnorm = calculateNormalisedValues(inpath, iY, iV)

    particlefidbinpath = inpath + "/" + iY + "/particlefidbinplots/"
    if not os.path.isdir(particlefidbinpath):
        os.system("mkdir -p " + particlefidbinpath)

    PlotParticleFidBinLevelResults(thedict, inpath, iY, iV, covmatnorm)
    return


def calculateNormalisedValues(inpath, iY, iV):
    varfile  = inpath + "/" + iY + "/" + iV + "/particleOutput.root"

    if not os.path.isfile(varfile):
        raise RuntimeError("FATAL: file {f} does not exist.".format(f = varfile))

    fvar = r.TFile.Open(varfile)
    dirfid = {}; dirvar = {};
    for key in fvar.GetListOfKeys():
        # print key.GetName()
        if iV not in key.GetName():
            if "nom" in key.GetName():
                newnam = iV + "_total" + ("Up" if "0" in key.GetName() else "Down")
                dirvar[newnam] = deepcopy(key.ReadObj().Clone(newnam))

                #for iB in range(1, dirvar[newnam].GetNbinsX() + 1):
                    #dirvar[newnam].SetBinContent(iB, dirvar[newnam].GetBinContent(iB) + dirvar[newnam].GetBinError(iB) if "Up" in newnam else dirvar[newnam].GetBinContent(iB) - dirvar[newnam].GetBinError(iB))
            pass

        else:
            dirvar[key.GetName()] = deepcopy(key.ReadObj().Clone(key.GetName()))
            # for iB in range(1, dirvar[key.GetName()].GetNbinsX() + 1):
            #     dirvar[key.GetName()].SetBinContent(iB, dirvar[key.GetName()].GetBinContent(iB) + dirvar[key.GetName()].GetBinError(iB) if "Up" in key.GetName() else dirvar[key.GetName()].GetBinContent(iB) - dirvar[key.GetName()].GetBinError(iB))
    covmatdict = {}
    # sys.exit()
    for key,histo in dirvar.iteritems():
        #print key, "CovMat_" + key.replace(iV + "_", "").replace(iV, "").replace("Up", "").replace("Down", "") if ("total" not in key and key != iV) else "CovMat"
        tmph = deepcopy(histo.Clone("tmpvar_" + key))
        for iB in range(1, histo.GetNbinsX() + 1):
            tmph.SetBinContent(iB, tmph.GetBinContent(iB) + (-tmph.GetBinError(iB) if "Down" in key else tmph.GetBinError(iB)))

        # covmatdict[key] = deepcopy(ep.getCovarianceFromVarv2(histo, tmph, iV))
        #covmatdict[key] = deepcopy(ep.getCovarianceFromVarv2(histo, tmph, iV)) if ((not "total" in key) or (key != iV)) else deepcopy(fvar.Get("CovMat").Clone("covtotal"))
        covmatdict[key] = deepcopy(fvar.Get("CovMat_" + key.replace(iV + "_", "").replace(iV, "").replace("Up", "").replace("Down", "") if ("total" not in key and key != iV) else "CovMat"))
        #covmatdict[key] = deepcopy(fvar.Get("CovMat").Clone("cov" + key.replace(iV, "")))
        fiduval = sum([histo.GetBinContent(i) for i in range(1, histo.GetNbinsX() + 1)])
        
        fiduunc = r.TMath.Sqrt(sum([covmatdict[key].GetBinContent(i, j) for i in range(1, covmatdict[key].GetNbinsX() + 1) for j in range(1, covmatdict[key].GetNbinsX() + 1) ]))
        
        htmp = r.TH1D(histo)
        for bin in range(1, htmp.GetNbinsX() + 1):
            htmp.SetBinContent(bin, fiduval)
            htmp.SetBinError(  bin, fiduunc)

        #if key == iV: print key, fiduval, fiduunc, fiduunc/fiduval
        dirfid[key.replace(iV, "Fiducial")] = deepcopy(htmp.Clone(key.replace(iV, "Fiducial")))
        del htmp
    fvar.Close(); del fvar

    #for key in dirvar:
        #print dirvar[key], dirfid[key]
    #sys.exit()

    #for key in covmatdict:
        #print key, covmatdict[key]
    #sys.exit()
    
    normcovmat = BibhuFunctionForCovs(covmatdict[iV], dirvar[iV], iV)
    dirfinal = {}
    for key in dirvar:
        name    = "Fiducial" if key == iV else "Fiducial" + key.replace(iV, "")
        dirfinal[key] = deepcopy(dirvar[key].Clone())
        dirfinal[key].Divide(dirfid[name])
        #print key
    
    BibhuFunction(dirvar, dirfid, dirfinal, covmatdict, iV)

    for key in dirfinal:
        name    = "Fiducial" if key == iV else "Fiducial" + key.replace(iV, "")
        # print "\n", key, name
        # for iB in range(1, dirfinal[key].GetNbinsX() + 1):
        #     print dirfinal[key].GetBinContent(iB), dirfinal[key].GetBinError(iB), dirfinal[key].GetBinError(iB)/dirfinal[key].GetBinContent(iB)
    #del dirfid

    allHistos     = { key.replace(iV + "_", "") : deepcopy(dirfinal[key]) for key in dirfinal}
    del allHistos[iV]
    allHistos[""] = deepcopy(dirfinal[iV])
    for iB in range(1, allHistos[""].GetNbinsX() + 1):
        allHistos[""].SetBinError(iB, 0.)

    savetfile = r.TFile(inpath + "/" + iY + "/" + iV + "/particlefidOutput.root", "recreate")
    for key in allHistos: allHistos[key].Write()
    savetfile.Close(); del savetfile

    normcovmat.Scale(1, "width")
    
    #savetfile3 = r.TFile(inpath + "/" + iY + "/" + varName + "/temp_allcovmats.root", "recreate")
    #covmatnorm.Write()
    #savetfile3.Close(); del savetfile3
    
    return allHistos, normcovmat


def giveMeOneComparisonNormalised(orighist, name, scalevalue, iV):
    outH = deepcopy(orighist.Clone(name if vl.spacingdict[name] == 0. else name + "tmp"))
    outH.SetLineWidth(0)
    # outH.SetLineColor(colour)
    outH.Scale(scalevalue)
    outH.SetMarkerSize(1.25)
    outH.SetMarkerColor(vl.comparisonColourDict[name])
    outH.SetMarkerStyle(vl.markersdict[name])
    htmp = deepcopy(outH.Clone("htmp"))
    fiduval = None; fiduunc = None

    fiduval = sum([htmp.GetBinContent(i) for i in range(1, htmp.GetNbinsX() + 1)])
    fiduunc = r.TMath.Sqrt(sum([htmp.GetBinError(i)**2 for i in range(1, htmp.GetNbinsX() + 1)]))

    for bin in range(1, htmp.GetNbinsX() + 1):
        htmp.SetBinContent(bin, fiduval)
        htmp.SetBinError(  bin, fiduunc)

    # if name == "tru": print fiduval

    thegooduncs = []
    for bin in range(1, outH.GetNbinsX() + 1):
        thegooduncs.append(r.TMath.Sqrt(outH.GetBinError(bin)**2 / fiduval**2 +
                                        outH.GetBinContent(bin)**2 * fiduunc**2 / fiduval**4 -
                                        2 * (outH.GetBinContent(bin) * outH.GetBinError(bin)**2) / fiduval**3 )
        )

    outH.Divide(htmp)

    for iB in range(1, outH.GetNbinsX() + 1):
        outH.SetBinError(iB, thegooduncs[iB - 1])

    outH.Scale(1, "width")
    del htmp

    if vl.spacingdict[name] != 0.:
        goodoutH = deepcopy(outH.Clone(name))
        smallestw = 9999
        for i in range(1, len(vl.varList[iV]["bins_particle"])):
            thedif = vl.varList[iV]["bins_particle"][i] - vl.varList[iV]["bins_particle"][i - 1]
            if thedif < smallestw:
                smallestw = thedif

        newbins = array("d", [el + vl.spacingdict[name]*smallestw/2 for el in vl.varList[iV]["bins_particle"]])
        goodoutH.SetBins(len(vl.varList[iV]["bins_particle"]) - 1, newbins)
        for iB in range(1, goodoutH.GetNbinsX() + 1):
            goodoutH.SetBinContent(iB, outH.GetBinContent(iB))
            goodoutH.SetBinError(iB, outH.GetBinError(iB))
        del outH
        return goodoutH
    else:
        return outH


def drawANormalisedCovMat(finalmat, inpath, iY, var):
    pathtothings = inpath + "/" + iY + "/" + var

    tdrstyle.setTDRStyle()
    finalmat.SetStats(False)
    finalmat.SetXTitle(vl.varList[var]['xaxis'])
    finalmat.SetYTitle(vl.varList[var]['xaxis'])
    finalmat.SetMarkerColor(r.kRed)
    finalmat.SetMarkerColor(r.kRed)

    if 'txtsize_covparticlefidbin' in vl.varList[var]:
        finalmat.SetMarkerSize(vl.varList[var]['txtsize_covparticlefidbin'])

    r.gStyle.SetPaintTextFormat("4.8f")
    r.gStyle.SetPadRightMargin(0.17)
    r.gStyle.SetPadTopMargin(0.05)
    r.gStyle.SetPadBottomMargin(0.1)
    r.gStyle.SetPadLeftMargin(0.12)
    finalmat.GetYaxis().SetTitleOffset(1.5)
    finalmat.GetXaxis().SetTitleOffset(1.1)
    finalmat.GetXaxis().SetTitleFont(43)
    finalmat.GetXaxis().SetTitleSize(22)
    finalmat.GetXaxis().SetLabelFont(43)
    finalmat.GetXaxis().SetLabelSize(22)
    finalmat.GetYaxis().SetTitleFont(43)
    finalmat.GetYaxis().SetTitleSize(22)
    finalmat.GetYaxis().SetLabelFont(43)
    finalmat.GetYaxis().SetLabelSize(22)
    r.gStyle.SetLabelFont(43, "XYZ")
    r.gStyle.SetLabelSize(22, "XYZ")

    c = r.TCanvas('c', "", 200, 10, 600, 600)

    txtanglestring = ""
    if ('txtangle_covparticlefidbin' in vl.varList[var]):
        txtanglestring = vl.varList[var]['txtangle_covparticlefidbin']

#    finalmat.Draw("colz text{s}".format(s = txtanglestring))
    finalmat.Draw("colz")
    r.gPad.Update()

    thelumi = vl.TotalLumi if iY == "run3" else vl.LumiDict[iY]
    if thelumi >= 10**2:
        CMS_lumi.lumi_13TeV = "%.0f fb^{-1}" %(thelumi)
    else:
        CMS_lumi.lumi_13TeV = "%.1f fb^{-1}" %(thelumi)

    CMS_lumi.extraText  = "Supplementary" + ' Preliminary' * vl.doPre
    CMS_lumi.lumi_sqrtS = '#sqrt{s} = 13 TeV'
    #CMS_lumi.cmsTextSize += 0.1
    CMS_lumi.CMS_lumi(r.gPad, 4, 0, 0.05)

    if doArXiv:
        arXiv = r.TLatex()
        arXiv.SetNDC()
        arXiv.SetTextAngle(0)
        arXiv.SetTextColor(r.kBlack)
        
        arXiv.SetTextFont(42)
        arXiv.SetTextAlign(31)
        arXiv.SetTextSize(0.6 * r.gPad.GetTopMargin())
        arXiv.DrawLatex(0.61, 1 - r.gPad.GetTopMargin() + 0.2 * r.gPad.GetTopMargin(), vl.arXivtext)

    plotsoutputpath = pathtothings + "/CovMatplots"
    if not os.path.isdir(plotsoutputpath):
        os.system('mkdir -p ' + plotsoutputpath)

    c.SaveAs(plotsoutputpath + "/Cov_{vr}_fidbin.png" .format(vr = var))
    c.SaveAs(plotsoutputpath + "/Cov_{vr}_fidbin.pdf" .format(vr = var))
    c.SaveAs(plotsoutputpath + "/Cov_{vr}_fidbin.eps" .format(vr = var))
    c.SaveAs(plotsoutputpath + "/Cov_{vr}_fidbin.root".format(vr = var))
    c.Close();
    return


#### POR BIN Y NORMALIZAOS
def PlotParticleFidBinLevelResults(thedict, inpath, iY, varName, covmatnorm):
    print("> Obtaining plots normalised to fiducial cross section and to bin's width for variable " + varName)
    plotspath  = inpath + "/" + iY + "/particlefidbinplots"
    folderpath = inpath + "/" + iY + "/" + varName

    for key in thedict:
        thedict[key].Scale(1, "width")

    savetfile2 = r.TFile(inpath + "/" + iY + "/" + varName + "/particlefidbinOutput.root", "recreate")
    for key in thedict: thedict[key].Write()
    savetfile2.Close(); del savetfile2


    addRelUncs = False
    if os.path.isdir(folderpath + "/sigextr_fit_combine/relativeuncs"):
        print("\t- Found relative uncertainty results. I'll plot them.")
        addRelUncs = True
        # indRelUncs = {}; gloRelUncs = {}
        # for iB in range(1, thedict[""].GetNbinsX() + 1):
        #     indRelUncs[iB] = vl.parseRelUncs(folderpath + "/sigextr_fit_combine/relativeuncs/individual/outputfit_r_tW_{i}.txt".format(i = iB - 1))
        #     gloRelUncs[iB] = vl.parseRelUncs(folderpath + "/sigextr_fit_combine/relativeuncs/global/outputfit_r_tW_{i}.txt".format(i = iB - 1))


    #nominal_withErrors = ep.propagateHisto(thedict, doSym = vl.doSym)
    nominal_withErrors = [deepcopy(thedict["totalUp"]), deepcopy(thedict["totalDown"])]
    #nominal_withErrors = [deepcopy(thedict[""]), deepcopy(thedict[""])]
    statOnlyList       = [deepcopy(thedict["statUp"]),  deepcopy(thedict["statDown"])]
    
    for iB in range(1, thedict[""].GetNbinsX() + 1):
        # if nominal_withErrors[1].GetBinError(iB) > nominal_withErrors[1].GetBinContent(iB):
        #     nominal_withErrors[1].SetBinError(iB, nominal_withErrors[1].GetBinContent(iB))

        #print nominal_withErrors[0].GetBinError(iB), nominal_withErrors[1].GetBinError(iB), statOnlyList[0].GetBinError(iB), statOnlyList[1].GetBinError(iB)
        thedict[""].SetBinError(iB, 1e-5)
#        thedict[""].SetBinError(iB, 0)

    #sys.exit()

    plot               = bp.beautifulUnfPlot(varName + "_particlefidbin", varName)
    plot.doRatio       = True
    plot.doFit         = False
    plot.plotspath     = inpath + "/" + iY + "/particlefidbinplots/"
    plot.displayedLumi = vl.TotalLumi if iY == "run3" else vl.LumiDict[iY]
    plot.doPreliminary = vl.doPre

    thedict[""].SetMarkerStyle(r.kFullCircle)
    thedict[""].SetLineColor(r.kBlack)
    thedict[""].SetMarkerSize(1)
    nominal_withErrors[0].SetFillColorAlpha(r.kOrange, 1.0)
    nominal_withErrors[0].SetLineColor(0)
    nominal_withErrors[0].SetFillStyle(1001)
    statOnlyList[0].SetFillColorAlpha(r.kGray + 1, 1.0)
    statOnlyList[0].SetLineColor(0)
    statOnlyList[0].SetFillStyle(1001)

    tex.saveLaTeXfromhisto(thedict[""], varName, path = inpath + "/" + iY + "/tables", errhisto = nominal_withErrors[0], ty = "particlefidbin")

    if "yaxismax_particlefidbin" in vl.varList[varName]:
        plot.yaxisuplimit = vl.varList[varName]["yaxismax_particlefidbin"]

    if "legpos_particlefidbin" in vl.varList[varName]:
        legloc = vl.varList[varName]["legpos_particlefidbin"]
    else:
        legloc = "TR"

    thelumi = vl.TotalLumi if iY == "run3" else vl.LumiDict[iY]
    scaleval = 1/thelumi/1000

    if not os.path.isfile(inpath + "/" + iY + "/" + varName + '/particle.root'):
        raise RuntimeError('The rootfile with the generated information does not exist.')
    

    tmptfile = r.TFile.Open(inpath + "/" + iY + "/" + varName + '/particle.root')


    tru                = giveMeOneComparisonNormalised(tmptfile.Get('x_tw'),                         "tru",                scaleval, varName)
    tru_DS             = giveMeOneComparisonNormalised(tmptfile.Get('x_twds'),                       "tru_DS",             scaleval, varName)
    #tru_herwig         = giveMeOneComparisonNormalised(tmptfile.Get('x_twherwig'),                   "tru_herwig",         scaleval, varName)
    tru_aMC_dr         = giveMeOneComparisonNormalised(tmptfile.Get('x_twamcatnlo_dr'),              "tru_aMC_dr",         scaleval, varName)
    tru_aMC_dr2        = giveMeOneComparisonNormalised(tmptfile.Get('x_twamcatnlo_dr2'),             "tru_aMC_dr2",        scaleval, varName)
    tru_aMC_ds         = giveMeOneComparisonNormalised(tmptfile.Get('x_twamcatnlo_ds'),              "tru_aMC_ds",         scaleval, varName)
    tru_aMC_ds_runn    = giveMeOneComparisonNormalised(tmptfile.Get('x_twamcatnlo_ds_runningBW'),    "tru_aMC_ds_runn",    scaleval, varName)
    tru_aMC_ds_IS      = giveMeOneComparisonNormalised(tmptfile.Get('x_twamcatnlo_ds_is'),           "tru_aMC_ds_IS",      scaleval, varName)
    tru_aMC_ds_IS_runn = giveMeOneComparisonNormalised(tmptfile.Get('x_twamcatnlo_ds_is_runningBW'), "tru_aMC_ds_IS_runn", scaleval, varName)

    tmptfile.Close(); del tmptfile

    savetfile2 = r.TFile(inpath + "/" + iY + "/" + varName + "/particlefidbinOutput.root", "update")
    tru.Write()
    tru_DS.Write()
    #tru_herwig.Write()
    tru_aMC_dr        .Write()
    tru_aMC_dr2       .Write()
    tru_aMC_ds        .Write()
    tru_aMC_ds_runn   .Write()
    tru_aMC_ds_IS     .Write()
    tru_aMC_ds_IS_runn.Write()
    savetfile2.Close()

    #themaxs = []
    #for el in [tru, tru_DS, tru_herwig, thedict[""], nominal_withErrors[0], nominal_withErrors[1]]:
        #themaxs.append(vl.getAConservativeMaximum(el))
    #tmpval = max(themaxs)
    """
    print tmpval
    nominal_withErrors[0].SetMaximum(tmpval)
    nominal_withErrors[1].SetMaximum(tmpval)
    statOnlyList[0].SetMaximum(tmpval)
    statOnlyList[1].SetMaximum(tmpval)
    """

    plot.addHisto(nominal_withErrors, 'A2',     'Total unc.',    'F', 'total')
    plot.addHisto(statOnlyList,       '2,same', 'Stat unc.',     'F', 'stat')
    
    plot.addHisto(tru,                'P,same', 'tW PH DR + P8',          'P', 'mc')
    plot.addHisto(tru_DS,             'P,same', 'tW PH DS + P8',          'P', 'mc')
    #plot.addHisto(tru_herwig,         'P,same', 'tW PH DR + H7',          'P', 'mc')
    plot.addHisto(tru_aMC_dr,         'P,same', 'tW aMC DR + P8',         'P', 'mc')
    plot.addHisto(tru_aMC_dr2,        'P,same', 'tW aMC DR2 + P8',        'P', 'mc')
    plot.addHisto(tru_aMC_ds,         'P,same', 'tW aMC DS + P8',         'P', 'mc')
    plot.addHisto(tru_aMC_ds_runn,    'P,same', 'tW aMC DS dyn. + P8',    'P', 'mc')
    plot.addHisto(tru_aMC_ds_IS,      'P,same', 'tW aMC DS IS + P8',      'L', 'mc')
    plot.addHisto(tru_aMC_ds_IS_runn, 'P,same', 'tW aMC DS IS dyn. + P8', 'L', 'mc')



    plot.addHisto(thedict[""],        'P,E,same{s}'.format(s = ",X0" if "equalbinsunf" in vl.varList[varName] else ""), vl.labellegend, 'P,same' + (",L" if not "equalbinsunf" in vl.varList[varName] else "") , 'data', redrawaxis = True)

    #print varName, 'P,same{s}'.format(s = ",X0" if "equalbinsunf" in vl.varList[varName] else "")
    #plot.addHisto(nominal_withErrors, 'A,same',     'Total unc.',    'F', 'bullshit')

#    plot.saveCanvasv2(legloc, leg = False)
    plot.saveCanvasv2(legloc)

    #plot.addHisto(nominal_withErrors, 'E2',     'Uncertainty',   'F')
    #plot.addHisto(thedict[""], 'P,same',vl.labellegend,'PE', 'data')
    #plot.saveCanvas(legloc)
    del plot

    plot2       = bp.beautifulUnfPlot(varName + 'uncs_particlefidbin', varName)
    plot2.doFit         = False
    plot2.displayedLumi = vl.TotalLumi if iY == "run3" else vl.LumiDict[iY]
    plot2.doPreliminary = vl.doPre
    plot2.plotspath     = inpath + "/" + iY + "/particlefidbinplots/"

    yaxismax_particlefidbinunc = 1.4
    if "yaxismax_particlefidbinunc" in vl.varList[varName]:
        yaxismax_particlefidbinunc = vl.varList[varName]["yaxismax_particlefidbinunc"]

    uncListorig, hincstat, hincsyst, hincmax = ep.drawTheRelUncPlotv2(nominal_withErrors, thedict, plot2, yaxismax_particlefidbinunc, doSym =  vl.doSym)

    if "legpos_particlefidbinunc" in vl.varList[varName]:
        unclegpos = vl.varList[varName]["legpos_particlefidbinunc"]
    else:
        unclegpos = "TR"

    plot2.saveCanvas(unclegpos)
    del plot2, nominal_withErrors#, tru, tru_DS, tru_herwig, tru_aMC_dr, tru_aMC_dr2, tru_aMC_ds, tru_aMC_ds_runn

    #### COVMAT
    drawANormalisedCovMat(covmatnorm, inpath, iY, varName)
    savetfile3 = r.TFile(inpath + "/" + iY + "/" + varName + "/CovMat_particlefidbin_v2.root", "recreate")
    covmatnorm.Write()
    savetfile3.Close(); del savetfile3
    return




if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage = "python nanoAOD_checker.py [options]", description = "Checker tool for the outputs of nanoAOD production (NOT postprocessing)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inpath',      '-i',  metavar = 'inpath',     dest = "inpath",      required = False, default = "./temp/differential")
    parser.add_argument('--year',        '-y',  metavar = 'year',       dest = "year",        required = False, default = "all")
    parser.add_argument('--variable',    '-v',  metavar = 'variable',   dest = "variable",    required = False, default = "all")
    parser.add_argument('--extraArgs',   '-e',  metavar = 'extra',      dest = "extra",       required = False, default = "")
    parser.add_argument('--nthreads',    '-j',  metavar = 'nthreads',   dest = "nthreads",    required = False, default = 0, type = int)
    parser.add_argument('--pretend',     '-p',  action  = "store_true", dest = "pretend",     required = False, default = False)


    args        = parser.parse_args()
    year        = args.year
    extra       = args.extra
    nthreads    = args.nthreads
    pretend     = args.pretend
    inpath      = args.inpath
    variable    = args.variable

    #print("\n===== Unfolding procedures: Response matrices & ROOT files production =====")
    #print("> Setting binning, paths, and other details...")

    #print "\n> Drawing matrices and writing ROOT file (old one will be overwritten!)."

    #### First, find the tasks
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
            tasks.append( (inpath, iY, iV) )

    #tasks = [ (inpath, "2016", "Lep1Lep2Jet1MET_Mt") ]

    if nthreads > 1:
        pool = Pool(nthreads)
        pool.map(NormaliseAndPlot, tasks)
        pool.close()
        pool.join()
    else:
        for tsk in tasks:
            NormaliseAndPlot(tsk)
