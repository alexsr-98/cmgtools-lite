from copy import deepcopy
import ROOT as r
import warnings as wr
import os, sys, math, array, argparse
import numpy as np
from scipy.stats import chisquare, chi2
from scipy.linalg import eigh

sys.path.append('{cmsswpath}/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/differential/'.format(cmsswpath = os.environ['CMSSW_BASE']))
import beautifulUnfoldingPlots as bp
import errorPropagator as ep
import getLaTeXtable as tex
import varList as vl

r.gROOT.SetBatch(True)
vl.SetUpWarnings()


def getNumpyArray(histo):
    ret = None
    if   isinstance(histo, r.TH2):
        ret = np.zeros( (histo.GetNbinsX(), histo.GetNbinsX()), dtype = np.double )
        for i in range(1, histo.GetNbinsX() + 1):
            for j in range(1, histo.GetNbinsX() + 1):
                ret[i-1, j-1] = histo.GetBinContent(i, j)

    elif isinstance(histo, r.TH1):
        ret = np.array([histo.GetBinContent(bin) for bin in range(1, histo.GetNbinsX() + 1)], dtype = np.double)
    elif isinstance(histo, r.TMatrixD):
        ret = np.zeros( (histo.GetNrows(), histo.GetNcols()), dtype = np.double )
        for i in range(histo.GetNrows()):
            for j in range(histo.GetNcols()):
                ret[i, j] = histo(i,j)
    else:
        raise RuntimeError("FATAL: unsupported type " + str(type(histo)) + ".")

    return ret


def Chi2TestForMultivarNormal(w1, V1, w2, V2):
    if not isinstance(w1, np.ndarray):
        w1 = getNumpyArray(w1)
    if not isinstance(w2, np.ndarray):
        w2 = getNumpyArray(w2)

    if isinstance(V1, r.TH1) and not isinstance(V1, r.TH2):
        V1 = np.diag(np.array([V1.GetBinError(bin)**2 for bin in range(1, V1.GetNbinsX() + 1)], dtype = np.double))
    elif not isinstance(V1, np.ndarray):
        V1 = getNumpyArray(V1)

        #### COSUQUES
        #V1_mat = r.TMatrixD(V1.GetNbinsX() + 2, V1.GetNbinsY() + 2, V1.GetArray(), "D")

        #V1_mat2 = V1_mat.GetSub(1, V1.GetNbinsX(), 1, V1.GetNbinsY()).Clone("V1_mat2")

        #del V1
        ##V1_mat2.SetTol(1e-23)
        #V1_mat2.Invert()
        #V1 = getNumpyArray(V1_mat2)

    if isinstance(V2, r.TH1) and not isinstance(V2, r.TH2):
        V2 = np.diag(np.array([V2.GetBinError(bin)**2 for bin in range(1, V2.GetNbinsX() + 1)], dtype = np.double))
    elif not isinstance(V2, np.ndarray):
        V2 = getNumpyArray(V2)


    if len(w1) != len(w2):
        raise RuntimeError("The length of both arrays is not the same.")
    if V1.shape != V2.shape:
        raise RuntimeError("The shape of the covariance matrices is not the same.")
    if V1.shape[0] != V1.shape[0] or V2.shape[1] != V2.shape[1]:
        raise RuntimeError("One (or both) covariance matrix is (are) not squared.")
    if len(w1) != V1.shape[0] or len(w2) != V2.shape[0]:
        raise RuntimeError("The length of one (or both) array is not the same as the number of lines or columns of its respective covariance matrix.")

    #### FORCED DIAGONALISATION
    #V1 = np.diag([ V1[el, el] for el in range(len(w1)) ])


    #### PA ELIMINAR COSAS
    #delInd = -3

    #w1 = np.delete(w1, delInd, 0)
    #w2 = np.delete(w2, delInd, 0)
    #V1 = np.delete(np.delete(V1, delInd, 0), delInd, 1)
    #V2 = np.delete(np.delete(V2, delInd, 0), delInd, 1)


    #### TRANSFORMACION
    vals, cambmat = eigh(V1)

    #print vals, min(vals)
    #print cambmat

    if min(vals) < 1e-14:
        print("WARNING: singular matrix. Removing affected rows/cols.")
        icambmat = np.linalg.inv(cambmat)
        V1 = icambmat.dot(V1).dot(cambmat)
        V1 = np.where(V1 < 1e-14, 0., V1)
        V2 = icambmat.dot(V2).dot(cambmat)
        w1 = (icambmat.dot(w1.transpose())).transpose()
        w2 = (icambmat.dot(w2.transpose())).transpose()

        # Ahora quitamos lo que nos jode la marrana
        thewheres = np.where(vals < 1e-14)[0]
        if len(thewheres) > 1:
            raise RuntimeError("FATAL: more than one zero eigenvalue found! WE ARE NOT READY FOR THIIIISSSS!")
        delInd = thewheres[0]
        w1 = np.delete(w1, delInd, 0)
        w2 = np.delete(w2, delInd, 0)
        V1 = np.delete(np.delete(V1, delInd, 0), delInd, 1)
        V2 = np.delete(np.delete(V2, delInd, 0), delInd, 1)

    #sys.exit()

    dof = len(w1) - 1

    #### COSA:
    #if min(vals) < 1e-14:
        #dof += 1


    #V1 = 1e6 * V1
    
    #print np.linalg.cond(V1), np.linalg.cond(V2)
    if np.linalg.cond(V1) > 1e5:
        print("WARNING: condition number of data's covariance matrix high ({cd}). It might be non-invertible!".format(cd = np.linalg.cond(V1)))

    invV1 = np.linalg.inv(V1); invV2 = np.linalg.inv(V2)
    #invV1 = V1; invV2 = np.linalg.inv(V2)
    W1 = sum(w1); W2 = sum(w2)
    
    #print "W1:",  W1
    #print "W2:",  W2
    #print "dof:", dof
    #print "pruebina:", invV1.dot(V1)[2, 2]
    
    # 1) Calculation of MLE of the vector of probabilities of each bin
    p    = np.linalg.inv(W1**2 * invV1 + W2**2 * invV2).dot( W1 * invV1.dot(w1) + W2 * invV2.dot(w2) )

    # 2) Calculation of the Pearson's chi-squared test statistic
    est  = ((w1 - W1 * p).dot(invV1)).dot((w1 - W1 * p).transpose()) + ((w2 - W2 * p).dot(invV2)).dot((w2 - W2 * p).transpose())
    #est  = ((w1 - w2).dot(invV1)).dot((w1 - w2).transpose())
    
    # 3) Calculation of the corresponding p-value
    pval = 1 - chi2.cdf(est, dof)
    return pval, est


def GiveMeMyGoodGOFTests(tsk):
    inpath, iY, var, ty = tsk
    pathtothings = inpath + "/" + iY + "/" + var
    print("\n====== Performing tests for variable", var, "and with type", ty, "\n")
    if   ty == "particle":
        f1 = r.TFile.Open(pathtothings + "/particleOutput.root",        "read")
        f2 = r.TFile.Open(pathtothings + "/particleOutput.root",       "read")
        covname = "FullCovMat"
    elif ty == "particlefidbin":
        f1 = r.TFile.Open(pathtothings + "/particlefidbinOutput.root",  "read")
        f2 = r.TFile.Open(pathtothings + "/CovMat_particlefidbin_v2.root", "read")
        covname = "CovMat_fidbin"
    else:
        raise RuntimeError("FATAL: type {t} is not implemented.".format(t = ty))
    
    hData     = deepcopy(f1.Get(var).Clone("hData"))
    htW        = deepcopy(f1.Get("tw").            Clone("htW"))
    hDS          = deepcopy(f1.Get("twds").         Clone("hDS"))
    hHerwig      = deepcopy(f1.Get("twherwig").     Clone("hHerwig"))
    haMC_dr      = deepcopy(f1.Get("twamcatnlo_dr").     Clone("haMC_dr"))
    haMC_dr2     = deepcopy(f1.Get("twamcatnlo_dr2").    Clone("haMC_dr2"))
    haMC_ds      = deepcopy(f1.Get("twamcatnlo_ds").     Clone("haMC_ds"))
    haMC_ds_runn = deepcopy(f1.Get("twamcatnlo_ds_runningBW").Clone("haMC_ds_runn"))
        
    hDataCov     = deepcopy(f2.Get(covname).Clone("hDataCov"))
    
    f1.Close(); f2.Close()
    
    hData.SetFillColor(r.kWhite)
    hData.SetLineWidth(2)
    hData.SetLineColor(r.kBlack)
    hData.SetLineStyle(2)

    thelumi = vl.TotalLumi if iY == "run3" else vl.LumiDict[iY]
    scaleval = 1/thelumi/1000

    htW.SetFillColor(r.kWhite)
    htW.SetLineWidth(2)
    htW.SetLineColor(r.kRed)
    
    hDS.SetFillColor(r.kWhite)
    hDS.SetLineWidth(2)
    hDS.SetLineColor(r.kGreen)
    hDS.GetXaxis().SetTitle(vl.varList[var]['xaxis'])
    #if ty == "particle":
        #for el in [hDR, hDS, hHerwig, haMC_dr, haMC_dr2, haMC_ds, haMC_ds_runn]:
            #el.Scale(scaleval)
    
    #if ty == "particlefidbin":
        #for iB in range(1, hData.GetNbinsX() + 1):
            #print "\nbin", iB
            #for el in [hData, hDR, hDS, hHerwig, haMC_dr, haMC_dr2, haMC_ds, haMC_ds_runn]:
                #print el.GetBinContent(iB)

        #sys.exit()
    #covmat = np.zeros( (hData.GetNbinsX(), hData.GetNbinsX()), dtype = np.double )
    
    #for i in range(1, hData.GetNbinsX() + 1):
        #for j in range(1, hData.GetNbinsX() + 1):
            ##if i != j: continue
            #covmat[i-1, j-1] = hDataCov.GetBinContent(i, j)
    
    tWcovmat          = np.diag(np.array([htW.GetBinError(bin)**2        for bin in range(1, htW.GetNbinsX() + 1)],        dtype = np.double))
    DScovmat          = np.diag(np.array([hDS.GetBinError(bin)**2          for bin in range(1, hDS.GetNbinsX() + 1)],          dtype = np.double))
    Herwigcovmat      = np.diag(np.array([hHerwig.GetBinError(bin)**2      for bin in range(1, hHerwig.GetNbinsX() + 1)],      dtype = np.double))
    aMC_drcovmat      = np.diag(np.array([haMC_dr.GetBinError(bin)**2      for bin in range(1, haMC_dr.GetNbinsX() + 1)],      dtype = np.double))
    aMC_dr2covmat     = np.diag(np.array([haMC_dr2.GetBinError(bin)**2     for bin in range(1, haMC_dr2.GetNbinsX() + 1)],     dtype = np.double))
    aMC_dscovmat      = np.diag(np.array([haMC_ds.GetBinError(bin)**2      for bin in range(1, haMC_ds.GetNbinsX() + 1)],      dtype = np.double))
    aMC_ds_runncovmat = np.diag(np.array([haMC_ds_runn.GetBinError(bin)**2 for bin in range(1, haMC_ds_runn.GetNbinsX() + 1)], dtype = np.double))
    
    
    #covmat         *= 10000**2
    #DRcovmat       *= 10000**2
    #DScovmat       *= 10000**2
    #aMCatNLOcovmat *= 10000**2
    #hDR.Scale(10000)
    #hDS.Scale(10000)
    #haMCatNLO.Scale(10000)
    #hData.Scale(10000)
    #print "\ncovmat de DR\n", DRcovmat
    
    coses = {}
    coses["tW"] = {
        "hist"   : deepcopy(htW),
        "covmat" : deepcopy(tWcovmat)
    }
    coses["DS"] = {
        "hist"   : deepcopy(hDS),
        "covmat" : deepcopy(DScovmat)
    }
    coses["Herwig"] = {
        "hist"   : deepcopy(hHerwig),
        "covmat" : deepcopy(Herwigcovmat)
    }
    coses["aMC_dr"] = {
        "hist"   : deepcopy(haMC_dr),
        "covmat" : deepcopy(aMC_drcovmat)
    }
    coses["aMC_dr2"] = {
        "hist"   : deepcopy(haMC_dr2),
        "covmat" : deepcopy(aMC_dr2covmat)
    }
    coses["aMC_ds"] = {
        "hist"   : deepcopy(haMC_ds),
        "covmat" : deepcopy(aMC_dscovmat)
    }
    coses["aMC_ds_runn"] = {
        "hist"   : deepcopy(haMC_ds_runn),
        "covmat" : deepcopy(aMC_ds_runncovmat)
    }
    #for i in range(1, coses["DR"]["hist"].GetNbinsX() + 1):
        #print coses["DR"]["hist"].GetBinContent(
    
    vdata = np.array([hData.GetBinContent(bin) for bin in range(1, hData.GetNbinsX() + 1)], dtype = np.double)
    
    #V1preinv = getNumpyArray(hDataCov)
    
    #for el in eig(V1preinv)[0]:
        #print el
    #print ""
    #tmpsum = 0.
    #for iB in range(1, hDataCov.GetNbinsX() + 1):
        #for jB in range(1, hDataCov.GetNbinsX() + 1):
            #tmpsum += hDataCov.GetBinContent(iB, jB)
    
    #print "LASUMA", tmpsum
    #print ""
    for key in coses:
        vtemp = np.array([coses[key]["hist"].GetBinContent(bin) for bin in range(1, coses[key]["hist"].GetNbinsX() + 1)], dtype = np.double)
        pv, es = Chi2TestForMultivarNormal(vdata, hDataCov, vtemp, coses[key]["covmat"])
        coses[key]["p-value"]   = pv
        coses[key]["statistic"] = es
    
    if not os.path.isdir(pathtothings + "/goftestplots"):
        os.system("mkdir -p " + pathtothings + "/goftestplots")
    
    plot           = bp.beautifulUnfPlot(var + "_GOF_" + ty, var)
    plot.doRatio   = False
    plot.doFit     = False
    plot.plotspath = pathtothings + "/goftestplots"
    
    #hDS.GetYaxis().SetRangeUser(0, 1.1 * max( [hDS.GetMaximum(), hDR.GetMaximum(), haMCatNLO.GetMaximum(), hData.GetMaximum()] ))
    htW.GetYaxis().SetRangeUser(0, 1.1 * max( [htW.GetMaximum(), hData.GetMaximum()] ))
    
    plot.addHisto(htW,          'hist',      'bb4l PH+P8',                   'L', 'mc')
    plot.addHisto(hDS,          'hist,same', 'tW DS + t#bar{t} PH+P8',       'L', 'mc')
    plot.addHisto(hHerwig,      'hist,same', 'tW DR + t#bar{t} PH+H7',       'L', 'mc')
    plot.addHisto(haMC_dr,      'hist,same', 'tW DR + t#bar{t} aMC+P8',      'L', 'mc')
    plot.addHisto(haMC_dr2,     'hist,same', 'tW DR2 + t#bar{t} aMC+P8',     'L', 'mc')
    plot.addHisto(haMC_ds,      'hist,same', 'tW DS + t#bar{t} aMC+P8',      'L', 'mc')
    plot.addHisto(haMC_ds_runn, 'hist,same', 'tW DS dyn. + t#bar{t} aMC+P8', 'L', 'mc')
    plot.addHisto(hData,        'hist,same', 'Data',                         'L', 'mc')
    
    thelatxs = []
    
    thelatxs.append(r.TLatex(0.65, 0.650, '#scale[0.4]{DR - p-val.:                %4.10f}'%coses["tW"]["p-value"]))
    thelatxs.append(r.TLatex(0.65, 0.625, '#scale[0.4]{DR - stat.:                 %4.10f}'%coses["tW"]["statistic"]))
    thelatxs.append(r.TLatex(0.65, 0.600, '#scale[0.4]{DS - p-val.:                %4.10f}'%coses["DS"]["p-value"]))
    thelatxs.append(r.TLatex(0.65, 0.575, '#scale[0.4]{DS - stat.:                 %4.10f}'%coses["DS"]["statistic"]))
    thelatxs.append(r.TLatex(0.65, 0.550, '#scale[0.4]{Herwig - p-val.:            %4.10f}'%coses["Herwig"]["p-value"]))
    thelatxs.append(r.TLatex(0.65, 0.525, '#scale[0.4]{Herwig - stat.:             %4.10f}'%coses["Herwig"]["statistic"]))
    thelatxs.append(r.TLatex(0.65, 0.500, '#scale[0.4]{aMC_dr - p-val.:            %4.10f}'%coses["aMC_dr"]["p-value"]))
    thelatxs.append(r.TLatex(0.65, 0.475, '#scale[0.4]{aMC_dr - stat.:             %4.10f}'%coses["aMC_dr"]["statistic"]))
    thelatxs.append(r.TLatex(0.65, 0.450, '#scale[0.4]{aMC_dr2 - p-val.:           %4.10f}'%coses["aMC_dr2"]["p-value"]))
    thelatxs.append(r.TLatex(0.65, 0.425, '#scale[0.4]{aMC_dr2 - stat.:            %4.10f}'%coses["aMC_dr2"]["statistic"]))
    thelatxs.append(r.TLatex(0.65, 0.400, '#scale[0.4]{aMC_ds - p-val.:            %4.10f}'%coses["aMC_ds"]["p-value"]))
    thelatxs.append(r.TLatex(0.65, 0.375, '#scale[0.4]{aMC_ds - stat.:             %4.10f}'%coses["aMC_ds"]["statistic"]))
    thelatxs.append(r.TLatex(0.65, 0.350, '#scale[0.4]{aMC_ds_runn - p-val.:       %4.10f}'%coses["aMC_ds_runn"]["p-value"]))
    thelatxs.append(r.TLatex(0.65, 0.325, '#scale[0.4]{aMC_ds_runn - stat.:        %4.10f}'%coses["aMC_ds_runn"]["statistic"]))
    
    for el in thelatxs:
        el.SetNDC(True)
    
    for el in thelatxs:
        el.Draw("same")
    
    plot.saveCanvas("TR", leg = False)
    del plot
    
    outfile = open(pathtothings + "/goftestplots/" + var + "_goftests_" + ty + ".txt", "w")
    outtxt = ""
    outtxt += "Results of GOF tests\n"
    outtxt += "Variable: {vr} \n".format(vr = var)
    outtxt += "=========================================\n"
    
    print("\n")
    for key in ["tW", "DS", "Herwig", "aMC_dr", "aMC_dr2", "aMC_ds", "aMC_ds_runn"]:
        outtxt += key + " / p-value: "        + str(coses[key]["p-value"])   + "\n"
        outtxt += key + " / test statistic: " + str(coses[key]["statistic"]) + "\n"
        print(key + ' - p-val.:', coses[key]["p-value"])
        print(key + ' - stat.:',  coses[key]["statistic"])

    outfile.write(outtxt)
    outfile.close(); del outfile;

    del coses
    return



if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage = "python nanoAOD_checker.py [options]", description = "Checker tool for the outputs of nanoAOD production (NOT postprocessing)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inpath',    '-i', metavar = 'inpath',     dest = "inpath",   required = False, default = "./temp/differential")
    parser.add_argument('--year',      '-y', metavar = 'year',       dest = "year",     required = False, default = "all")
    parser.add_argument('--variable',  '-v', metavar = 'variable',   dest = "variable", required = False, default = "all")
    parser.add_argument('--extraArgs', '-e', metavar = 'extra',      dest = "extra",    required = False, default = "")
    parser.add_argument('--nthreads',  '-j', metavar = 'nthreads',   dest = "nthreads", required = False, default = 0, type = int)
    parser.add_argument('--pretend',   '-p', action  = "store_true", dest = "pretend",  required = False, default = False)


    args     = parser.parse_args()
    year     = args.year
    extra    = args.extra
    nthreads = args.nthreads
    pretend  = args.pretend
    inpath   = args.inpath
    variable = args.variable

    #### First, find the tasks
    tasks = []
    if year == "all":
        if variable == "all":
            theyears = []
            presentyears = next(os.walk(inpath))[1]
            #if "2022" in presentyears:
            #    theyears.append("2022")
            #if "2022PostEE" in presentyears:
            #    theyears.append("2022PostEE")
            if "run3" in presentyears:
                theyears.append("run3")
            for iY in theyears:
                #if iY != "run2": continue
                thevars = next(os.walk(inpath + "/" + iY))[1]
                for iV in thevars:
                    if any(el in iV for el in vl.vetolist + ["Fiducial"]): continue

                    for t in ["particle","particlefidbin"]:
                    #for t in ["particlefidbin"]:
                    #for t in ["particle"]:
                        tasks.append( (inpath, iY, iV, t) )

    #tasks = [ (inpath, "2016", "Lep1_Pt") ]

    if nthreads > 1:
        pool = Pool(nthreads)
        pool.map(GiveMeMyGoodGOFTests, tasks)
        pool.close()
        pool.join()
    else:
        for tsk in tasks:
            GiveMeMyGoodGOFTests(tsk)

    #tex.getgoftestsLaTeXtable(["Lep1Lep2_DPhi"],
    #                            inpath + "/2016/tables/",
    #                            inpath + "/2016/")

    if year == "all" or year == "run3":
        print("\nCreating LaTeX table with all the information...\n")

        cleanedvars = []
        for var in vl.varList["Names"]["Variables"]:
            if var != "Fiducial":
                cleanedvars.append(var)
        tex.getgoftestsLaTeXtable(cleanedvars,
                                  inpath + "/run3/tables/",
                                  inpath + "/run3/")
        tex.getgoftestsLaTeXtable(cleanedvars,
                                  inpath + "/run3/tables/",
                                  inpath + "/run3/",
                                  ty="particle")
        
    