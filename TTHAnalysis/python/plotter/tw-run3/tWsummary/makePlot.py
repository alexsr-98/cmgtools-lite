import ROOT
from prediction import Prediction, PredictionValue
from measurement import Measurement
from measurements_LHC import measurementData
from helpers import makeText, makeLogo, makeSubtitle

# https://gitlab.cern.ch/cms-analysis/top/SummaryPlots/-/tree/master/tX/xsecCurve?ref_type=heads
# https://gitlab.cern.ch/lhctopwg/summary-plots/-/tree/TTXS_ratio_refs/TTbarXSecvsEnergy_ratio 

def main(internal = False):
    ROOT.gROOT.SetBatch(ROOT.kTRUE)

    ############################################################################
    # Get theory curves
    predictions = [
        #Prediction("LHC_pdf4lhc21.dat", "NNLO+NNLL, PDF4LHC21 (pp)", ROOT.kSpring,  range=(2.3, 14.5), range_ratio=(2.8, 14.5)),
        Prediction("kidonakisaNNNLO.dat" , "aNNLO+aN^{3}LL, PDF4LHC21 (pp), JHEP 05 (2021) 278"   , ROOT.kAzure+7, range=(2.3, 14.5), range_ratio=(2.8, 14.5)),
        #Prediction("LHC_pdf4lhc21_tW.dat", "NNLO+NNLL, PDF4LHC21??? (pp), PRD 82 (2010) 054018", ROOT.kSpring , range=(2.3, 14.5), range_ratio=(2.8, 14.5)),
        #Prediction("Tevatron_pdf4lhc21.dat", "NNLO+NNLL, PDF4LHC21 (p#bar{p})", ROOT.kAzure+7, range=(1.5, 14.5), range_ratio=(1.5, 2.7)),
        # Prediction("Tevatron_nnpdf3p0.dat", "NNLO+NNLL, NNPDF3.0 (p#bar{p})", ROOT.kAzure+7, range=(1.5, 14.5), range_ratio=(1.5, 2.7)),
    ]
    i_lhc = 0
    i_tevatron = 1


    ############################################################################
    # Get theory values for specific sqrt(s) for inset
    # p_13TeV = PredictionValue(
    #     name = "PDF4LHC21",
    #     xs = 833.97,
    #     scaleUp = 854.43-833.97,
    #     scaleDown = 803.82-833.97,
    #     pdfUp = 20.94,
    #     pdfDown = 20.94,
    #     alphasUp = 0., # already included in pdf uncert
    #     alphasDown = 0., # already included in pdf uncert
    #     xmin = 12.8,
    #     xmax = 13.2,
    #     color = ROOT.kGreen
    #     )
    #
    # p_13p6TeV = PredictionValue(
    #     name = "PDF4LHC21",
    #     xs = 923.64,
    #     scaleUp = 946.87-923.64,
    #     scaleDown = 890.80-923.64,
    #     pdfUp = 22.82,
    #     pdfDown = 22.82,
    #     alphasUp = 0., # already included in pdf uncert
    #     alphasDown = 0., # already included in pdf uncert
    #     xmin = 13.4,
    #     xmax = 13.8,
    #     color = ROOT.kGreen
    #     )
    ############################################################################
    # Read Measurements
    # put measurements in measurement class, make list of references, check if some measurement is preliminary (*)
    hasPreliminary = False
    measurements = []
    refs = []
    refnr = 0
    for data in measurementData:
        print("Read measurement data:", data)
        (name, sqrts, lumi, xs, uncertUp, uncertDown, markerstyle, markersize, color, xoffset, ref) = data
        if "*" in name:
            hasPreliminary = True
        if ref not in refs: # make sure to not list the same ref twice
            refnr += 1
            refs.append( ref )
        m = Measurement(name, sqrts, lumi, xs, uncertUp, uncertDown, markerstyle, markersize, color, xoffset, refnr, ref)
        measurements.append(m)

    ############################################################################
    # Set axes ranges
    xmin = 5.0
    xmax = 14.1
    # logY    
    #ymin = 5.0
    #ymax = 200.
    # Lin Y
    ymin = 0.001
    ymax = 110.

    ymin_ratio = 0.65
    ymax_ratio = 1.35

    xmin_inset = 12.7
    xmax_inset = 13.9
    ymin_inset = 760.0
    ymax_inset = 970.0

    ############################################################################
    # Create a dummy to get axes ranges right
    dummy = ROOT.TGraph(2)
    dummy.SetPoint(0, xmin, ymin)
    dummy.SetPoint(1, xmax, ymax)
    dummy.SetMarkerSize(0.0)
    dummy.SetTitle("")
    dummy.GetXaxis().SetLabelSize(0.0)
    dummy.GetXaxis().SetTitle(" ")
    dummy.GetYaxis().SetTitle("Inclusive tW cross section (pb)")
    dummy.GetYaxis().SetLabelSize(0.05)
    dummy.GetYaxis().SetTitleSize(0.06)
    dummy.GetYaxis().SetTitleOffset(0.8)

    dummy_ratio = ROOT.TGraph(2)
    dummy_ratio.SetPoint(0, xmin, ymin_ratio)
    dummy_ratio.SetPoint(1, xmax, ymax_ratio)
    dummy_ratio.SetMarkerSize(0.0)
    dummy_ratio.SetTitle("")
    dummy_ratio.GetXaxis().SetTitle("#sqrt{s} (TeV)")
    dummy_ratio.GetYaxis().SetTitle("#splitline{  Ratio to}{Prediction}")
    dummy_ratio.GetXaxis().SetTitleSize(0.16)
    dummy_ratio.GetXaxis().SetLabelSize(0.16)
    dummy_ratio.GetXaxis().SetLabelOffset(0.02)
    dummy_ratio.GetYaxis().SetTitleSize(0.15)
    dummy_ratio.GetYaxis().SetTitleOffset(0.28)
    dummy_ratio.GetYaxis().SetLabelSize(0.12)
    dummy_ratio.GetYaxis().CenterTitle()
    dummy_ratio.GetYaxis().SetNdivisions(505)
    dummy_ratio.GetXaxis().SetTickLength(0.07)

    # dummy_inset = ROOT.TGraph(2)
    # dummy_inset.SetPoint(0, xmin_inset, ymin_inset)
    # dummy_inset.SetPoint(1, xmax_inset, ymax_inset)
    # dummy_inset.SetMarkerSize(0.0)
    # dummy_inset.SetTitle("")
    # dummy_inset.GetXaxis().SetTitle("")
    # dummy_inset.GetYaxis().SetTitle("")
    # dummy_inset.GetXaxis().SetNdivisions(101)
    # dummy_inset.GetXaxis().SetLabelSize(0.0)
    # dummy_inset.GetYaxis().SetNdivisions(503)
    # dummy_inset.GetYaxis().SetLabelSize(0.1)

    ############################################################################
    # Now put together the plot

    canvas = ROOT.TCanvas("canvas", "canvas", 700, 600)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetEndErrorSize(3)
    y_ratio = 0.25
    ratio_seperation = 0.02
    left_margin = 0.1
    right_margin = 0.02

    # Main pad
    pad_main = ROOT.TPad("", "", 0., y_ratio+ratio_seperation, 1., 1.)
    pad_main.SetTopMargin(0.05)
    pad_main.SetBottomMargin(0.0)
    pad_main.SetLeftMargin(left_margin)
    pad_main.SetRightMargin(right_margin)
    #pad_main.SetLogy(True)
    pad_main.SetTickx()
    pad_main.SetTicky()
    pad_main.Draw()

    # Ratio pad
    pad_ratio = ROOT.TPad("", "", 0., 0., 1., y_ratio)
    pad_ratio.SetTopMargin(0.02)
    pad_ratio.SetBottomMargin(0.4)
    pad_ratio.SetLeftMargin(left_margin)
    pad_ratio.SetRightMargin(right_margin)
    pad_ratio.SetTickx()
    pad_ratio.SetTicky()
    pad_ratio.Draw()
    pad_main.Draw()

    # Inset pad
    # x1_inset = 0.535
    # y1_inset = 0.568
    # w_inset = 0.25
    # h_inset = 0.23
    # pad_inset = ROOT.TPad("", "", x1_inset, y1_inset, x1_inset+w_inset, y1_inset+h_inset)
    # left_margin_inset = 0.1
    # right_margin_inset = 0.01
    # pad_inset.SetTopMargin(0.01)
    # pad_inset.SetBottomMargin(0.1)
    # pad_inset.SetLeftMargin(left_margin_inset)
    # pad_inset.SetRightMargin(right_margin_inset)
    # pad_inset.SetTickx()
    # pad_inset.SetTicky()
    # pad_inset.Draw()

    pad_main.cd()
    dummy.Draw("AP")
    dummy.GetXaxis().SetRangeUser(xmin, xmax)
    dummy.GetYaxis().SetRangeUser(ymin, ymax)

    # Make legend for predictions
    #l_prediction = ROOT.TLegend(.12, .70, .48, .82)
    l_prediction = ROOT.TLegend(.34, .09, .7, .14)
    l_prediction.SetTextSize(0.038)
    for p in predictions:
        l_prediction.AddEntry(p.curve, p.legtext, "lf")
    l_prediction.Draw()

    # Additional info for predictions
    t_leg = [
        #makeText("Czakon, Fiedler, Mitov, PRL 110 (2013) 252004",            .13, .67, font=43, size=13, align=12, color=1),
        makeText("m_{top} = 172.5 GeV, #alpha_{s}(M_{Z}) = 0.118 #pm 0.001", .355, .055, font=43, size=16, align=12, color=1),
    ]
    for t in t_leg:
        t.Draw()

    # Make legend for measurements

    #l_measurement = ROOT.TLegend(.6, .17, .99, .42)
    l_measurement = ROOT.TLegend(.1, .5, .7, .8)
    l_measurement.SetTextSize(0.043)

    distance = "#kern[-8.0]{ }" # adding this to the title decreases the distance between marker and text in legend
    for m in measurements:
        l_measurement.AddEntry(m.graph, distance+m.legtext, "p")

    l_measurement.Draw()

    # Draw preliminary label if needed
    if hasPreliminary:
        t_prelim = makeText("* Preliminary", .87, .04, font=43, size=12, align=12, color=13)
        t_prelim.Draw()

    # Draw the prediction curves
    for p in reversed(predictions):
        p.curve.Draw("3 SAME")
        p.curve.Draw("lx SAME")

    # Draw measurement markers
    for m in measurements:
        m.graph.Draw("P SAME")

    # Draw references
    t_refs = []
    x_start = 0.65
    y_start = 0.135
    for i, ref in enumerate(refs):
        x = x_start
        y = y_start - i*0.025
        if i >= 5:
            x += 0.17
            y = y_start - (i-5) * 0.025
        refnr = i+1
        ref_string = "["+str(refnr)+"] "+ref
        t_refs.append( makeText(ref_string, x, y, font=43, size=10, align=14, color=13) )
    #for t in t_refs:
    #    t.Draw()

    # Draw ATLAS+CMS Logo
    titles = []
    x_title_left = 0.14
    x_title_right = x_title_left+0.275
    y_title = 0.86
    y_subtitle = 0.82
    prelimlabel = "Preliminary"
    if internal:
        prelimlabel = "Internal"
    else:
        x_title_right += 0.04

    #titles.append( makeText("ATLAS+CMS",     x_title_left,  y_title, font=62, size=0.045, align=11, color=1) )
    titles.append( makeText("CMS",     x_title_left,  y_title, font=62, size=0.06, align=11, color=1) )
    #titles.append( makeText(prelimlabel,   x_title_right-0.07, y_title, font=51, size=0.06, align=31, color=1) )
    #titles.append( makeText("LHC#it{top}WG", x_title_left,  y_subtitle, font=42, size=0.030, align=11, color=1) )
    #titles.append( makeText("March 2024"   , x_title_right-0.224, y_subtitle, font=42, size=0.030, align=31, color=1) )

    for t in titles:
        t.Draw()

    # Last step, redraw axes
    pad_main.RedrawAxis()

    ################################################################################
    # Draw ratio pad
    pad_ratio.cd()
    dummy_ratio.Draw("AP")
    dummy_ratio.GetXaxis().SetRangeUser(xmin, xmax)
    dummy_ratio.GetYaxis().SetRangeUser(ymin_ratio, ymax_ratio)

    # Draw predictions
    for p in reversed(predictions):
        p.curve_ratio.Draw("3 SAME")
        p.curve_ratio.Draw("lx SAME")

    # Draw measurement markers
    ratios = []
    for m in measurements:
        theoryValue = predictions[i_lhc].curve.Eval(m.sqrts)
        if "Tevatron" in m.name:
            theoryValue = predictions[i_tevatron].curve.Eval(m.sqrts)
        ratio = m.getRatio(theoryValue).Clone(m.name+"_ratio")
        ratios.append(ratio)

    for r in ratios:
        r.Draw("P SAME")

    # Last step, redraw axes
    pad_ratio.RedrawAxis()

    ################################################################################
    # Draw inset
    # pad_inset.cd()
    # dummy_inset.Draw("AP")
    # dummy_inset.GetXaxis().SetRangeUser(xmin_inset, xmax_inset)
    # dummy_inset.GetYaxis().SetRangeUser(ymin_inset, ymax_inset)
    #
    # p_13TeV.graph.Draw("3 SAME")
    # p_13TeV.graph.Draw("lx SAME")
    # p_13p6TeV.graph.Draw("3 SAME")
    # p_13p6TeV.graph.Draw("lx SAME")
    #
    # for m in measurements:
    #     m.graph.Draw("P SAME")
    #
    # # Draw special bin labels on X axis
    # BinLabels = [
    #     (13.0, "13 TeV"),
    #     (13.6, "13.6 TeV"),
    # ]
    # xaxis_labels = []
    # for (value, label) in BinLabels:
    #     x_length_rel = 1.0 - left_margin_inset -right_margin_inset
    #     x_length_tev = xmax_inset - xmin_inset
    #     x_length_perGeV = x_length_rel/x_length_tev
    #     x = left_margin_inset+(value-xmin_inset)*x_length_perGeV
    #     y = 0.08
    #     xaxis_labels.append( makeText(label, x, y, font=42, size=0.08, align=23, color=1) )
    # for t in xaxis_labels:
    #     t.Draw()
    #
    # pad_inset.RedrawAxis()

    ################################################################################
    # Save Canvas
    savename = "tWvsSqrtS_CMS"
    if internal:
        savename = "tWvsSqrtS_LHC_internal"
    canvas.SaveAs(savename+".pdf")
    canvas.SaveAs(savename+".png")




if __name__ == "__main__":
    print("Make plot with 'Internal' label...")
    main(internal=True)
    print("Make plot with 'Preliminary' label...")
    main(internal=False)
