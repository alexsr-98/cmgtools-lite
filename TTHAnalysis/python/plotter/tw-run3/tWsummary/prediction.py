import ROOT
from math import sqrt, pow

class Prediction:
    def __init__(self, filename, legname, color, alpha=1.0, range=None, range_ratio=None):
        self.range = range
        self.range_ratio = range_ratio
        self.curve, self.curve_ratio = self.get_curve(filename)
        self.curve.SetLineColor(1)
        self.curve.SetLineWidth(1)
        self.curve.SetTitle("")
        self.curve.SetFillColor(color)
        self.curve.SetFillColorAlpha(color, alpha)
        self.curve_ratio.SetLineColor(1)
        self.curve_ratio.SetLineWidth(1)
        self.curve_ratio.SetTitle("")
        self.curve_ratio.SetFillColorAlpha(color, alpha)
        self.legtext = legname


    def get_curve(self, filename):
        values = {
            "central": {},
            "scalesUp": {},
            "scalesDown": {},
        }

        f = open(filename, "r")
        for line in f:
            sqrts = float(line.split()[0]) * 0.001 # convert to TeV
            mt = float(line.split()[1])
            muR = float(line.split()[2])
            muF = float(line.split()[3])
            alphaS = float(line.split()[4])
            xs = float(line.split()[5])
            err = float(line.split()[6])
            pdfUp = float(line.split()[7]) - xs
            pdfDown = xs - float(line.split()[8])

            # convention is that "Up" refers to direction of XS
            # thus, the "Up" error come from setting the scales to a smaller value
            if muR == mt and muF == mt:
                values["central"][sqrts] = (xs, pdfUp, pdfDown)
            elif muR > mt and muF > mt:
                values["scalesDown"][sqrts] = xs
            elif muR < mt and muF < mt:
                values["scalesUp"][sqrts] = xs
            else:
                raise Exception("Unexpected format of prediction file")


        # Go through values again and remove those out of range
        valuesForCurve = {
            "central": {},
            "scalesUp": {},
            "scalesDown": {},
        }
        valuesForRatio = {
            "central": {},
            "scalesUp": {},
            "scalesDown": {},
        }
        if self.range is not None:
            sqrts_min, sqrts_max = self.range
        else:
            sqrts_min, sqrts_max = -1, 9999999.

        if self.range_ratio is not None:
            sqrts_min_ratio, sqrts_max_ratio = self.range_ratio
        else:
            sqrts_min_ratio, sqrts_max_ratio = -1, 9999999.

        for sqrts in sorted(values["central"].keys()):
            if sqrts > sqrts_min and sqrts < sqrts_max:
                valuesForCurve["central"][sqrts] = values["central"][sqrts]
                valuesForCurve["scalesUp"][sqrts] = values["scalesUp"][sqrts]
                valuesForCurve["scalesDown"][sqrts] = values["scalesDown"][sqrts]
            if sqrts > sqrts_min_ratio and sqrts < sqrts_max_ratio:
                valuesForRatio["central"][sqrts] = values["central"][sqrts]
                valuesForRatio["scalesUp"][sqrts] = values["scalesUp"][sqrts]
                valuesForRatio["scalesDown"][sqrts] = values["scalesDown"][sqrts]

        # Fill Curve
        curve = ROOT.TGraphAsymmErrors(len(valuesForCurve["central"]))
        for i, sqrts in enumerate(sorted(valuesForCurve["central"].keys())):
            (central, pdfUp, pdfDown) = valuesForCurve["central"][sqrts]
            scaleErrUp = valuesForCurve["scalesUp"][sqrts] - central
            scaleErrDown = central - valuesForCurve["scalesDown"][sqrts]

            totUp = sqrt(pow(scaleErrUp,2)+pow(pdfUp,2))
            totDown = sqrt(pow(scaleErrDown,2)+pow(pdfDown,2))

            curve.SetPoint(i, sqrts, central)
            curve.SetPointError(i, 0., 0., totDown, totUp)

        # Fill Ratio
        curve_ratio = ROOT.TGraphAsymmErrors(len(valuesForRatio["central"]))
        for i, sqrts in enumerate(sorted(valuesForRatio["central"].keys())):
            (central, pdfUp, pdfDown) = valuesForRatio["central"][sqrts]
            scaleErrUp = valuesForRatio["scalesUp"][sqrts] - central
            scaleErrDown = central - valuesForRatio["scalesDown"][sqrts]

            totUp = sqrt(pow(scaleErrUp,2)+pow(pdfUp,2))
            totDown = sqrt(pow(scaleErrDown,2)+pow(pdfDown,2))

            curve_ratio.SetPoint(i, sqrts, central/central)
            curve_ratio.SetPointError(i, 0., 0., totDown/central, totUp/central)

        return curve,curve_ratio

class PredictionValue:
    def __init__(self, name, xs, scaleUp, scaleDown, pdfUp, pdfDown, alphasUp, alphasDown, xmin, xmax, color):
        self.name = name
        self.graph = self.getGraph(xs, scaleUp, scaleDown, pdfUp, pdfDown, alphasUp, alphasDown, xmin, xmax)
        self.graph.SetLineColor(1)
        self.graph.SetLineWidth(1)
        self.graph.SetTitle("")
        self.graph.SetFillColor(color)

    def getGraph(self, xs, scaleUp, scaleDown, pdfUp, pdfDown, alphasUp, alphasDown, xmin, xmax):
        graph = ROOT.TGraphAsymmErrors(2)
        totUp = sqrt(pow(scaleUp,2)+pow(pdfUp,2)+pow(alphasUp,2))
        totDown = sqrt(pow(scaleDown,2)+pow(pdfDown,2)+pow(alphasDown,2))
        graph.SetPoint(0, xmin, xs)
        graph.SetPointError(0, 0., 0., totDown, totUp)
        graph.SetPoint(1, xmax, xs)
        graph.SetPointError(1, 0., 0., totDown, totUp)
        return graph
