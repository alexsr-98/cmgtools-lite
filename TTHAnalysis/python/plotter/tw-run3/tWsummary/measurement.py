import ROOT

class Measurement:
    def __init__(self, name, sqrts, lumi, xs, uncertUp, uncertDown, markerstyle, markersize, color, xoffset, refnr, ref):
        self.name = name
        self.ref = ref
        self.lumi = lumi
        self.xs = xs
        self.sqrts = sqrts
        self.xoffset = xoffset
        self.markersize = markersize
        self.markerstyle = markerstyle
        self.color = color
        self.graph = self.getGraph(sqrts, xs, uncertDown, uncertUp, xoffset)
        self.graph.SetTitle("")
        self.graph.SetMarkerStyle(markerstyle)
        self.graph.SetMarkerSize(1.7*markersize)
        self.graph.SetMarkerColor(color)
        self.graph.SetLineColor(color)
        sqrts_string = self.convertToSqrtSString(sqrts)
        #self.legtext = name+" ("+sqrts_string+" TeV, "+lumi+" fb^{-1}_{}) ["+str(refnr)+"]"
        if ref == "":
            self.legtext = name+" ("+sqrts_string+" TeV, "+lumi+" fb^{-1}_{})"
        else:
            self.legtext = name+" ("+sqrts_string+" TeV, "+lumi+" fb^{-1}_{}), "+str(ref)


    def getGraph(self, sqrts, xs, uncertDown, uncertUp, xoffset):
        g = ROOT.TGraphAsymmErrors(1)
        g.SetPoint(0, sqrts+3*xoffset, xs)
        g.SetPointError(0, 0.0, 0.0, uncertDown, uncertUp)
        return g

    def convertToSqrtSString(self, sqrts):
        if sqrts in [13.0, 7.0, 8.0]:
            return "%i" %sqrts
        elif sqrts in [13.6]:
            return "%.1f" %sqrts
        else:
            return "%.2f" %sqrts

    def getRatio(self, theoryValue):
        ratio = ROOT.TGraphAsymmErrors(1)
        errYUp = self.graph.GetErrorYhigh(0)/theoryValue
        errYDown = self.graph.GetErrorYlow(0)/theoryValue
        ratio.SetPoint(0, self.sqrts+3*self.xoffset, self.xs/theoryValue)
        ratio.SetPointError(0, 0.0, 0.0, errYDown, errYUp)
        ratio.SetTitle("")
        ratio.SetMarkerStyle(self.markerstyle)
        ratio.SetMarkerSize(1.4*self.markersize)
        ratio.SetMarkerColor(self.color)
        ratio.SetLineColor(self.color)
        return ratio
