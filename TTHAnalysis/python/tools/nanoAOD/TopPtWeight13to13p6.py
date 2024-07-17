from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import ROOT as r
import os


class TopPtWeight13to13p6(Module):
    # =================== ### Main methods
    def __init__(self):
        return


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch('TopPtWeight13to13p6', 'F')
        return


    def analyze(self, event):
        topptw = 1.0
        if event.isTop and event.Top2_pt >= 0:
            # From a fit to the top pt distribution in MC 13TeV vs 13p6 TeV
            topptw *= 0.99056 + 7.496e-5 * event.Top1_pt
            topptw *= 0.99056 + 7.496e-5 * event.Top2_pt
            topptw = r.TMath.Sqrt(topptw)

        self.out.fillBranch('TopPtWeight13to13p6', topptw)
        return True
