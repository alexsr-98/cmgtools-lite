from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection as Collection
from CMGTools.TTHAnalysis.tools.nanoAOD.friendVariableProducerTools import declareOutput, writeOutput
import os
from copy import deepcopy
import ROOT

class puWeighter( Module ):
    def __init__(self, filename):
        # Open the rootfile and load the weights 
        self.filename = filename
        self.load_histograms()
        self.nvtxVar = "Pileup_nTrueInt"
        self.branches = [("puWeight", "F"), ("puWeight_up", "F"),("puWeight_dn", "F")]
        return

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        declareOutput(self, wrappedOutputTree, self.branches)
        return

    def load_histograms(self):
        self.weights = {}
        histname = "puWeights"
        tf = ROOT.TFile.Open( self.filename )
        self.weights["nom"] = deepcopy(tf.Get(histname))
        self.weights["up"]  = deepcopy(tf.Get(histname + "Up"))
        self.weights["dn"]  = deepcopy(tf.Get(histname + "Down"))
        tf.Close()
        return

    def analyze(self, event):
        """ Method to apply the PU weights """
        ret = {} 
        # Get number of vertices
        nvtx = int(getattr(event, self.nvtxVar)) 
       
        # Fetch weight 
        nombin = self.weights["nom"].FindBin(nvtx) 
        upbin = self.weights["up"].FindBin(nvtx) 
        dnbin = self.weights["dn"].FindBin(nvtx) 

        nomweight = self.weights["nom"].GetBinContent(nombin) 
        upweight = self.weights["up"].GetBinContent(upbin) 
        dnweight = self.weights["dn"].GetBinContent(dnbin) 
       
        # Save in dictionary
        ret["puWeight"] = nomweight
        ret["puWeight_up"] = upweight
        ret["puWeight_dn"] = dnweight

        # Write
        writeOutput(self, ret)
        return True

autoPuWeight = lambda : autoPuWeightModule(puAutoWeight_2016, puAutoWeight_2017, puAutoWeight_2018)
