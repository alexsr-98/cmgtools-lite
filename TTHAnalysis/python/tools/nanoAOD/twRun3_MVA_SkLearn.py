from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module 
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection as NanoAODCollection 
from CMGTools.TTHAnalysis.treeReAnalyzer import EventLoop

from CMGTools.TTHAnalysis.treeReAnalyzer import Collection as CMGCollection
from CMGTools.TTHAnalysis.tools.nanoAOD.friendVariableProducerTools import declareOutput, writeOutput

from CMGTools.TTHAnalysis.tools.nanoAOD.TopRun3_modules import ch

import onnxruntime as rt
import numpy as np
import ROOT

class tW_MVA(Module):
    def __init__(self, label = "", MVApath1j1b = "", MVApath2j1b = "", MVApath1j1b_mm = "", MVApath1j1b_ee = "", isMC = True, jecvars = ["jesTotal", "jer"], lepvars = ["mu"]):
        ''' 
        Module to compute the MVA for tW analysis.
        '''
        if (MVApath1j1b == ""):
            raise RuntimeError("No model was provided.")

        #self.vars1j1b = ["nJetSel20{v}_Recl", 
        #                 "Jet1_Pt{v}",
        #                 "JetLoose1_Pt{v}",
        #                 "Lep1Lep2Jet1MET_M{v}",
        #                 "Lep1Lep2Jet1_C{v}",
        #                 "Lep1Lep2Jet1_Pt{v}"]
        
        self.vars1j1b = ["Lep1Jet1_Pt{v}",
                         "Lep1Lep2Jet1MET_M{v}",
                         "Lep1Lep2Jet1_Pt{v}",
                         "Lep1Lep2Jet1MET_PtOverHTtot{v}",
                         "Lep1_PtLep2_PtOverHTtot{v}",
                         "Lep1Jet1_DR{v}",
                         "Lep1Lep2Jet1MET_Mt{v}"]

        # -- Variables to be used in the MVA
        self.vars2j1b = ["Jet2_Pt{v}",
                         "Lep1Jet1_DR{v}",
                         "Lep12Jet12_DR{v}"]

        # -- Variables to be used in the MVA - Multiclass
        self.vars1j1b_jecs = ["JetLoose1_Pt{v}", "Lep1Lep2Jet1_Pt{v}", "Mll", "Lep1Lep2_DPhi", "Lep1Jet1_Pt{v}", "Lep1_Pt", "Jet1_Pt{v}", "Lep1Lep2Jet1_M{v}"]

        self.vars1j1b_leps = ["JetLoose1_Pt{v}", "Lep1Lep2Jet1_Pt{v}", "Mll{v}", "Lep1Lep2_DPhi{v}", "Lep1Jet1_Pt{v}", "Lep1_Pt{v}", "Jet1_Pt{v}", "Lep1Lep2Jet1_M{v}"]

        # -- Variables to be used in the MVA - Multiclass
        self.vars2j1b_jecs = ["Mll", "Lep12Jet12_DR{v}", "Lep1Lep2Jet1_Pt{v}", "Lep1Jet1_DR{v}", "Lep1Lep2_DR", "Lep2_Pt", "Jet2_Pt{v}", "Lep1Lep2Jet1_C{v}"]

        self.vars2j1b_leps = ["Mll{v}", "Lep12Jet12_DR{v}", "Lep1Lep2Jet1_Pt{v}", "Lep1Jet1_DR{v}", "Lep1Lep2_DR{v}", "Lep2_Pt{v}", "Jet2_Pt{v}", "Lep1Lep2Jet1_C{v}"]
        
        # -- Set ONNX options to be single threaded
        session_options = rt.SessionOptions()
        session_options.intra_op_num_threads = 1
        session_options.inter_op_num_threads = 1
        session_options.execution_mode = rt.ExecutionMode.ORT_SEQUENTIAL
 
        # -- Load the model and save it into an attribute. 
        self.MVApath1j1b = MVApath1j1b
        self.model_1j1b = rt.InferenceSession(self.MVApath1j1b, providers = ["CPUExecutionProvider"], sess_options = session_options)
            
        self.MVApath2j1b = MVApath2j1b
        self.model_2j1b = rt.InferenceSession(self.MVApath2j1b, providers = ["CPUExecutionProvider"], sess_options = session_options)

        self.MVApath1j1b_mm = MVApath1j1b_mm
        self.model_1j1b_mm = rt.InferenceSession(self.MVApath1j1b_mm, providers = ["CPUExecutionProvider"], sess_options = session_options)

        self.MVApath1j1b_ee = MVApath1j1b_ee
        self.model_1j1b_ee = rt.InferenceSession(self.MVApath1j1b_ee, providers = ["CPUExecutionProvider"], sess_options = session_options)

        self.input_name_1j1b = self.model_1j1b.get_inputs()[0].name
        self.label_name_1j1b = self.model_1j1b.get_outputs()[0].name
        self.input_name_2j1b = self.model_2j1b.get_inputs()[0].name
        self.label_name_2j1b = self.model_2j1b.get_outputs()[0].name
        self.input_name_1j1b_mm = self.model_1j1b_mm.get_inputs()[0].name
        self.label_name_1j1b_mm = self.model_1j1b_mm.get_outputs()[0].name
        self.input_name_1j1b_ee = self.model_1j1b_ee.get_inputs()[0].name
        self.label_name_1j1b_ee = self.model_1j1b_ee.get_outputs()[0].name
        # -- Nominal variables (MVA outputs)
        self.MVAbranches = ["mvaRF_1j1b", "mvaRF_2j1b", "mvaRF_1j1b_mm", "mvaRF_1j1b_ee"]
        self.branches = []
        
        # -- Variations
        self.label    = "" if (label in ["", None]) else ("_" + label)

        self.systsJEC   = {0: ""} # 0 is nominal
        self.systsLepEn = {} # we don't include the nominal in this case
        self.isMC       = isMC
        if not self.isMC:
            jecvars = []
            if "elscale" in lepvars:
                lepvars = ["elscale"]
            else:
                lepvars = []

        if len(jecvars):
            for i, var in enumerate(jecvars):
                self.systsJEC[i+1]    = "_%sUp"%var
                self.systsJEC[-(i+1)] = "_%sDown"%var

        if len(lepvars):
            for i, var in enumerate(lepvars):
                self.systsLepEn[i+1]    = "_%sUp"%var
                self.systsLepEn[-(i+1)] = "_%sDown"%var
        
        # Extend the branches for systematic variations
        for delta,var in self.systsJEC.items():
            self.branches.extend([br + self.label + var for br in self.MVAbranches])
        for delta,var in self.systsLepEn.items():
            self.branches.extend([br + self.label + var for br in self.MVAbranches])

        return 

    # old interface (CMG)
    def listBranches(self):
        return self.branches[:]
    def __call__(self,event):
        return self.run(event, CMGCollection, "met")

    # new interface (nanoAOD-tools)
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        declareOutput(self, wrappedOutputTree, self.branches)
    def analyze(self, event):
        writeOutput(self, self.run(event, NanoAODCollection))
        return True

    # Common processing
    def run(self, event, Collection):
        # ============================ Definitions and declarations
        allret = {}
        
        # -- Get the object collections used during the training -- #
        #leps = [l for l in Collection(event, "LepGood")]
        #nLepGood = len(leps)

        # ============================ Initialisations
        for var in self.MVAbranches:
            for delta,sys in self.systsJEC.items(): # Remember that the nominal is included in the systsJEC
                allret[var + sys] = -99

        for var in self.MVAbranches:
            for delta,sys in self.systsLepEn.items():
                allret[var + sys] = -99
        
        # ============================ Calculations
        # JECs variations
        if (event.nLepGood >= 2):
            for delta,sys in self.systsJEC.items():
                if event.channel == ch.ElMu:
                    # -- 1j1b -- #
                    if getattr(event, 'nJetSel30{v}_Recl'.format(v = sys if "unclustEn" not in sys else "")) > 0 and getattr(event, 'nBJetSelMedium30{v}_Recl'.format(v = sys if "unclustEn" not in sys else "")) >= 0: # we require 1j1b, the bjet requirement helps to improve the process time
                        # We pass the value of the input variables (nominal and varied) as a numpy array 
                        inputVars1j1b = np.array([getattr(event, var.format(v = sys)) for var in self.vars1j1b_jecs], dtype = np.float32).reshape(1, -1) # we need to reshape the array to have the correct shape
                        allret["mvaRF_1j1b" + sys] = self.model_1j1b.run(None, {self.input_name_1j1b: inputVars1j1b})[1][0][1] # Run returns: (prediction, [class probabilities]) so we take class probabilities and from there the signal probability
                    # -- 2j1b -- #
                    if getattr(event, 'nJetSel30{v}_Recl'.format(v = sys if "unclustEn" not in sys else "")) > 1 and getattr(event, 'nBJetSelMedium30{v}_Recl'.format(v = sys if "unclustEn" not in sys else "")) > 0: # we require 2j1b
                        inputVars2j1b = np.array([getattr(event, var.format(v = sys)) for var in self.vars2j1b_jecs], dtype = np.float32).reshape(1, -1)
                        allret["mvaRF_2j1b" + sys] = self.model_2j1b.run(None, {self.input_name_2j1b: inputVars2j1b})[1][0][1]
                
                if event.channel == ch.Muon:
                    # -- 1j1b mm -- #
                    if getattr(event, 'nJetSel30{v}_Recl'.format(v = sys if "unclustEn" not in sys else "")) > 0 and getattr(event, 'nBJetSelMedium30{v}_Recl'.format(v = sys if "unclustEn" not in sys else "")) > 0: # we require 1j1b, the bjet requirement helps to improve the process time
                        # We pass the value of the input variables (nominal and varied) as a numpy array 
                        inputVars1j1b = np.array([getattr(event, var.format(v = sys)) for var in self.vars1j1b_jecs], dtype = np.float32).reshape(1, -1) # we need to reshape the array to have the correct shape
                        allret["mvaRF_1j1b_mm" + sys] = self.model_1j1b_mm.run(None, {self.input_name_1j1b_mm: inputVars1j1b})[1][0][1] # Run returns: (prediction, [class probabilities]) so we take class probabilities and from there the signal probability
                if event.channel == ch.Elec:
                    # -- 1j1b ee -- #
                    if getattr(event, 'nJetSel30{v}_Recl'.format(v = sys if "unclustEn" not in sys else "")) > 0 and getattr(event, 'nBJetSelMedium30{v}_Recl'.format(v = sys if "unclustEn" not in sys else "")) > 0: # we require 1j1b, the bjet requirement helps to improve the process time
                        # We pass the value of the input variables (nominal and varied) as a numpy array 
                        inputVars1j1b = np.array([getattr(event, var.format(v = sys)) for var in self.vars1j1b_jecs], dtype = np.float32).reshape(1, -1) # we need to reshape the array to have the correct shape
                        allret["mvaRF_1j1b_ee" + sys] = self.model_1j1b_ee.run(None, {self.input_name_1j1b_ee: inputVars1j1b})[1][0][1] # Run returns: (prediction, [class probabilities]) so we take class probabilities and from there the signal probability

        # Lepton energy scale variations
        for delta,sys in self.systsLepEn.items():
            if getattr(event, "nLepGood" + sys[1:]) >= 2:
                if getattr(event, "channel" + sys) == ch.ElMu:
                    # -- 1j1b -- #
                    if getattr(event, 'nJetSel30{v}_Recl'.format(v = sys)) > 0 and getattr(event, 'nBJetSelMedium30{v}_Recl'.format(v = sys)) >= 0:
                        inputVars1j1b = np.array([getattr(event, var.format(v = sys)) for var in self.vars1j1b_leps], dtype = np.float32).reshape(1, -1)
                        allret["mvaRF_1j1b" + sys] = self.model_1j1b.run(None, {self.input_name_1j1b: inputVars1j1b})[1][0][1]
                    # -- 2j1b -- #
                    if getattr(event, 'nJetSel30{v}_Recl'.format(v = sys)) > 1 and getattr(event, 'nBJetSelMedium30{v}_Recl'.format(v = sys)) > 0:
                        inputVars2j1b = np.array([getattr(event, var.format(v = sys)) for var in self.vars2j1b_leps], dtype = np.float32).reshape(1, -1)
                        allret["mvaRF_2j1b" + sys] = self.model_2j1b.run(None, {self.input_name_2j1b: inputVars2j1b})[1][0][1]
                
                if getattr(event, "channel" + sys) == ch.Muon:
                    # -- 1j1b mm -- #
                    if getattr(event, 'nJetSel30{v}_Recl'.format(v = sys)) > 0 and getattr(event, 'nBJetSelMedium30{v}_Recl'.format(v = sys)) > 0:
                        inputVars1j1b = np.array([getattr(event, var.format(v = sys)) for var in self.vars1j1b_leps], dtype = np.float32).reshape(1, -1)
                        allret["mvaRF_1j1b_mm" + sys] = self.model_1j1b_mm.run(None, {self.input_name_1j1b_mm: inputVars1j1b})[1][0][1]
                if getattr(event, "channel" + sys) == ch.Elec:
                    # -- 1j1b ee -- #
                    if getattr(event, 'nJetSel30{v}_Recl'.format(v = sys)) > 0 and getattr(event, 'nBJetSelMedium30{v}_Recl'.format(v = sys)) > 0:
                        inputVars1j1b = np.array([getattr(event, var.format(v = sys)) for var in self.vars1j1b_leps], dtype = np.float32).reshape(1, -1)
                        allret["mvaRF_1j1b_ee" + sys] = self.model_1j1b_ee.run(None, {self.input_name_1j1b_ee: inputVars1j1b})[1][0][1]        
        # ============================ Return
        return allret

        
