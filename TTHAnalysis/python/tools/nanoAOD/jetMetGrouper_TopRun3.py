import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
# group jec components
# taken from https://docs.google.com/spreadsheets/d/1Feuj1n0MdotcPq19Mht7SUIgvkXkA4hiB0BxEuBShLw/edit?ouid=111820255692530107608&usp=sheets_home&ths=true

import os
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import math

class jetmetGrouper_TopRun3(Module):
    def __init__(self, groups, jetbranch, metbranch, dumpMore = []):
        self.groups    = groups
        self.jetbranch = jetbranch
        self.metbranch = metbranch
        self.dumpMore  = dumpMore
        return

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.wrappedOutputTree = wrappedOutputTree
        for group in self.groups:
            for sign in ['Up','Down']:
                self.wrappedOutputTree.branch("%s_pt_jes%s%s"%(  self.jetbranch,group,sign), "F", lenVar="n%s"%self.jetbranch)
                self.wrappedOutputTree.branch("%s_mass_jes%s%s"%(self.jetbranch,group,sign), "F", lenVar="n%s"%self.jetbranch)
                self.wrappedOutputTree.branch("%s_pt_jes%s%s"%(  self.metbranch,group,sign), "F")
                self.wrappedOutputTree.branch("%s_phi_jes%s%s"%( self.metbranch,group,sign), "F")
        for br in self.dumpMore:
            self.wrappedOutputTree.branch(**br)
        return

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass


    def analyze(self, event):
        """Nominal variables"""
        jets   = Collection(event, '%s'%self.jetbranch)
        met    = Object(event,     self.metbranch)
        metPt  = getattr(met, "pt",  "pt")
        metPhi = getattr(met, "phi", "phi")
        metPx  = metPt*math.cos(metPhi)
        metPy  = metPt*math.sin(metPhi)
        for group in self.groups:
            for sign in ['Up','Down']:
                jetVarPt = []
                jetVarM  = []
                metPxVar = metPx
                metPyVar = metPy
                for j in jets:
                    thePt   = getattr(j, "pt_nom",   "pt")
                    theMass = getattr(j, "mass_nom", "mass")
                    thePhi  = getattr(j, "phi")
                    jetVarPt.append(0.)
                    jetVarM.append(0.)
                    for comp in self.groups[group]:
                        jetVarPt[-1] = (jetVarPt[-1]**2 + (getattr(j, "pt_jes"   + comp + sign) - thePt)**2)**0.5
                        jetVarM[-1]  = (jetVarM[-1]**2  + (getattr(j, "mass_jes" + comp + sign) - theMass)**2)**0.5
                    
                    # To split in flavour (option 1)
                    if group == "FlavorQCD_bottom" and j.hadronFlavour != 5:
                        jetVarPt[-1] = 0
                        jetVarM[-1]  = 0
                    if group == "FlavorQCD_charm" and j.hadronFlavour != 4:
                        jetVarPt[-1] = 0
                        jetVarM[-1]  = 0
                    if group == "FlavorQCD_light" and j.hadronFlavour != 0:
                        jetVarPt[-1] = 0
                        jetVarM[-1]  = 0
                    
                    # To split in flavour (option 2)
                    if group == "FlavorPure_Bottom" and abs(j.partonFlavour) != 5:
                        jetVarPt[-1] = 0
                        jetVarM[-1]  = 0
                    if group == "FlavorPure_Charm" and abs(j.partonFlavour) != 4:
                        jetVarPt[-1] = 0
                        jetVarM[-1]  = 0
                    if group == "FlavorPure_Quark" and abs(j.partonFlavour) != 1 and abs(j.partonFlavour) != 2 and abs(j.partonFlavour) != 3:
                        jetVarPt[-1] = 0
                        jetVarM[-1]  = 0
                    if group == "FlavorPure_Gluon" and abs(j.partonFlavour) != 21:
                        jetVarPt[-1] = 0
                        jetVarM[-1]  = 0
                    
                    metPxVar     = metPxVar - (1 if sign == 'Up' else -1) * jetVarPt[-1] * math.cos(thePhi)
                    metPyVar     = metPyVar - (1 if sign == 'Up' else -1) * jetVarPt[-1] * math.sin(thePhi)
                    jetVarPt[-1] = thePt    + (1 if sign == 'Up' else -1) * jetVarPt[-1]
                    jetVarM[-1]  = theMass  + (1 if sign == 'Up' else -1) * jetVarM[-1]
                    
                self.wrappedOutputTree.fillBranch("%s_pt_jes%s%s"%(  self.jetbranch, group, sign), jetVarPt)
                self.wrappedOutputTree.fillBranch("%s_mass_jes%s%s"%(self.jetbranch, group, sign), jetVarM)
                self.wrappedOutputTree.fillBranch("%s_pt_jes%s%s"%(  self.metbranch, group, sign), (metPxVar**2 + metPyVar**2)**0.5)
                self.wrappedOutputTree.fillBranch("%s_phi_jes%s%s"%( self.metbranch, group, sign), math.atan2(metPyVar,  metPxVar))
        
        for br in self.dumpMore:
            if "lenVar" in br:
                tmplist = []
                for i in range(getattr(event, br["lenVar"])):
                    tmplist.append(getattr(event, br['name'])[i])
                self.wrappedOutputTree.fillBranch(br['name'], tmplist)
            else:
                self.wrappedOutputTree.fillBranch(br['name'], getattr(event, br['name']))
        return True


moreVars = [
    {'name' : 'Jet_pt_nom'           ,      'rootBranchType' : 'F','lenVar': 'nJet'},
    {'name' : 'Jet_pt_jerUp'         ,      'rootBranchType' : 'F','lenVar': 'nJet'},
    {'name' : 'Jet_pt_jerDown'       ,      'rootBranchType' : 'F','lenVar': 'nJet'},
    {'name' : 'Jet_pt_jesTotalUp',          'rootBranchType' : 'F','lenVar': 'nJet'},
    {'name' : 'Jet_pt_jesTotalDown',        'rootBranchType' : 'F','lenVar': 'nJet'},
    {'name' : 'Jet_mass_nom'         ,      'rootBranchType' : 'F','lenVar': 'nJet'},
    {'name' : 'Jet_mass_jerUp'       ,      'rootBranchType' : 'F','lenVar': 'nJet'},
    {'name' : 'Jet_mass_jerDown'     ,      'rootBranchType' : 'F','lenVar': 'nJet'},
    {'name' : 'Jet_mass_jesTotalUp',        'rootBranchType' : 'F','lenVar': 'nJet'},
    {'name' : 'Jet_mass_jesTotalDown',      'rootBranchType' : 'F','lenVar': 'nJet'},
    {'name' : 'PuppiMET_T1_pt'               ,   'rootBranchType' : 'F'},
    {'name' : 'PuppiMET_T1_pt_jesTotalUp',       'rootBranchType' : 'F'},
    {'name' : 'PuppiMET_T1_pt_jesTotalDown',     'rootBranchType' : 'F'},
    {'name' : 'PuppiMET_T1_pt_unclustEnUp'   ,   'rootBranchType' : 'F'},
    {'name' : 'PuppiMET_T1_pt_unclustEnDown' ,   'rootBranchType' : 'F'},
    {'name' : 'PuppiMET_T1_phi'              ,   'rootBranchType' : 'F'},
    {'name' : 'PuppiMET_T1_phi_jesTotalUp',      'rootBranchType' : 'F'},
    {'name' : 'PuppiMET_T1_phi_jesTotalDown',    'rootBranchType' : 'F'},
    {'name' : 'PuppiMET_T1_phi_unclustEnUp',     'rootBranchType' : 'F'},
    {'name' : 'PuppiMET_T1_phi_unclustEnDown',   'rootBranchType' : 'F'},
]

###for jer in list(range(6)):
###    for var in ['Up','Down']:
###        moreVars.extend( [
###            {'name' : 'Jet_pt_jer%s%s'%(jer,var) ,   'rootBranchType' : 'F','lenVar': 'nJet'},
###            {'name' : 'Jet_mass_jer%s%s'%(jer,var) , 'rootBranchType' : 'F','lenVar': 'nJet'},
###            {'name' : 'MET_T1_pt_jer%s%s'%(jer,var) ,'rootBranchType' : 'F'},
###            {'name' : 'MET_T1_phi_jer%s%s'%(jer,var),'rootBranchType' : 'F'},
###        ])

### Aqui hay cosas que se pueden quitar, no son necesarias, pero bueno las dejo por no tocar mucho el asunto.
groups = {'HF'                 : ['PileUpPtHF', 'RelativeJERHF', 'RelativePtHF'],
          'BBEC1_year'         : ['RelativeJEREC1', 'RelativePtEC1', 'RelativeStatEC'],
          'FlavorQCD'          : ['FlavorQCD'],
          'FlavorQCD_bottom'   : ['FlavorQCD'],
          'FlavorQCD_charm'    : ['FlavorQCD'],
          'FlavorQCD_light'    : ['FlavorQCD'],
          #'FlavorPure_Quark'    : ['FlavorPureQuark'],
          #'FlavorPure_Charm'    : ['FlavorPureCharm'],
          #'FlavorPure_Bottom'   : ['FlavorPureBottom'],
          #'FlavorPure_Gluon'    : ['FlavorPureGluon'],
          'RelativeSample_year': ['RelativeSample'],
          'EC2'                : ['PileUpPtEC2'],
          'HF_year'            : ['RelativeStatHF'],
          'RelativeBal'        : ['RelativeBal'],
          'Absolute_year'      : ['AbsoluteStat', 'RelativeStatFSR', 'TimePtEta'],
          'BBEC1'              : ['PileUpPtBB', 'PileUpPtEC1', 'RelativePtBB'],
          'EC2_year'           : ['RelativeJEREC2', 'RelativePtEC2'],
          'Absolute'           : ['AbsoluteMPFBias', 'AbsoluteScale', 'Fragmentation', 'PileUpDataMC',
                                  'PileUpPtRef', 'RelativeFSR', 'SinglePionECAL', 'SinglePionHCAL'],
}

jetMetCorrelate2022 = lambda  : jetmetGrouper_TopRun3(groups, "Jet", "PuppiMET_T1", dumpMore = moreVars)
