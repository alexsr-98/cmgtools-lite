import ROOT as r
import sys, os
from array import array

r.gROOT.SetBatch(True)
print("===== Setting binning of the BDT procedure\n")
print("> Preliminaries...")
printthings = 1



friendspath = "/pool/phedexrw/userstorage/vrbouza/proyectos/tw_run2/productions"
mcpath      = "/pool/ciencias/nanoAODv6/29jan2020_MC"
prodname    = "2020-09-20"
friendfolders = ["0_yeartag", "1_lepmerge_roch", "2_cleaning", "3_varstrigger", "4_scalefactors", "5_mvas"]

#### 2016
year        = 2016
ttbarname   = "Tree_TT_TuneCUETP8M2T4_PSweights_"
xsec        = 831.76
nttbarfiles = 1
#### 2017


#### 2018



#print "> Opening files and creating chains\n"
#thetree = r.TChain("Events", "thetree")
#thechainfriendlist = []
#for friendfolder in friendfolders:
    #thechainfriendlist.append(r.TChain("Friends", friendfolder[2:]))


#def getsumw(fpath):
    #f = r.TFile.Open(fpath, "READ")
    #val = 0
    #for ev in f.Runs:
        #val = ev.genEventSumw_
    #f.Close(); del f
    #return val

#gensumwlist = []
#for ifil in range(nttbarfiles):
    #thetree.Add(mcpath + "/" + str(year) + "/" + ttbarname + str(ifil) + ".root")

    ##print getsumw(mcpath + "/" + str(year) + "/" + ttbarname + str(ifil) + ".root")
    #gensumwlist.append(getsumw(mcpath + "/" + str(year) + "/" + ttbarname + str(ifil) + ".root"))

    #for idx,fchain in enumerate(thechainfriendlist):
        #fchain.Add(friendspath + "/" + prodname + "/" + str(year) + "/" + friendfolders[idx] + "/" + ttbarname + str(ifil) + "_Friend.root")


#for fchain in thechainfriendlist:
    #thetree.AddFriend(fchain)

#thew = xsec/sum(gensumwlist)
#cut1j1b =  "(genWeight * MuonSF * ElecSF * TrigSF * puWeight * bTagWeight * PrefireWeight * {w})".format(w = thew)
#cut1j1b += "* (pass_trigger * Flag_goodVertices * Flag_globalSuperTightHalo2016Filter * Flag_HBHENoiseFilter * Flag_HBHENoiseIsoFilter * Flag_EcalDeadCellTriggerPrimitiveFilter * Flag_BadPFMuonFilter)"
#cut1j1b += "* ((channel == 1) && (nLepGood >= 2) && (LepGood_pt_corrAll[0] > 25) && (LepGood_pt_corrAll[1] > 20) && isSS == 0 && (minMllAFAS_Recl > 20) )"
#cut2j1b = cut1j1b + "* ((nJetSel30_Recl == 2) && (nBJetSelMedium30_Recl == 1))"
#cut1j1b += "* ((nJetSel30_Recl == 1) && (nBJetSelMedium30_Recl == 1))"


#print cut1j1b, cut2j1b
#sys.exit()

#def getBinningForThatVariable(tree, var, thecut, nq):
    #xq          = array('d', [0]*nq)
    #yq          = array('d', [0]*nq)
    #Base  ='''
    ##ifndef ADDSTUFF_H
    ##define ADDSTUFF_H
    ##include<iostream>
    #using namespace std;
    #'''

    #for i in range(0,nq):
        #xq[i] = float(i+1)/nq

    #count = 0
    #print "\n> Constructing C++ file of the binning"
    #c       = r.TCanvas()
    #count   = count + 1
    #tree.Draw(var, thecut)
    #c.Update()
    #print "\n> Obtaining quantiles..."
    #tree.GetHistogram().GetQuantiles(nq, yq, xq)
    #del c
    #Base = Base + '''\n /////////////////////// \n
#Float_t theBDT(Double_t BDT){
  #if      (BDT < %f   ) return 1;
#'''%(yq[0])

    #if (printthings == 1):
        #print "Bin 1:", yq[0]
        #raw_input('Copy this string into your .C function and press enter.')
    #subBin = 2
    #for i in range(1, nq - 1):
        #Base = Base + '  else if (BDT < %f) return %d;\n'%(yq[i],subBin)
        #subBin = subBin + 1
        #if (printthings == 1):
            #print "Bin {i}:".format(i = i + 1), yq[i]
            #raw_input('Press enter.')

    #Base = Base + '  else                       return %d;\n} \n'%subBin
    #Base = Base + '#endif'

#print "> Saving file"
#if not os.path.isdir('temp/{var}_'.format(var = varName)):
    #os.system('mkdir -p temp/{var}_'.format(var = varName))
#output = open('./temp/{var}_/'.format(var = varName) + varName + '.C', 'w')
#output.write(Base)
#output.close()

#print "> Doing sanity checks..."
#r.gROOT.LoadMacro('./temp/{var}_/'.format(var = varName) + varName + '.C+')

#for k in range(1, len(bins)):
    #if not hasattr(r, 'theBDt_bin%d'%k): print 'Something went wrong'
    #else: print '    - Bin', k, 'OK'


#print "> Sanity check correct\n"
#print "> BDT binning done!"

#getBinningForThatVariable(thetree, "SergioBDT_1j1b", cut1j1b, 10)
#getBinningForThatVariable(thetree, "SergioBDT_2j1b", cut2j1b, 6)


def getBinningForThatVariable(fpath, histosname, nq):
    xq          = array('d', [0]*nq)
    yq          = array('d', [0]*nq)
    Base  ='''
    #ifndef ADDSTUFF_H
    #define ADDSTUFF_H
    #include<iostream>
    using namespace std;
    '''

    for i in range(0,nq):
        xq[i] = float(i+1)/nq

    count = 0
    thefile = r.TFile.Open(fpath, "READ")
    theHistos = histosname.split("+")
    thehisto = thefile.Get(theHistos[0])
    if len(theHistos) > 1:
        for histoname in theHistos[1:]:
            thehisto.Add(thefile.Get(histoname))
    
    thehisto.GetQuantiles(nq, yq, xq)
    thefile.Close(); del thefile
    Base = Base + '''\n /////////////////////// \n
Float_t theBDT(Double_t BDT){
  if      (BDT < %f   ) return 1;
'''%(yq[0])


    print("> Histoname", histosname, "with", nq, "bins")
    if (printthings == 1):
        print("Bin 1:", yq[0])
        input('Copy this string into your .C function and press enter.')
    subBin = 2
    for i in range(1, nq - 1):
        Base = Base + '  else if (BDT < %f) return %d;\n'%(yq[i],subBin)
        subBin = subBin + 1
        if (printthings == 1):
            print("Bin {i}:".format(i = i + 1), yq[i])
            input('Press enter.')

    Base = Base + '  else                       return %d;\n} \n'%subBin
    Base = Base + '#endif'
    print(Base)
    return

#########To rebin with more processes add name and + symbol in the second argument of the function, example: "tmvaBDT_1j1b_smallb_ttbar+tmvaBDT_1j1b_smallb_dy"
#getBinningForThatVariable("/nfs/fanae/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/temp_Run3_plots/2023_02_01_Test_RF/2022/1j1t/MVAtrain/plots-tw-1j1t_MVAtrain.root", "tmvaBDT_1j1b_smallb_ttbar", 10)
#getBinningForThatVariable("/nfs/fanae/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/temp_Run3_plots/2023_02_01_Test_RF/2022/2j1t/MVAtrain/plots-tw-2j1t_MVAtrain.root", "tmvaBDT_2j1b_smallb_ttbar", 6)

# For the RF
getBinningForThatVariable("/nfs/fanae/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/temp_Run3_plots/2023_09_30_smallb/2022PostEE/1j1t/MVAtrain/plots-tw-1j1t_MVAtrain.root", "mvaRF_1j1b_smallb_ttbar", 20)
getBinningForThatVariable("/nfs/fanae/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/temp_Run3_plots/2023_09_30_smallb/2022PostEE/2j1t/MVAtrain/plots-tw-2j1t_MVAtrain.root", "mvaRF_2j1b_smallb_ttbar", 12)
getBinningForThatVariable("/nfs/fanae/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/temp_Run3_plots/2023_09_26_smallB/2022PostEE/1j0t/MVAtrain/plots-tw-1j1t_MVAtrain.root", "mvaRF_1j0b_smallb_ttbar", 10)
# RF for mm and ee
getBinningForThatVariable("/nfs/fanae/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/temp_Run3_plots/2023_09_30_smallb/2022PostEE/1j1t-mm/MVAtrain/plots-tw-1j1t_MVAtrain.root", "mvaRF_1j1b_mm_smallb_ttbar", 10)
getBinningForThatVariable("/nfs/fanae/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/temp_Run3_plots/2023_09_30_smallb/2022PostEE/1j1t-ee/MVAtrain/plots-tw-1j1t_MVAtrain.root", "mvaRF_1j1b_ee_smallb_ttbar", 10)




