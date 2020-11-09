systsGroup = {
    #### Statistical
    'mc_stat' : [
        "prop_binch1_bin0",
        "prop_binch1_bin1",
        "prop_binch1_bin2",
        "prop_binch1_bin3",
        "prop_binch1_bin4",
        "prop_binch1_bin5",
        "prop_binch1_bin6",
        "prop_binch1_bin7",
        "prop_binch1_bin8",
        "prop_binch1_bin9",
        "prop_binch2_bin0",
        "prop_binch2_bin1",
        "prop_binch2_bin2",
        "prop_binch2_bin3",
        "prop_binch2_bin4",
        "prop_binch2_bin5",
        "prop_binch3_bin0",
        "prop_binch3_bin1",
        "prop_binch3_bin2",
        "prop_binch3_bin3",
        "prop_binch3_bin4",
        "prop_binch3_bin5",
        "prop_binch3_bin6",
        "prop_binch3_bin7",
    ],

    #### Systematic
    # Experimental
    'jecs': [
        "jer_2016",
        "jer_2017",
        "jer_2018",
        "jes",
    ],
    'trigger': [
        "triggereff_2016",
        "triggereff_2017",
        "triggereff_2018",
    ],
    'pileup': [
        "pileup",
    ],
    'lep': [
        "elecidsf",
        "elecrecosf",
        "muonen_2016",
        "muonen_2017",
        "muonen_2018",
        "muonidsf_stat_2016",
        "muonidsf_stat_2017",
        "muonidsf_stat_2018",
        "muonidsf_syst",
        "muonisosf_stat_2016",
        "muonisosf_stat_2017",
        "muonisosf_stat_2018",
        "muonisosf_syst",
    ],
    'btag': [
        "btagging",
        "mistagging",
    ],
    'lumi': [
        "lumi_2016",
        "lumi_2017",
        "lumi_2018",
        "lumi_BBD",
        "lumi_BCC",
        "lumi_DB",
        "lumi_GS",
        "lumi_LS",
        "lumi_XY",
    ],
    'prefiring' : [
        "prefiring_2016",
        "prefiring_2017",
    ],

    # Normalisation
    'norm' : [
        "ttbar_norm",
        "nonworz_norm",
        "dy_norm",
        "vvttv_norm",
    ],

    # Modelling
    'pdf' : [
        "pdf",
    ],
    'matching' : [
        "ttbar_matching",
    ],
    #'scales' : [
    #],
    'ps' : [
        "fsr_ttbar",
        "isr_ttbar",
    ],
    'colour' : [
        "colour_rec",
    ],
    'ue' : [
        "ue",
    ],
    'toppt' : [
        "topptrew",
    ],
    #'mtop' : [
    #],
    #'ds' : [
    #],
}

dowhat = 'step1'
POIs   = ["r"]

thecard = "./temp_2020_10_30/cards/combinada.root"

groupList   = ['mc_stat', 'jecs', 'trigger', 'pileup', 'lep', 'btag', 'lumi', 'prefiring', 'norm', "pdf", "matching", "ps", "colour", "ue"]

basecommand = '\ncombineTool.py -M MultiDimFit --algo grid --points 100 --rMin 0 --rMax 3 --floatOtherPOIs=1 -m 125  --split-points 1 --setParameters r=1 -t -1 --expectSignal=1 --job-mode SGE --saveInactivePOI 1 '

if dowhat == 'step1':
    for poi in POIs:
        cumulative      = [x for x in ["r"] if x != poi]

        print basecommand + '-n nominal_%s %s --task-name nominal_%s -P %s %s'%(poi, thecard, poi, poi, ",".join(cumulative))

        print basecommand.replace('--algo grid','--algo none').replace("--points 100","").replace("--job-mode SGE","")+ '-n bestfit_%s --saveWorkspace tw-run2_WS.root -P %s '%(poi,poi)

        for group in groupList:
            cumulative += systsGroup[group]
            print basecommand + ' -P %s '%poi + '-n ' + group + '_%s'%poi + ' higgsCombinebestfit_%s.MultiDimFit.mH125.root --snapshotName MultiDimFit  --freezeParameters %s'%(poi,",".join(cumulative)) + ' --task-name %s_%s'%(group,poi)

if dowhat == 'step2':
    for poi in POIs:
        fileList        = [ ]
        print 'hadd higgsCombinenominal_%s.MultiDimFit.mH125.root higgsCombinenominal_%s.POINTS.*.MultiDimFit.mH125.root'%(poi,poi)
        for gr in groupList:
            print 'hadd higgsCombine%s.MultiDimFit.mH125.root higgsCombine%s.POINTS.*.MultiDimFit.mH125.root'%(gr+'_'+poi,gr+'_'+poi)
            fileList.append( "'higgsCombine%s.MultiDimFit.mH125.root:Freeze += %s:%d'"%(gr+'_'+poi,gr,groupList.index(gr)))

        print 'plot1DScan.py higgsCombinenominal_%s.MultiDimFit.mH125.root --others '%poi +' '.join( fileList ) +' --breakdown '  + ','.join(groupList) +',stat'  + " --POI %s "%poi
        print 'mv scan.pdf scan_%s.pdf'%poi
        print 'mv scan.png scan_%s.png'%poi
