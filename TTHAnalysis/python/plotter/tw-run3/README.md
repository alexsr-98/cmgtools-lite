# tW Run3

Last updated: **18 Sept 2023**

(Documentation work in progress)
  * [Framework setup](#setup)
  * [Friend trees](#ftrees)
  * [tW Inclusive cross section measurement](#inclusive)
  * [tW differential cross section measurement (bkg. substraction)](#differential)
  * [tW differential cross section measurement (combine fit)](#differentialfit)
  * [Samples](#samples)
  * [Important twikis](#twikis)

<a name="setup"></a>
## Framework setup
Three commands:
source /cms/cmsset_default.sh
cmsenv
alias python=python3

To run the commands in this readme that need combine, you need to setup combine in its respective release.


<a name="ftrees"></a>
## Friend trees
With prepareEventVariablesFriendTree.py we create the neccesary friend trees for the analysis. 

We have created a helper to automatise the friend tree creation process: . You should always follow these four steps:

  * Create the friend trees chunks: `python produceFriendTrees_TopRun3.py -y 2022 -s 0 -d all -q batch -n 8`
  * Check if all chunks are created: `python produceFriendTrees_TopRun3.py -y 2022 -s 0 -d all -c`
  * Merge all chunks: `python produceFriendTrees_TopRun3.py -y 2022 -s 0 -d all -m`
  * Check if the merged chunks are ok: `python produceFriendTrees_TopRun3.py -y 2022 -s 0 -d all -m -c`

<a name="inclusive"></a>
## tW Inclusive cross section measurement
 
 * To produce plots:
    python tw-run3/plotterHelper.py -P 2023-06-02 -y 2022PostEE -o ./temp_Run3_plots/2023_06_08 -r nojets -u -q batch -j 32
 * To produce cards:
    python tw-run3/cardsHelper.py -P 2023-06-02 -y 2022PostEE -o ./temp_Run3_cards/2023-06-02_Test -r all -v all -a -q batch -j 16
    * Validate cards:
       ValidateDatacards.py temp_Run3_cards/2023-07-10_newSFs/2022PostEE/1j1t/cuts-tw-1j1t.txt
 * To make the fit:
    python tw-run3/fitsHelper.py -y 2022PostEE -i temp_Run3_cards/2023-06-02_Test/ -r 1j1t,2j1t,2j2t
 * To make impacts:
    python tw-run3/getInclusiveImpacts.py -y 2022PostEE -i temp_Run3_cards/2023-06-02_Test/ -r 1j1t,2j1t,2j2t -j 12

 * To compute btag eff:
    python mcEfficiencies.py --tree NanoAOD  -P /beegfs/data/nanoAODv11/tw-run3/productions/2023-06-02/2022PostEE/ --split-factor=-1 --year 2022 --FMCs {P}/x_btageff  tw-run3/mca-tw-includes/mca-2022PostEE-tw-btageff.txt tw-run3/cuts-tw-btageff.txt  tw-run3/plots-tw/plots-tw-nojets_btageffsels.txt tw-run3/plots-tw/plots-tw-nojets_btageffvars.txt  -o temp_Run3_plots/2022_11_22_addedLeptonSFs_correctedJson_addedUnc/2022/eff/output.root -j 12

 * To produce 2D SFs plots:
    python tw-run3/plotSFhistograms.py ./temp_Run3_plots/2023_01_21_Plots_de_SFs_ttbarRun3
 * To produce uncertainty variations plots:
    python tw-run3/plotUncsVariations.py temp_Run3_cards/2023-06-28_newSFsBarbara/2022PostEE/2j2t/ -j 8
 * To produce the electron and jet veto maps:
    python tw-run3/utils/checkVetoMaps.py
 
 * To produce trigger SFs:
    python tw-run3/plotterHelper.py -P 2023-06-02 -y 2022PostEE -o ./temp_Run3_plots/2023_06_05_trigger -r trigger_METandLeptonTrig -q batch -j 64
    python tw-run3/plotterHelper.py -P 2023-06-02 -y 2022PostEE -o ./temp_Run3_plots/2023_06_05_trigger -r trigger_METTrig -q batch -j 64

    python tw-run3/triggerSFs/makeEffandSF.py -i temp_Run3_plots/2023_06_07_Trigger/2022PostEE/ -o tw-run3/triggerSFs/SFs-plots/ -c trigger

 ## Train MVA

The script to train is under `MVA-training` folder. To use your model in the analysis, you need to do the following:
 * Train and save the model with the training script. The model should be saved in joblib format.
 * Convert the joblib model to ONNX (there is a dedicated script to do this). It is located in `MVA-training/onnxConverter/`. Follow instructions the readme instructions in that folder.
 * Use the ONNX model to create the friend trees.
 * To train the models: `python tw-run3/MVA-Training/tW_MVA_trainer.py -M MultiRF -r 1j1b -o MultiRF_1j1b_2023_06_04 -j 12`

<a name="differential"></a>
## tW differential cross section measurement (bkg. substraction)

 * Create the cards:
    python tw-run3/differential/cardsForDifferentialStudies.py -P 2023-06-02 -y 2022PostEE -q batch -j 32 -o ./temp_cards_diff/2023-06-09 -v 'Lep1_Pt' -a
 * Compute response matrices:
    python tw-run3/differential/getMatrices.py -i ./temp_cards_diff/2023-06-09 -y 2022PostEE -j 8
 * Extract signal:
    python tw-run3/differential/signalExtracter.py -i ./temp_cards_diff/2023-06-09 -y 2022PostEE -j 8
 * Unfold:
    python tw-run3/differential/unfoldHelper.py -i ./temp_cards_diff/2023-06-09 -y 2022PostEE -j 8
 * Normalise:
    python tw-run3/differential/doFiducial.py -i ./temp_cards_diff/2023-06-09 -y 2022PostEE -j 8

There is a subfolder under tw-run3 called utils. It contains several scripts useful for different tasks, they are:
 * copyTo_www.py: copy a folder to your website and adds .php files for visualisation in the browser.
 * replacepdftext.sh: it takes a plot in pdf style and substitutes one string by another. Useful for removing preliminary from plots.

<a name="differentialfit"></a>
## tW differential cross section measurement (combine fit)
WIP

 <a name="samples"></a>
## Samples
 [link](https://docs.google.com/spreadsheets/d/1ajj4HEToMIxzCNrgyGBdPO0e-UniPjZSXWibZdBAO00/edit#gid=0)

 <a name="twikis"></a>
## Important twikis
Here we will list the important twikis for the analysis:
 * Analysis [twiki](https://twiki.cern.ch/twiki/bin/view/CMS/TWRun3).
 * Top systematics [twiki](https://twiki.cern.ch/twiki/bin/view/CMS/TopSystematics).
 * Cross sections [twiki](https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopNNLORef).
 * Top HEPData [twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/TOPHepData).
 * PdMV Run3 [twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis).
 * JetID Run3 [twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13p6TeV).
 * NanoAODv11 [twiki](https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv11/2022postEE/doc_WZ_TuneCP5_13p6TeV_pythia8_Run3Summer22EENanoAODv11-126X_mcRun3_2022_realistic_postEE_v1-v1.html).
 * Json integration POGs: [gitlab](https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG).
 * Lumi Recommendation: [twiki](https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun3).
 * EGamma SFs: [twiki](https://twiki.cern.ch/twiki/bin/view/CMS/EgammSFandSSRun3).
 * BTV POG: [web](https://btv-wiki.docs.cern.ch/ScaleFactors/).
 * JSON LUMI (PU): [twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData).
 * Stat. Comm. recommendations on Unfolding: [twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/ScrecUnfolding).