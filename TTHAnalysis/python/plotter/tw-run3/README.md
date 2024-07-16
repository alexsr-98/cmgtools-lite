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
```bash
source /cms/cmsset_default.sh
cmsenv
alias python=python3
```

Some commands in this readme need combine, you have to setup combine in its appropiate release.


<a name="ftrees"></a>
## Friend trees
With prepareEventVariablesFriendTree.py we create the neccesary friend trees for the analysis. 

We have created a helper to automatise the friend tree creation process. You should always follow these four steps:

  * **Create** the friend trees chunks: `python produceFriendTrees_TopRun3.py -y 2022PostEE -s 0 -d all -q batch -n 8`
  * **Check** if all chunks are created: `python produceFriendTrees_TopRun3.py -y 2022PostEE -s 0 -d all -c`
  * **Merge** all chunks: `python produceFriendTrees_TopRun3.py -y 2022PostEE -s 0 -d all -m`
  * **Check** if the merged chunks are ok: `python produceFriendTrees_TopRun3.py -y 2022PostEE -s 0 -d all -m -c`

<a name="inclusive"></a>
## tW Inclusive cross section measurement
 
 * To produce plots:
   ```bash
    python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o ./temp_Run3_plots/2023_06_08 -r nojets -u -q batch -j 16
   ```
 * To produce cards:
   ```bash
    python tw-run3/cardsHelper.py -P 2023-12-23 -y run3 -o ./temp_Run3_cards/2024-02-12_cardsIncl -r all -v all -a -q batch -j 16
   ```
    * Validate cards:
      ```bash
       ValidateDatacards.py temp_Run3_cards/2023-07-10_newSFs/2022PostEE/1j1t/cuts-tw-1j1t.txt
      ```
 * To make the fit:
   ```bash
    python tw-run3/fitsHelper.py -y 2022PostEE -i temp_Run3_cards/2023-06-02_Test/ -r 1j1t,2j1t,2j2t
   ```
 * To make impacts:
   ```bash
    python tw-run3/getInclusiveImpacts.py -y 2022PostEE -i temp_Run3_cards/2023-06-02_Test/ -r 1j1t,2j1t,2j2t -j 12
   ```
 * To make the uncertainty table:
   ```bash
   python ../../../../tw-run3/produceUncertaintyTable.py -i ../combcard_1j1t2j1t2j2t.root -s 0 -j 12
   ```
   ```bash
   python ../../../../tw-run3/produceUncertaintyTable.py -i ../combcard_1j1t2j1t2j2t.root -s 1 -j 12
   ```
 * To run GOF test:
   * Create the cards for the input variables and run the combine command: python tw-run3/createCardsInputVariables.py

 * To compute btag eff:
   ```bash
   python mcEfficiencies.py --tree NanoAOD  -P /lustrefs/hdd_pool_dir/nanoAODv12/tw-run3/productions/2023-12-23/2022PostEE/ --split-factor=-1 --year 2022PostEE --FMCs {P}/x_btageff  tw-run3/mca-tw-includes/mca-2022PostEE-tw-btageff.txt tw-run3/cuts-tw-btageff.txt  tw-run3/plots-tw/plots-tw-nojets_btageffsels.txt tw-run3/plots-tw/plots-tw-nojets_btageffvars.txt  -o temp_Run3_plots/2024_02_28_btagEff_tW_comb/2022PostEE/eff/output.root -j 12
   ```
   ```bash
   python mcEfficiencies.py --tree NanoAOD  -P /lustrefs/hdd_pool_dir/nanoAODv12/tw-run3/productions/2023-12-23/2022/ --split-factor=-1 --year 2022 --FMCs {P}/x_btageff  tw-run3/mca-tw-includes/mca-2022-tw-btageff.txt tw-run3/cuts-tw-btageff.txt  tw-run3/plots-tw/plots-tw-nojets_btageffsels.txt tw-run3/plots-tw/plots-tw-nojets_btageffvars.txt  -o temp_Run3_plots/2024_02_28_btagEff_tW_comb/2022/eff/output.root -j 12
   ```
 * To produce 2D SFs plots:
   ```bash
    python tw-run3/plotSFhistograms.py ./temp_Run3_plots/2023_01_21_Plots_de_SFs_ttbarRun3
   ```
 * To produce uncertainty variations plots:
   ```bash
    python tw-run3/plotUncsVariations.py temp_Run3_cards/2024-02-12_cardsIncl/run3/1j1t/ -j 8
   ```
 * To produce the electron and jet veto maps:
   ```bash
    python tw-run3/utils/checkVetoMaps.py
   ```
 
 * To produce trigger SFs:
   ```bash
    python tw-run3/plotterHelper.py -P 2023-10-25 -y 2022PostEE -o ./temp_Run3_plots/2023_06_05_trigger -r trigger_METandLeptonTrig -q batch -j 16
   ```
   ```bash
   python tw-run3/plotterHelper.py -P 2023-10-25 -y 2022PostEE -o ./temp_Run3_plots/2023_06_05_trigger -r trigger_METTrig -q batch -j 16
   ```
   ```bash
   python tw-run3/plotterHelper.py -P 2023-10-25 -y 2022PostEE -o ./temp_Run3_plots/2023_12_10_trigger -r trigger_LeptonTrig -q batch -j 16
   ```
   ```bash
   python tw-run3/plotterHelper.py -P 2023-10-25 -y 2022PostEE -o ./temp_Run3_plots/2023_12_10_trigger -r trigger_None -q batch -j 16
   ```
   ```bash
   python tw-run3/triggerSFs/makeEffandSF.py -i temp_Run3_plots/2023_12_10_trigger/2022PostEE/ -o tw-run3/triggerSFs/SFs-plots_2023_12_10/ -c trigger -y 2022PostEE
   ```

 ## Train MVA

The script to train is under `MVA-training` folder. To use your model in the analysis, you need to do the following:
 * Train and save the model with the training script. The model should be saved in joblib format.
 * Convert the joblib model to ONNX (there is a dedicated script to do this). It is located in `MVA-training/onnxConverter/`. Follow instructions the readme instructions in that folder.
 * Use the ONNX model to create the friend trees.
 * To train the models: `python tw-run3/MVA-Training/tW_MVA_trainer.py -M MultiRF -r 1j1b -o MultiRF_1j1b_2023_06_04 -j 12`

<a name="differential"></a>
## tW differential cross section measurement (bkg. substraction)

 * Create the cards:
   ```bash
    python tw-run3/differential/cardsForDifferentialStudies.py -P 2023-12-23 -y run3 -q batch -j 16 -o ./temp_cards_diff/2023-11-22 -v all -a
   ```
 * Compute response matrices:
   ```bash
    python tw-run3/differential/getMatrices.py -i ./temp_cards_diff/2023-06-09 -y run3 -j 8
   ```
 * Extract signal:
   ```bash
    python tw-run3/differential/signalExtracter.py -i ./temp_cards_diff/2023-06-09 -y run3 -j 8
   ```
 * Unfold:
   ```bash
    python tw-run3/differential/unfoldHelper.py -i ./temp_cards_diff/2023-06-09 -y run3 -j 8
   ```
 * Normalise:
   ```bash
    python tw-run3/differential/doFiducial.py -i ./temp_cards_diff/2023-06-09 -y run3 -j 8
   ```
 * Get latexTableCondNum:
   ```bash
   python tw-run3/differential/getLaTeXtable.py -i ./temp_cards_diff/2023-11-28_differential/  -t condnumtable
   ```
 * Get the GOF:
   ```bash
   python tw-run3/differential/goftests.py -i ./temp_cards_diff/2024-03-09_cardsDiff_ARCv1_test/
   ```

There is a subfolder under tw-run3 called utils. It contains several scripts useful for different tasks, they are:
 * copyTo_www.py: copy a folder to your website and adds .php files for visualisation in the browser.
 * replacepdftext.sh: it takes a plot in pdf style and substitutes one string by another. Useful for removing preliminary from plots.

<a name="differentialfit"></a>
## tW differential cross section measurement (combine fit)

 * Create the cards. There are two steps, first you have to create the cards using the same script as above but now as the signal extraction will be done with combine, we don't need the `forExtr` region. Second, you have to run the script `cardsForDifferentialStudies_WithCombine.py` to create the card for the signal extraction. Example:
   ```bash
    python tw-run3/differential/cardsForDifferentialStudies.py -P 2023-12-23 -y run3 -q batch -j 32 -o ./temp_cards_diff/2024-02-24_WithCombine -v all -a
   python tw-run3/differential/cardsForDifferentialStudies_WithCombine.py -P 2023-12-23 -y run3 -q batch -j 32 -o ./temp_cards_diff/2024-02-24_WithCombine -v all -a  
   ```
   Now, we also don't need the `Fiducial` variable.
  
 * Extract signal and unfold (done at the same time with the fit). Remember to setup combine (in its appropiate release).
   ```bash
     python tw-run3/differential/unfoldHelper_WithCombine.py -i ./temp_cards_diff/2023-10-03_WithCombine/ -y run3 -V
   ```
 
 * Estimate the impact of each source separately, and the run again the unfolding. 
   ```bash
     python tw-run3/differential/toolkitForRelativeUncs_WithCombine.py -i ./temp_cards_diff/2023-10-03_WithCombine/ -y run3 -j 1
   ```
 
 * Normalise to the fiducial cross section.
    ```bash
      python tw-run3/differential/doFiducial_WithCombine.py -i ./temp_cards_diff/2023-10-03_WithCombine/ -y run3
    ```

 <a name="samples"></a>
## Samples
 [GoogleSpreadSheet](https://docs.google.com/spreadsheets/d/1ajj4HEToMIxzCNrgyGBdPO0e-UniPjZSXWibZdBAO00/edit#gid=0)

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
 * Muon SFs: [twiki](https://twiki.cern.ch/twiki/bin/view/CMS/MuonRun32022).
 * BTV POG: [web](https://btv-wiki.docs.cern.ch/ScaleFactors/).
 * JSON LUMI (PU): [twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData).
 * Stat. Comm. recommendations on Unfolding: [twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/ScrecUnfolding)(https://indico.cern.ch/event/1311191/timetable/#20230927). 
 * bb4l gridpack instructions: [codiMD](https://codimd.web.cern.ch/58DeTshTTheTlRcNGk81bQ?both).
 * Compute the diff between two versions of the paper: [tdrDiff](https://cms-tdr-diff.web.cern.ch/).
 * Checklist and links for PAS preparation: [twiki-PASchecklist](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PCPasPaperChecklists), [mathJax](http://genkuroki.web.fc2.com/MathJax/LivePreviewMathJax-jquery.html) and [checkLatexBuild](https://icms.cern.ch/tools-api/cadi/checkLatexBuild).
 * ARC review steps: [twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/ArcCollaborationWideReview#Overview_of_the_publication_revi).
