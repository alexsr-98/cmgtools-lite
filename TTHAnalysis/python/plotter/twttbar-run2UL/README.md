# bb4l Run2

Last updated: **13 May 2025**

Index:
  * [Framework setup](#setup)
  * [Friend trees](#ftrees)
  * [bb4l differential cross section measurement (bkg. substraction)](#differential)
  * [bb4l differential cross section measurement (combine fit)](#differentialfit)


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
Friend trees are with the scripts in the folder `CMGTools/TTHAnalysis/macros`.

With prepareEventVariablesFriendTree.py we create the neccesary friend trees for the analysis. 

We have created a helper to automatise the friend tree creation process. You should always follow these four steps:

  * **Create** the friend trees chunks: `python produceFriendTrees_TopRun2UL.py -y 2018 -s 0 -d all -q batch -n 8`
  * **Check** if all chunks are created: `python produceFriendTrees_TopRun2UL.py -y 2018 -s 0 -d all -c`
  * **Merge** all chunks: `python produceFriendTrees_TopRun2UL.py -y 2018 -s 0 -d all -m`
  * **Check** if the merged chunks are ok: `python produceFriendTrees_TopRun2UL.py -y 2018 -s 0 -d all -m -c`


<a name="differential"></a>
## Plots:
 * Read the bash script, modify it to produce the desired plots and run it:
    ```bash
    source twttbar-run2UL/runANplots.sh
    ```


## bb4l differential cross section measurement (bkg. substraction)

 * Create the cards:
   ```bash
    python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o ./temp_Run2_cards_diff/2024-09-12_bb4l -v all -a
   ```
   Also possible to create them with the `runANplots.sh` script as for the plots.
 * Compute response matrices:
   ```bash
    python twttbar-run2UL/differential/getMatrices.py -i ./temp_Run2_cards_diff/2024-09-12_bb4l -y run2 -j 8
   ```
 * Extract signal:
   ```bash
    python twttbar-run2UL/differential/signalExtracter.py -i ./temp_Run2_cards_diff/2024-09-12_bb4l -y run2 -j 8
   ```
 * Unfold:
   ```bash
    python twttbar-run2UL/differential/unfoldHelper.py -i ./temp_Run2_cards_diff/2024-09-12_bb4l -y run2 -j 8
   ```
 * Normalise to the fiducial cross section and bin width:
   ```bash
    python twttbar-run2UL/differential/doFiducial.py -i ./temp_Run2_cards_diff/2024-09-12_bb4l -y run2 -j 8
   ```
 * Get the condition numbers for each matrix (for the Analysis Note):
   ```bash
   python twttbar-run2UL/differential/getLaTeXtable.py -i ./temp_Run2_cards_diff/2024-09-12_bb4l/  -t condnumtable
   ```
 * Get the GOF:
   ```bash
   python twttbar-run2UL/differential/goftests.py -i ./temp_Run2_cards_diff/2024-09-12_bb4l/
   ```


<a name="differentialfit"></a>
## bb4l differential cross section measurement (combine fit)
We need the cards created by the bkg. substraction scripts also. 

 * Create the cards:
   ```bash 
   python twttbar-run2UL/differential/cardsForDiffWithCombine.py -P 2024-09-03 -y run2 -q batch -j 100 -o ./temp_Run2_cards_diff/2025-03-04_combine/ -v Lep1Lep2_DPhi -a
   ```
   Also possible to create them with the `runANplots.sh` script as for the plots.
 * Extract signal and unfold. This will draw one uncertainty band extracted from combine:
   ```bash
   python twttbar-run2UL/differential/unfoldHelperWithCombine.py -i ./temp_Run2_cards_diff/2025-03-04_combine/ -y run2 -V -v Lep1Lep2_DPhi
   ```
 * To separate the uncertainties in individual sources (or simply in stat and sys) it is necessary to run individual scan freezing the systematics:
   ```bash
   python twttbar-run2UL/differential/toolkitForRelativeUncs_fitdiag.py -i ./temp_Run2_cards_diff/2025-03-04_combine/ -y run2 -j 1 -v Lep1Lep2_DPhi
   ```
 * Unfold again to draw the splitted uncertainties:
   ```bash
   python twttbar-run2UL/differential/unfoldHelperWithCombine.py -i ./temp_Run2_cards_diff/2025-03-04_combine/ -y run2 -V -v Lep1Lep2_DPhi
   ```
 * Normalise to the fiducial cross section and bin width:
   ```bash
   python twttbar-run2UL/differential/doFiducial_WithCombine.py -i ./temp_Run2_cards_diff/2025-03-04_combine/ -y run2 -v Lep1Lep2_DPhi
   ```