<a name="differential"></a>
## tW differential cross section measurement (bkg. substraction)

 * Create the cards:
   ```bash
    python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o ./temp_Run2_cards_diff/2024-09-12_bb4l -v all -a
   ```
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
 * Normalise:
   ```bash
    python twttbar-run2UL/differential/doFiducial.py -i ./temp_Run2_cards_diff/2024-09-12_bb4l -y run2 -j 8
   ```
 * Get latexTableCondNum:
   ```bash
   python twttbar-run2UL/differential/getLaTeXtable.py -i ./temp_Run2_cards_diff/2024-09-12_bb4l/  -t condnumtable
   ```
 * Get the GOF:
   ```bash
   python twttbar-run2UL/differential/goftests.py -i ./temp_Run2_cards_diff/2024-09-12_bb4l/
   ```