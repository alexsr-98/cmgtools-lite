OUTDIR="./temp_Run2_plots/2024_06_24_forAN/"
OUTDIR2="./temp_Run2_plots/2024_06_24_noUnc/"
OUTDIR3="./temp_Run2_plots/2024_07_11_noUncbb4l/"
OUTDIRDIFF="./temp_Run2_cards_diff/2024-07-15_newbb4l_nojets/"


### Cards for the differential measurement
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1Lep2_DPhi -f -a
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v minimax_ATLAS -f -a
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Fiducial -f -a
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Jet1_Eta -f -a
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Jet1_Pt -f -a
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1_Pt -f -a
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1_Eta -f -a
#
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1Lep2_DPhi -f -a -njnt 1j1t
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Fiducial -f -a -njnt 1j1t
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Jet1_Eta -f -a -njnt 1j1t
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Jet1_Pt -f -a -njnt 1j1t
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1_Pt -f -a -njnt 1j1t
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1_Eta -f -a -njnt 1j1t

#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1Lep2_DPhi -f -a -njnt nojets
python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Fiducial -f -a -njnt nojets
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Jet1_Pt -f -a -njnt nojets
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1_Pt -f -a -njnt nojets
python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v NJets -f -a -njnt nojets
python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2022-05-06 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v NBJets -f -a -njnt nojets

# With uncertainties
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y run2 -o $OUTDIR -r nojets -q batch -j 16 -u -f --sP nJetnBJet
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y run2 -o $OUTDIR -r 2j2t -q batch -j 16 -u -f --sP minimax --sP jet1_pt --sP jet2_pt --sP lep1_pt --sP lep2_pt

# Without uncertainties
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y run2 -o $OUTDIR2 -r nojets -q batch -j 16 -f
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y run2 -o $OUTDIR2 -r nojets-ee -q batch -j 16 -f
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y run2 -o $OUTDIR2 -r nojets-mm -q batch -j 16 -f
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y run2 -o $OUTDIR2 -r 1j1t -q batch -j 16 -f
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y run2 -o $OUTDIR2 -r 2j1t -q batch -j 16 -f
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y run2 -o $OUTDIR2 -r 2j2t -q batch -j 16 -f --sP minimax

# Comparison with bb4l
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y run2 -o $OUTDIR3 -r nojets -q batch -j 32 -f --comparison
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y run2 -o $OUTDIR3 -r 1j1t -q batch -j 32 -f --comparison
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y run2 -o $OUTDIR3 -r 2j1t -q batch -j 32 -f --comparison
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y run2 -o $OUTDIR3 -r 2j2t -q batch -j 32 -f --comparison

#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y 2016 -o $OUTDIR3 -r 2j2t -q batch -j 16 -f --comparison --sP minimax
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y 2016apv -o $OUTDIR3 -r 2j2t -q batch -j 16 -f --comparison --sP minimax
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y 2017 -o $OUTDIR3 -r 2j2t -q batch -j 16 -f --comparison --sP minimax
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y 2018 -o $OUTDIR3 -r 2j2t -q batch -j 16 -f --comparison --sP minimax
#python twttbar-run2UL/plotterHelper.py -P 2022-05-06 -y run2 -o $OUTDIR3 -r 2j2t -q batch -j 16 -f --comparison --sP minimax

