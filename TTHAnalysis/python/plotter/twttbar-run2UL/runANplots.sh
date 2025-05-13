OUTDIR="./temp_Run2_plots/2024_10_04_forAN/"
OUTDIR2="./temp_Run2_plots/2024_09_10_noUnc/"
OUTDIR3="./temp_Run2_plots/2025_02_14/"
OUTDIRDIFF="./temp_Run2_cards_diff/2025-02-12_cardsDiff/"

#### Cards for the differential with combine
python twttbar-run2UL/differential/cardsForDiffWithCombine.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1_Pt -a
python twttbar-run2UL/differential/cardsForDiffWithCombine.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1Lep2_DPhi -a
python twttbar-run2UL/differential/cardsForDiffWithCombine.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Jet1_Pt -a
python twttbar-run2UL/differential/cardsForDiffWithCombine.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v minimax_ATLAS -a

### Cards for the differential measurement
python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1Lep2_DPhi -a 
python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v minimax_ATLAS -a 
python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Fiducial -a 
python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Jet1_Pt -a 
python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1_Pt -a 
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1_Eta -a
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Jet1_Eta -a





#
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1Lep2_DPhi -a -njnt 1j1t
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Fiducial -a -njnt 1j1t
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Jet1_Eta -a -njnt 1j1t
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Jet1_Pt -a -njnt 1j1t
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1_Pt -a -njnt 1j1t
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1_Eta -a -njnt 1j1t

#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1Lep2_DPhi -a -njnt nojets
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Fiducial -a -njnt nojets
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Jet1_Pt -a -njnt nojets
#python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v Lep1_Pt -a -njnt nojets
python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v NJets -a -njnt nojets
python twttbar-run2UL/differential/cardsForDifferentialStudies.py -P 2024-09-03 -y run2 -q batch -j 16 -o $OUTDIRDIFF -v NBJets  -a -njnt nojets

# With uncertainties
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR -r nojets -q batch -j 64 -u --sP nJetnBJet --sP nJets --sP nBJets --comparison
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR -r nojets -q batch -j 64 -u --sP nJets --sP nBJets --comparison
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR -r 1j1t -q batch -j 64 -u --sP jet1_pt --sP lep1_pt --sP lep2_pt --sP lep1lep2_dphi --sP lep1lep2_m --sP tot_weight --comparison
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR -r 2j1t -q batch -j 64 -u --sP jet1_pt --sP jet2_pt --sP lep1_pt --sP lep2_pt --sP lep1lep2_dphi --sP lep1lep2_m --sP tot_weight --comparison
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR -r 2j2t -q batch -j 64 -u --sP jet1_pt --sP jet2_pt --sP lep1_pt --sP lep2_pt --sP lep1lep2_dphi --sP lep1lep2_m --sP tot_weight --comparison

# Comparison with bb4l
python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR3 -r nojets -q batch -j 64  --comparison
python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR2 -r nojets-ee -q batch -j 16 
python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR2 -r nojets-mm -q batch -j 16
python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR3 -r 1j1t -q batch -j 64  --comparison
python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR3 -r 2j1t -q batch -j 16  --comparison
python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR3 -r 2j2t -q batch -j 64  --comparison
python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR3 -r 2j2t-ee -q batch -j 64  --comparison
python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR3 -r 2j2t-mm -q batch -j 64  --comparison
#
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2016apv -o $OUTDIR3 -r nojets -q batch -j 16  --comparison --sP nJetnBJet
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2016apv -o $OUTDIR3 -r 2j2t -q batch -j 16  --comparison
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2016 -o $OUTDIR3 -r nojets -q batch -j 16  --comparison --sP nJetnBJet
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2016 -o $OUTDIR3 -r 2j2t -q batch -j 16  --comparison
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2017 -o $OUTDIR3 -r nojets -q batch -j 16  --comparison --sP nJetnBJet
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2017 -o $OUTDIR3 -r 2j2t -q batch -j 16  --comparison
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2018 -o $OUTDIR3 -r nojets -q batch -j 16  --comparison --sP nJetnBJet
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2018 -o $OUTDIR3 -r 2j2t -q batch -j 16  --comparison

#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2016apv -o $OUTDIR3 -r 2j2t -q batch -j 16 --comparison --sP minimax 
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2016 -o $OUTDIR3 -r 2j2t -q batch -j 16 --comparison  --sP minimax
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2017 -o $OUTDIR3 -r 2j2t -q batch -j 16 --comparison  --sP minimax 
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2018 -o $OUTDIR3 -r 2j2t -q batch -j 16  --comparison --sP minimax 
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR3 -r 2j2t -q batch -j 16 --comparison --sP minimax

#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2016apv -o $OUTDIR3 -r 2j2t-mm -q batch -j 16 --comparison --sP minimax 
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2016 -o $OUTDIR3 -r 2j2t-mm -q batch -j 16 --comparison  --sP minimax
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2017 -o $OUTDIR3 -r 2j2t-mm -q batch -j 16 --comparison  --sP minimax 
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2018 -o $OUTDIR3 -r 2j2t-mm -q batch -j 16  --comparison --sP minimax 
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR3 -r 2j2t-mm -q batch -j 16 --comparison --sP minimax

#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2016apv -o $OUTDIR3 -r 2j2t-ee -q batch -j 16 --comparison --sP minimax 
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2016 -o $OUTDIR3 -r 2j2t-ee -q batch -j 16 --comparison  --sP minimax
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2017 -o $OUTDIR3 -r 2j2t-ee -q batch -j 16 --comparison  --sP minimax 
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2018 -o $OUTDIR3 -r 2j2t-ee -q batch -j 16  --comparison --sP minimax 
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR3 -r 2j2t-ee -q batch -j 16 --comparison --sP minimax



#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y 2016 -o $OUTDIR3 -r 3j3t -q batch -j 16 -f --sP minimax



##### MET Study
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR3 -r 2j2t -q batch -j 16  --comparison --sP met --sP met_phi --sP puppimet --sP puppimet_phi
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR3 -r 2j2t-ee -q batch -j 16  --comparison --sP met --sP met_phi --sP puppimet --sP puppimet_phi --sP lep1lep2_m
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR3 -r 2j2t-mm -q batch -j 16  --comparison --sP met --sP met_phi --sP puppimet --sP puppimet_phi --sP lep1lep2_m


#### Study of Nvertex 
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR3 -r nojets -q batch -j 64  --normaliseToData
#python twttbar-run2UL/plotterHelper.py -P 2024-09-03 -y run2 -o $OUTDIR3 -r 2j2t -q batch -j 64  --normaliseToData