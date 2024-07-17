OUTDIR="./temp_Run3_plots/2024_06_13_forAN_CWR/"
OUTDIR2="./temp_Run3_plots/2024_04_09_testCarlosWithoutSFs/"
OUTDIR3="./temp_Run3_plots/2024_03_18_PUvars_FullWeightFixTestMoreVars/"

#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR -r nojets -q batch -j 16 -u --sP nJetnBJet
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR -r 1j1t -q batch -j 16 -u --sP jet1_pt --sP lep1_pt --sP lep1lep2_pt --sP lep1lep2_ptsum --sP lep1lep2_dphi --sP lep1lep2_m --sP tot_weight --sP nloosejets
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR -r 2j1t -q batch -j 16 -u --sP jet1_pt --sP lep1_pt --sP lep1lep2_pt --sP lep1lep2_ptsum --sP lep1lep2_dphi --sP lep1lep2_m --sP tot_weight
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR -r 2j2t -q batch -j 16 -u --sP jet1_pt --sP lep1_pt --sP lep1lep2_pt --sP lep1lep2_ptsum --sP lep1lep2_dphi --sP lep1lep2_m --sP tot_weight
##python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR -r 1j1t_MVAtrain -q batch -j 16 -u --sP loosejet1_pt --sP lep1_pt --sP lep1lep2jet1_pt --sP lep1lep2_m --sP lep1lep2_dphi --sP lep1lep2jet1_m --sP lep1jet1_pt --sP jet1_pt
##python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR -r 2j1t_MVAtrain -q batch -j 16 -u --sP lep1lep2_m --sP lep2_pt --sP lep12jet12_dr --sP lep1lep2jet1_pt --sP lep1lep2_dr --sP lep1jet1_dr --sP lep1lep2jet1_c --sP jet2_pt
##python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR -r 1j1t_differential -q batch -j 16 -u

#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR -r nojets -q batch -j 16 -u --sP nJetnBJet
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR -r 1j1t -q batch -j 16 -u --sP jet1_pt --sP lep1_pt --sP lep1lep2_pt --sP lep1lep2_ptsum --sP lep1lep2_dphi --sP lep1lep2_m --sP tot_weight --sP nloosejets
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR -r 2j1t -q batch -j 16 -u --sP jet1_pt --sP lep1_pt --sP lep1lep2_pt --sP lep1lep2_ptsum --sP lep1lep2_dphi --sP lep1lep2_m --sP tot_weight
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR -r 2j2t -q batch -j 16 -u --sP jet1_pt --sP lep1_pt --sP lep1lep2_pt --sP lep1lep2_ptsum --sP lep1lep2_dphi --sP lep1lep2_m --sP tot_weight
##python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR -r 1j1t_MVAtrain -q batch -j 16 -u --sP loosejet1_pt --sP lep1_pt --sP lep1lep2jet1_pt --sP lep1lep2_m --sP lep1lep2_dphi --sP lep1lep2jet1_m --sP lep1jet1_pt --sP jet1_pt
##python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR -r 2j1t_MVAtrain -q batch -j 16 -u --sP lep1lep2_m --sP lep2_pt --sP lep12jet12_dr --sP lep1lep2jet1_pt --sP lep1lep2_dr --sP lep1jet1_dr --sP lep1lep2jet1_c --sP jet2_pt
##python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR -r 1j1t_differential -q batch -j 16 -u

#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR -r nojets -q batch -j 16 -u --sP nJetnBJet
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR -r 1j1t -q batch -j 16 -u --sP jet1_pt --sP lep1_pt --sP lep1lep2_pt --sP lep1lep2_ptsum --sP lep1lep2_dphi --sP lep1lep2_m --sP tot_weight --sP nloosejets
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR -r 2j1t -q batch -j 16 -u --sP jet1_pt --sP lep1_pt --sP lep1lep2_pt --sP lep1lep2_ptsum --sP lep1lep2_dphi --sP lep1lep2_m --sP tot_weight
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR -r 2j2t -q batch -j 16 -u --sP jet1_pt --sP lep1_pt --sP lep1lep2_pt --sP lep1lep2_ptsum --sP lep1lep2_dphi --sP lep1lep2_m --sP tot_weight
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR -r 1j1t_MVAtrain -q batch -j 32 -u --sP loosejet1_pt --sP lep1_pt --sP lep1lep2jet1_pt --sP lep1lep2_m --sP lep1lep2_dphi --sP lep1lep2jet1_m --sP lep1jet1_pt --sP jet1_pt --sP mvaRF_1j1b
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR -r 2j1t_MVAtrain -q batch -j 32 -u --sP lep1lep2_m --sP lep2_pt --sP lep12jet12_dr --sP lep1lep2jet1_pt --sP lep1lep2_dr --sP lep1jet1_dr --sP lep1lep2jet1_c --sP jet2_pt --sP mvaRF_2j1b
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR -r 1j1t_differential -q batch -j 16 -u

# Without uncertainties
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR2 -r nojets -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR2 -r 1j1t -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR2 -r 2j1t -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR2 -r 2j2t -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR2 -r 1j1t-mm -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR2 -r 1j1t-ee -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR2 -r nojets-ee -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR2 -r nojets-mm -q batch -j 16
#
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR2 -r nojets -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR2 -r 1j1t -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR2 -r 2j1t -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR2 -r 2j2t -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR2 -r 1j1t-mm -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR2 -r 1j1t-ee -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR2 -r nojets-ee -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR2 -r nojets-mm -q batch -j 16
#
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR2 -r nojets -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR2 -r 1j1t -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR2 -r 2j1t -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR2 -r 2j2t -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR2 -r 1j1t-mm -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR2 -r 1j1t-ee -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR2 -r nojets-ee -q batch -j 16
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR2 -r nojets-mm -q batch -j 16

# Test of the PU reweighting
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR3 -r nojets -q batch -j 32 --sP nvertex --sP nvertexnotgood --sP rho_pu --sP rho_calo --sP rho_neutral --sP nJetnBJet
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR3 -r nojets -q batch -j 32 --sP nvertex --sP nvertexnotgood --sP rho_pu --sP rho_calo --sP rho_neutral --sP nJetnBJet
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR3 -r nojets -q batch -j 16 -u --sP nvertex --sP nvertexnotgood --sP rho_pu --sP rho_calo --sP rho_neutral

#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR3 -r 1j1t_MVAtrain -q batch -j 64 -u --sP loosejet1_pt --sP lep1_pt --sP lep1lep2jet1_pt --sP lep1lep2_m --sP lep1lep2_dphi --sP lep1lep2jet1_m --sP lep1jet1_pt --sP jet1_pt
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR3 -r nojets -q batch -j 64 -u --sP jet1_pt --sP lep1_pt --sP lep1lep2_pt --sP lep1lep2_ptsum --sP lep1lep2_dphi --sP lep1lep2_m


# Plots for paper only
python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR -r nojets -q batch -j 16 -u --sP nJetnBJet -f --cmsText ''
python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR -r 1j1t -q batch -j 16 -u --sP nloosejets -f --cmsText ''
python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR -r 1j1t_MVAtrain -q batch -j 32 -u --sP loosejet1_pt --sP lep1_pt --sP lep1lep2jet1_pt --sP lep1lep2_m -f --cmsText ''
python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR -r 1j1t_differential -q batch -j 32 -u -f --cmsText ''
python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR -r 2j1t_MVAtrain -q batch -j 32 -u --sP lep1lep2_m --sP lep2_pt --sP lep12jet12_dr --sP lep1lep2jet1_pt -f --cmsText 'Supplementary'




#### Test for carlos
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR2 -r nojets-ee -q batch -j 16 --sP lep1lep2_m_picu --sP lep1_pt --sP lep1_eta --sP lep2_pt --sP lep2_eta
#python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR2 -r nojets-ee -q batch -j 16 --sP  lep1lep2_m_picu --sP lep1_pt --sP lep1_eta --sP lep2_pt --sP lep2_eta
#python tw-run3/plotterHelper.py -P 2023-12-23 -y run3 -o $OUTDIR2 -r nojets-ee -q batch -j 16 --sP lep1lep2_m_picu --sP lep1_pt --sP lep1_eta --sP lep2_pt --sP lep2_eta