# First step is to run the plotterHelper to get the plots:
OUTDIR="./temp_Run3_plots/2024_01_06_trigger_allErasFixed/"
# Create a list to loop over with the different regions
declare -a regions=("trigger_METandLeptonTrig" "trigger_METTrig" "trigger_METTrigMore3j" "trigger_METTrigLess3j" "trigger_METandLeptonTrigMore3j" "trigger_METandLeptonTrigLess3j" "trigger_METandLeptonTrigLess35vtx" "trigger_METandLeptonTrigMore35vtx" "trigger_METTrigLess35vtx" "trigger_METTrigMore35vtx" "trigger_LeptonTrig" "trigger_None" "trigger_DoubleLeptonTrigNoMETcut" "trigger_SingleLeptonTrigNoMETcut" "trigger_LeptonTrigNoMETcut" "trigger_NoneNoMETcut" "trigger-ee_METandLeptonTrig" "trigger-ee_METTrig" "trigger-ee_METandLeptonTrigMore3j" "trigger-ee_METandLeptonTrigLess3j" "trigger-ee_METTrigMore3j" "trigger-ee_METTrigLess3j" "trigger-ee_METandLeptonTrigLess35vtx" "trigger-ee_METandLeptonTrigMore35vtx" "trigger-ee_METTrigLess35vtx" "trigger-ee_METTrigMore35vtx" ""trigger-ee_LeptonTrig"" "trigger-ee_None" "trigger-mm_METandLeptonTrig" "trigger-mm_METTrig" "trigger-mm_METandLeptonTrigMore3j" "trigger-mm_METandLeptonTrigLess3j" "trigger-mm_METTrigMore3j" "trigger-mm_METTrigLess3j" "trigger-mm_METandLeptonTrigLess35vtx" "trigger-mm_METandLeptonTrigMore35vtx" "trigger-mm_METTrigLess35vtx" "trigger-mm_METTrigMore35vtx" "trigger-mm_LeptonTrig" "trigger-mm_None")

# Loop over the regions
for region in "${regions[@]}"
do
    python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022PostEE -o $OUTDIR -r $region -q batch -j 16
    python tw-run3/plotterHelper.py -P 2023-12-23 -y 2022 -o $OUTDIR -r $region -q batch -j 16
done