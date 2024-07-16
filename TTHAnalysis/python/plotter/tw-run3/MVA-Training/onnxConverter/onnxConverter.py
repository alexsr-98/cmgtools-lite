from skl2onnx import convert_sklearn
from skl2onnx.common.data_types import FloatTensorType
import joblib

#path = '/nfs/fanae/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/MVA-Training/models/1j1b_RF/RandomForest/'
#
#model=joblib.load(path + 'RandomForest.joblib')
#vars = ["train_nloosejets", 
#        "train_jet1_pt",
#        "train_loosejet1_pt",
#        "train_lep1lep2jet1met_m",
#        "train_lep1lep2jet1_c",
#        "train_lep1lep2jet1_pt"]
#vars = ["train_jet2_pt",
#        "train_lep1jet1_dr",
#        "train_lep12jet12_dr"]

### New 1j1b vars (tests)
#path = '/nfs/fanae/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/MVA-Training/models/1j1b_RF_testDifferentVars/RandomForest/'
#vars = ["train_lep1jet1_pt", "train_lep1lep2jet1met_m", "train_lep1lep2jet1_pt", "train_lep1lep2jet1met_ptOVERhttot", "train_htlepOVERhttot", "train_lep1jet1_dr", "train_lep1lep2jet1met_mt"]

#### Multi-class (CURRENT)
path = '/mnt_pool/c3_users/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/MVA-Training/models/MultiRF_1j1b_2024_01_07/1j1b/'
vars = ["train_loosejet1_pt", "train_lep1lep2jet1_pt", "train_lep1lep2_m", "train_lep1lep2_dphi", "train_lep1jet1_pt", "train_lep1_pt", "train_jet1_pt", "train_lep1lep2jet1_m"]

path = '/mnt_pool/c3_users/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/MVA-Training/models/MultiRF_2j1b_2024_01_07/2j1b/'
vars = ["train_lep1lep2_m", "train_lep12jet12_dr", "train_lep1lep2jet1_pt", "train_lep1jet1_dr", "train_lep1lep2_dr", "train_lep2_pt", "train_jet2_pt", "train_lep1lep2jet1_c"]

# Multi-class (mm channel)
#path = '/mnt_pool/c3_users/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/MVA-Training/models/MultiRF_1j1b-mm_2023_09_28/'
#vars = ["train_loosejet1_pt", "train_lep1lep2jet1_pt", "train_lep1lep2_m", "train_lep1lep2_dphi", "train_lep1jet1_pt", "train_lep1lep2_dr", "train_lep1_pt", "train_jet1_pt"]
#
## Multi-class (ee channel)
#path = '/mnt_pool/c3_users/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/MVA-Training/models/MultiRF_1j1b-ee_2023_09_28/'
#vars = ["train_loosejet1_pt", "train_lep1lep2jet1_pt", "train_lep1lep2_m", "train_lep1lep2_dphi", "train_lep1jet1_pt", "train_lep1lep2_dr", "train_lep1_pt", "train_jet1_pt"]

model = joblib.load(path + 'MultiRF.joblib')

initial_type = [('float_input', FloatTensorType([None, len(vars)]))]
onx = convert_sklearn(model, initial_types=initial_type)
with open("rf2j1b_full2022.onnx", "wb") as f:
    f.write(onx.SerializeToString())