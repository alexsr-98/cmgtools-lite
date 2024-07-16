# -*- coding: utf-8 -*-

import os, sys, argparse

# -- Data handling -- #
import pandas as pd
import numpy as np
import ROOT as r
from root_numpy import tree2array
import math
import pickle
import joblib
from scipy import stats
# -- ML -- #
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.decomposition import PCA
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score 
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import GradientBoostingClassifier
import xgboost as xgb
from tensorflow import keras
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import label_binarize
from sklearn.metrics import ConfusionMatrixDisplay
import seaborn as sns

# -- Plotting -- #
import matplotlib
matplotlib.use('template')
import matplotlib.pyplot as plt


################# Customization #################

# -- Paths -- #
inputPath = "/lustrefs/hdd_pool_dir/nanoAODv12/tw-run3/productions/2023-12-23/{year}/x_mvatrain/"
years = ["2022", "2022PostEE"]
outputPath = "/nfs/fanae/user/asoto/Proyectos/tW-Run3/CMSSW_12_4_12/src/CMGTools/TTHAnalysis/python/plotter/tw-run3/MVA-Training/models/"

# -- Samples -- #
twFiles = ["tw.root", "tbarw.root"] 
ttbarFiles = ["ttbar.root"]
nonworzFiles = ["ttbarsemilep.root", "ttbarsemilepext1.root"]
dyFiles = ["dy_10to50.root", "dy_50.root"]
dataFiles = ["data.root"]

# -- Use Data to plot correlation plots -- #
useData = True

# -- Variables to be used in the training -- #
vars = ["train_nloosejets", "train_jet1_pt", "train_loosejet1_pt", "train_lep1lep2jet1met_m", "train_lep1lep2jet1_c", "train_lep1lep2jet1_pt",
"train_nbloosejets", "train_lep1lep2jet1met_ptOVERhttot", "train_htlepOVERhttot", "train_httot", "train_lep1jet1_dr", "train_lep1lep2jet1met_mt", "train_lep1lep2jet1_e",
"train_lep1jet1_pt", "train_lep1jet1_m", "train_lep2jet1_pt", "train_lep1lep2jet1met_pt"]

# All vars
vars = ['njets', 'nbjets', 'train_nloosejets', 'train_nbloosejets', 'train_jet1_pt', 'train_loosejet1_pt', 'train_lep1lep2jet1_pt', 'train_lep1lep2jet1met_pt', 'train_lep1lep2jet1met_ptOVERhttot', 'train_lep1lep2jet1met_m', 'train_htlepOVERhttot', 'train_lep1lep2jet1_c', 'train_lep1lep2jet1_cscalar', 'train_httot', 'train_jet2_pt', 'train_lep1jet1_dr', 'train_lep12jet12_dr', 'train_lep12jet12met_dr', 'train_lep1lep2jet1met_mt', 'train_lep1lep2jet1_pz', 'train_lep1lep2jet1_m', 'train_lep1lep2jet1_e', 'train_lep1jet1_pt', 'train_lep1jet1_m', 'train_lep1jet2_m', 'train_lep2jet1_pt', 'train_lep2jet1_m', 'train_lepjet11lep2_dr', 'train_lepjet11lep2_dphi', 'train_lepjet11lep2_deta', 'train_lep1jetjet21_dr', 'train_lep1jetjet21_dphi', 'train_lep1jetjet21_deta', 'train_met_pt', 'train_met_phi', 'train_nfwdjet', 'train_nloosefwdjet', 'train_lep1lep2_pt', 'train_lep1lep2_ptsum', 'train_lep1lep2_dr', 'train_lep1lep2_dphi', 'train_lep1lep2_deta', 'train_lep1lep2_m', 'train_lep1_pt', 'train_lep1_eta', 'train_lep1_phi', 'train_lep1_m', 'train_lep2_pt', 'train_lep2_eta', 'train_lep2_phi', 'train_lep2_m', 'train_jet1_eta', 'train_jet1_phi', 'train_jet1_m', 'train_jet2_eta', 'train_jet2_phi', 'train_jet2_m', 'train_loosejet1_eta', 'train_loosejet1_phi', 'train_loosejet1_m', 'train_fwdjet1_pt', 'train_fwdjet1_eta', 'train_fwdjet1_phi', 'train_fwdjet1_m', 'train_fwdloosejet1_pt', 'train_fwdloosejet1_eta', 'train_fwdloosejet1_phi', 'train_fwdloosejet1_m']

# Vars for 1j1b (TOP-21-010)
#vars = ["train_nloosejets", 
#        "train_jet1_pt",
#        "train_loosejet1_pt",
#        "train_lep1lep2jet1met_m",
#        "train_lep1lep2jet1_c",
#        "train_lep1lep2jet1_pt"]
# Vars for 2j1b (TOP-21-010)
#vars_2j1b = ["train_jet2_pt",
#        "train_lep1jet1_dr",
#        "train_lep12jet12_dr"]


# ----- Test different variables
# Vars for 1j1b (until 2023-09-28)
vars_1j1b = ["train_lep1jet1_pt", "train_lep1lep2jet1met_m", "train_lep1lep2jet1_pt", "train_lep1lep2jet1met_ptOVERhttot", "train_htlepOVERhttot", "train_lep1jet1_dr", "train_lep1lep2jet1met_mt", "train_lep1lep2_m", "train_jet1_pt", "train_lep1lep2jet1_c","train_lep1lep2jet1_pz","train_lep1lep2_dr","train_lep1lep2_dphi", "train_met_pt","train_loosejet1_pt", "train_lep1_pt"]

# Vars for 2j1b (until 2023-09-28)
vars_2j1b = ["train_jet2_pt",
        "train_lep1jet1_dr",
        "train_lep12jet12_dr",
        "train_lep1jet1_pt", "train_lep1lep2jet1met_m", "train_lep1lep2jet1_pt", "train_lep1lep2jet1met_ptOVERhttot", "train_htlepOVERhttot", "train_lep1lep2jet1met_mt", "train_lep1lep2_m", "train_jet1_pt", "train_lep1lep2jet1_c","train_lep1lep2jet1_pz","train_lep1lep2_dr","train_lep1lep2_dphi", "train_met_pt"]

# Vars for 1j1b (until v3 note)
vars_1j1b = ["train_loosejet1_pt", "train_lep1lep2jet1_pt", "train_lep1lep2_m", "train_lep1lep2_dphi",
             "train_lep1jet1_pt", "train_lep1lep2_dr", "train_lep1_pt", "train_jet1_pt",]
#             "train_lep1lep2_ptsum", "train_lep1_eta", "train_lep2_eta"]
             #"train_nbloosejets"]
            # "train_lep1lep2jet1met_ptOVERhttot", "train_lep1lep2jet1_cscalar", "train_httot", "train_htlepOVERhttot", "train_lep1lep2jet1met_pt", "train_lep1lep2jet1met_m"]
             #"train_lep1lep2jet1met_ptOVERhttot", "train_lep1lep2jet1met_mt", "train_met_pt", "train_lep1lep2jet1_cscalar", "train_lep2jet1_m", "train_lep1jet1_m", 'train_lepjet11lep2_dr', 'train_lepjet11lep2_dphi', 'train_lepjet11lep2_deta']

# Vars for 2j1b (until v3 note)
vars_2j1b = ["train_lep1lep2_m", "train_lep12jet12_dr", 
             "train_lep1lep2jet1_pt", "train_lep1jet1_dr", "train_lep1lep2_dr", "train_lep2_pt", "train_jet2_pt",]
             #"train_minimax"]


# Vars for 1j1b 
vars_1j1b = ["train_loosejet1_pt", "train_lep1lep2jet1_pt", "train_lep1lep2_m", "train_lep1lep2_dphi", "train_lep1jet1_pt", "train_lep1_pt", "train_jet1_pt",
             "train_lep1lep2jet1_m"]

# Vars for 2j1b 
vars_2j1b = ["train_lep1lep2_m", "train_lep12jet12_dr", "train_lep1lep2jet1_pt", "train_lep1jet1_dr", "train_lep1lep2_dr", "train_lep2_pt", "train_jet2_pt",
             "train_lep1lep2jet1_c"]

# -- Training parameters -- #
test_size = 0.3


################# Customization end #################

def step_decay(epoch):
    '''
    Step decay function
    '''
    initial_lrate = 0.002
    drop = 0.5
    epochs_drop = 3.0 #reduce a la mitad el lrate cada n epocas
    lrate = initial_lrate * math.pow(drop, math.floor((1+epoch)/epochs_drop))
    return lrate




class modelConstructor:
    def __init__(self, modelId, region, n_jobs):
        self.modelId = modelId
        self.region = region
        self.n_jobs = n_jobs
    
    def getModel(self):
        if self.modelId == "RF":
            if self.region == "1j1b":
                return RandomForestClassifier(n_estimators=2000,max_depth=5,class_weight='balanced',oob_score=True,random_state=0, n_jobs = self.n_jobs, verbose = 1)
            elif self.region == "2j1b":
                return RandomForestClassifier(n_estimators=200,max_depth=3,class_weight='balanced',oob_score=True,random_state=0, n_jobs = self.n_jobs, verbose = 1)
            else:
                print("Region not implemented yet")
                return None
        
        elif self.modelId == "BDT":
            ##model = GradientBoostingClassifier(n_estimators=2000, max_depth=4,random_state=0,learning_rate=0.01,subsample=0.5)
            #model=xgb.XGBClassifier(objective="binary:logistic",n_estimators=500,max_depth=2)
            return AdaBoostClassifier(DecisionTreeClassifier(max_depth=3), algorithm="SAMME", n_estimators=200, learning_rate=0.5)
        
        elif self.modelId == "NN":
            #######
            def swish(x, beta = 1):
                '''
                Swish activation function
                ''' 
                return (x * keras.backend.sigmoid(beta * x))
            #######
            # -- Create a neural network with keras -- #  ###TODO: Add a validation set
            keras.utils.get_custom_objects().update({'swish': swish}) # Add the swish activation function to keras
            model = keras.models.Sequential()
            normal_ini=keras.initializers.glorot_normal(seed=None)
            ## Model 1:
            model.add(keras.layers.Dense(128,  input_shape = (len(vars),), activation='swish'))
            model.add(keras.layers.Dense(64, activation='swish'))
            model.add(keras.layers.Dense(32, activation='swish'))
            model.add(keras.layers.Dense(16, activation='swish'))
            model.add(keras.layers.Dense(8, activation='swish'))
            model.add(keras.layers.Dense(2, activation='softmax'))
            model.summary()
            
            adam=keras.optimizers.Adam(lr=0.0001)
            sgd = keras.optimizers.SGD(lr=0.01)
            adamax=keras.optimizers.Adamax(lr=0.002, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0)
            model.compile(loss='binary_crossentropy', optimizer="adam", metrics = ["accuracy"])
            return model
            
        
        elif self.modelId == "MLP":
            return MLPClassifier(hidden_layer_sizes=(200, 50, 200), max_iter=1000, activation='logistic', solver='adam', random_state=0, learning_rate_init=0.001, verbose=True)        
        
        elif self.modelId == "MultiRF":
            if self.region == "1j1b":
                return RandomForestClassifier(n_estimators=2000,max_depth=5,class_weight='balanced',oob_score=True,random_state=0, n_jobs = self.n_jobs, verbose = 1)
            elif self.region == "2j1b":
                return RandomForestClassifier(n_estimators=2000,max_depth=4,class_weight='balanced',oob_score=True,random_state=0, n_jobs = self.n_jobs, verbose = 1)
            else:
                print("Region not implemented yet")
                return None

# Let's create a class to handle the data, it will return a pandas dataframe with the label
class DataHandler:
    def __init__(self, inputPath, years, files, label, cuts, start = 0, stop = -1):
        self.inputPath = inputPath
        self.files = files
        self.years = years
        self.label = label  
        self.cuts = cuts
        self.dataFrame = self.getData()

    def getData(self):
        # Create a list of TChain objects
        chain = r.TChain("Friends")
        for year in self.years:
            for file in self.files:
                chain.Add(self.inputPath.format(year=year) + file)
        
        ttree = chain
        arr = tree2array(ttree, start = start, stop = stop)
        
        # Convert to dataframe
        df_out = pd.DataFrame(arr)
        df_out["label"] = self.label

        # Apply cuts
        df_out = df_out.query(self.cuts)
        return df_out


# Let's create a class to handle the models
class ModelHandler:
    def __init__(self, model, name):
        self.model = model
        self.name = name
        # Set n_jobs to 1 to avoid problems when running in the cluster and verbose to 0 
        try:
            self.model.set_params(n_jobs = 1, verbose = 0)
        except:
            pass


    def save(self, outputPath, asPickle=False, asJobLib=False):
        # Save the model
        if asPickle:
            with open(outputPath + self.name + ".pkl", 'wb') as f:
                pickle.dump(self.model, f)
        elif asJobLib:
            joblib.dump(self.model, outputPath + self.name + ".joblib")
        else:
            self.model.save(outputPath + self.name + ".h5")
    
    def load(self, outputPath):
        # Load the model
        self.model = load_model(outputPath + self.name + ".h5")

# Create class to evaluate the performance of the model, plot the ROC curve and evaluate overtraining
class ModelEvaluator:
    def __init__(self, model, X_train, X_test, y_train, y_test, vars, isMultiClass, region):
        self.model = model
        self.X_train = X_train
        self.X_test = X_test
        self.y_train = y_train
        self.y_test = y_test
        self.vars = vars
        self.isMultiClass = isMultiClass
        self.region = region
        self.y_pred_train = self.model.predict(self.X_train)
        self.y_pred_test = self.model.predict(self.X_test)
        # Predict probabilities (in keras predict and predict_proba are the same)
        self.y_pred_train_prob = self.model.predict_proba(self.X_train)
        self.y_pred_test_prob = self.model.predict_proba(self.X_test)
        if self.isMultiClass:
            if self.region == "1j1b":
                self.classDict = {1: "tW", 0: r"$t\bar{t}$", 2: "DY"}
            elif self.region == "2j1b":
                self.classDict = {1: "tW", 0: r"$t\bar{t}$", 2: "Non-W/Z"}
        else:
            self.classDict = {1: "tW", 0: r"$t\bar{t}$"}
#        # Comment this part when it is not needed
#        # Measure the time it takes to evaluate the model (in seconds) with a for loop over each event
#        print(X_test.values[0])
#        import time
#        start = time.time()
#        for i in range(len(self.X_test)):
#            self.model.predict_proba([self.X_test.values[i]])
#        end = time.time()
#        print("Time to evaluate the model: ", end - start, "s")

    def plotROCforClassProb(self, outputPath):
        # Plot ROC curve 
        if not self.isMultiClass:
            # Compute ROC curve and ROC area for training and test sets
            fpr_train, tpr_train, _ = roc_curve(self.y_train, self.y_pred_train_prob[:,1])
            roc_auc_train = roc_auc_score(self.y_train, self.y_pred_train_prob[:,1])
            fpr_test, tpr_test, _ = roc_curve(self.y_test, self.y_pred_test_prob[:,1])
            roc_auc_test = roc_auc_score(self.y_test, self.y_pred_test_prob[:,1])

            # Plot ROC curve
            plt.figure()
            lw = 2
            plt.plot(fpr_train, tpr_train, color='darkorange', lw=lw, label='ROC curve train (area = %0.4f)' % roc_auc_train)
            plt.plot(fpr_test, tpr_test, color='green', lw=lw, linestyle='--', label='ROC curve test (area = %0.4f)' % roc_auc_test)
            plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title('Receiver operating characteristic')
            plt.legend(loc="lower right")
            plt.savefig(outputPath + "ROC_prob.png")
            plt.savefig(outputPath + "ROC_prob.pdf")
            plt.close()
        else:
            # Compute the One-vs-Rest multiclass ROC for tW vs rest (training and test)
            # Convert the labels to binary
            y_train_binary = label_binarize(self.y_train, classes=[0, 1, 2]) # This conserves the order of the classes
            y_test_binary = label_binarize(self.y_test, classes=[0, 1, 2])
            for i in range(3):
                fpr_train, tpr_train, _ = roc_curve(pd.DataFrame(y_train_binary[:, i]), self.y_pred_train_prob[:, i])
                roc_auc_train = roc_auc_score(pd.DataFrame(y_train_binary[:, i]), self.y_pred_train_prob[:, i])
                fpr_test, tpr_test, _ = roc_curve(pd.DataFrame(y_test_binary[:, i]), self.y_pred_test_prob[:, i])
                roc_auc_test = roc_auc_score(pd.DataFrame(y_test_binary[:, i]), self.y_pred_test_prob[:, i])
                # Plot ROC curve
                plt.figure()
                lw = 2
                plt.plot(fpr_train, tpr_train, color='darkorange', lw=lw, label='ROC curve train (area = %0.4f)' % roc_auc_train)
                plt.plot(fpr_test, tpr_test, color='green', lw=lw, linestyle='--', label='ROC curve test (area = %0.4f)' % roc_auc_test)
                plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('Receiver operating characteristic ({process} vs rest)'.format(process=self.classDict[i]))
                plt.legend(loc="lower right")
                plt.savefig(outputPath + "ROC_prob_{process}.png".format(process=i))
                plt.savefig(outputPath + "ROC_prob_{process}.pdf".format(process=i))
                plt.close()
            
            # Convert the labels to binary
            # Compute the One-vs-One multiclass ROC for tW vs ttbar (training and test)
            from itertools import combinations
            pair_list = list(combinations([1, 0, 2], 2))
            for ix, (label1, label2) in enumerate(pair_list):
                mask_1 = (self.y_train == label1) | (self.y_train == label2)
                mask_2 = (self.y_test == label1) | (self.y_test == label2)
                y_train_binary = label_binarize(self.y_train[mask_1], classes=[label2, label1])
                y_test_binary = label_binarize(self.y_test[mask_2], classes=[label2, label1])
                fpr_train, tpr_train, _ = roc_curve(pd.DataFrame(y_train_binary[:, 0]), self.y_pred_train_prob[mask_1, label1])
                roc_auc_train = roc_auc_score(pd.DataFrame(y_train_binary[:, 0]), self.y_pred_train_prob[mask_1, label1])
                fpr_test, tpr_test, _ = roc_curve(pd.DataFrame(y_test_binary[:, 0]), self.y_pred_test_prob[mask_2, label1])
                roc_auc_test = roc_auc_score(pd.DataFrame(y_test_binary[:, 0]), self.y_pred_test_prob[mask_2, label1])
                # Plot ROC curve
                plt.figure()
                lw = 2
                plt.plot(fpr_train, tpr_train, color='darkorange', lw=lw, label='ROC curve train (area = %0.4f)' % roc_auc_train)
                plt.plot(fpr_test, tpr_test, color='green', lw=lw, linestyle='--', label='ROC curve test (area = %0.4f)' % roc_auc_test)
                plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                if label1 == 0:
                    plt.xlabel('False Positive Rate')
                    plt.ylabel('True Positive Rate')
                    plt.title('Receiver operating characteristic ({process1} vs {process2})'.format(process1=self.classDict[label1], process2=self.classDict[label2]))
                else:
                    plt.xlabel('False Positive Rate')
                    plt.ylabel('True Positive Rate')
                    plt.title('Receiver operating characteristic ({process1} vs {process2})'.format(process1=self.classDict[label1], process2=self.classDict[label2]))
                plt.legend(loc="lower right") 
                plt.savefig(outputPath + "ROC_prob_{process1}_{process2}.png".format(process1=label1, process2=label2))
                plt.savefig(outputPath + "ROC_prob_{process1}_{process2}.pdf".format(process1=label1, process2=label2))
                plt.close()
                

            






    def evaluate(self, outpath):
        # Evaluate the model
        print("Training accuracy: ")
        print(accuracy_score(self.y_train, self.y_pred_train))
        print("Test accuracy: ")
        print(accuracy_score(self.y_test, self.y_pred_test))
        print("Confusion matrix: ")
        print(confusion_matrix(self.y_test, self.y_pred_test))
        print("Normalized confusion matrix: ")
        print(confusion_matrix(self.y_test, self.y_pred_test, normalize = "true"))
        print("Classification report: ")
        print(classification_report(self.y_test, self.y_pred_test))
        # Plot the confusion matrix
        orderedProcesses = []
        for i in range(len(self.classDict)):
            orderedProcesses.append(self.classDict[i])
        plt.figure()
        disp = ConfusionMatrixDisplay.from_estimator(self.model, self.X_test, self.y_test, display_labels=orderedProcesses, cmap=plt.cm.Blues, normalize = "true")
        disp.ax_.set_title("Normalized confusion matrix")
        disp.plot()
        plt.savefig(outputPath + "confusionMatrix.png")
        plt.savefig(outputPath + "confusionMatrix.pdf")
        plt.close()


    def overtraining(self):
        # Evaluate overtraining
        print(("Overtraining: ", accuracy_score(self.y_train, self.y_pred_train) - accuracy_score(self.y_test, self.y_pred_test)))

    def plotVarImportance(self, outputPath):
        # Plot variable importance
        plt.figure()
        # Order the variables by importance
        importance = self.model.feature_importances_
        indices = np.argsort(importance)[::-1]
        # Select 15 most important variables
        indices = indices[:15]
        top_vars = np.array([self.vars[i].split("train_")[1] for i in indices])
        plt.title("Feature importances")
        plt.bar(list(range(len(top_vars))), importance[indices], color="r", align="center")
        plt.xticks(list(range(len(top_vars))), top_vars, rotation=90)
        plt.xlim([-1, len(top_vars)])
        # Increase bottom margin to fit the x-axis labels
        plt.subplots_adjust(bottom=0.3)
        plt.savefig(outputPath + "varImportance.png")
        plt.savefig(outputPath + "varImportance.pdf")
        plt.close()

    def plotOutputModel(self, outputPath):
        if self.isMultiClass:
            processes = [0, 1, 2]
        else:
            processes = [0, 1]
        
        for i in processes:
            y_pred_train_prob_signal = self.y_pred_train_prob[:,i][self.y_train == 1]
            y_pred_test_prob_signal = self.y_pred_test_prob[:,i][self.y_test == 1]
            y_pred_train_prob_background = self.y_pred_train_prob[:,i][self.y_train == 0]
            y_pred_test_prob_background = self.y_pred_test_prob[:,i][self.y_test == 0]
            if self.isMultiClass:
                y_pred_train_prob_background2 = self.y_pred_train_prob[:,i][self.y_train == 2]
                y_pred_test_prob_background2 = self.y_pred_test_prob[:,i][self.y_test == 2]
 
            # Plot the output of the model
            plt.figure(figsize=(8,5))
            plt.rcParams.update({'font.size': 15})
            plt.hist(y_pred_test_prob_signal, 40, range=(0,1), color='orange', edgecolor='orange', lw=2, label='{process} test'.format(process=self.classDict[1]), alpha=0.2, density=True)
            plt.hist(y_pred_test_prob_background, 40, range=(0,1), color='red', edgecolor='red', lw=2, label='{process} test'.format(process=self.classDict[0]), alpha=0.2, density=True)
            if self.isMultiClass:
                plt.hist(y_pred_test_prob_background2, 40, range=(0,1), color='blue', edgecolor='blue', lw=2, label='{process} test'.format(process=self.classDict[2]), alpha=0.2, density=True)
            plt.hist(y_pred_train_prob_signal, 40, range=(0,1), color='brown',histtype='step', lw=2, label='{process} train'.format(process=self.classDict[1]),  density=True)
            plt.hist(y_pred_train_prob_background, 40, range=(0,1), color='green',histtype='step', lw=2, label='{process} train'.format(process=self.classDict[0]), density=True)
            if self.isMultiClass:
                plt.hist(y_pred_train_prob_background2, 40, range=(0,1), color='purple',histtype='step', lw=2, label='{process} train'.format(process=self.classDict[2]), density=True)
#            # Add poisson error bars to the histograms but don't join them with a line, ls='none'
#            plt.errorbar(np.linspace(0,1,40), np.histogram(y_pred_test_prob_background, 40, range=(0,1), density=True)[0], yerr=np.sqrt(np.histogram(y_pred_test_prob_background, 40, range=(0,1), density=True)[0]), color='red', lw=1, ls='none')
#            plt.errorbar(np.linspace(0,1,40), np.histogram(y_pred_test_prob_signal, 40, range=(0,1), density=True)[0], yerr=np.sqrt(np.histogram(y_pred_test_prob_signal, 40, range=(0,1), density=True)[0]), color='orange', lw=1, ls='none')
#            if self.isMultiClass:
#                plt.errorbar(np.linspace(0,1,40), np.histogram(y_pred_test_prob_background2, 40, range=(0,1), density=True)[0], yerr=np.sqrt(np.histogram(y_pred_test_prob_background2, 40, range=(0,1), density=True)[0]), color='blue', lw=1, ls='none')
#            # And the train histograms
#            plt.errorbar(np.linspace(0,1,40), np.histogram(y_pred_train_prob_background, 40, range=(0,1), density=True)[0], yerr=np.sqrt(np.histogram(y_pred_train_prob_background, 40, range=(0,1), density=True)[0]), color='green', lw=1, ls='none')
#            plt.errorbar(np.linspace(0,1,40), np.histogram(y_pred_train_prob_signal, 40, range=(0,1), density=True)[0], yerr=np.sqrt(np.histogram(y_pred_train_prob_signal, 40, range=(0,1), density=True)[0]), color='brown', lw=1, ls='none')
#            if self.isMultiClass:
#                plt.errorbar(np.linspace(0,1,40), np.histogram(y_pred_train_prob_background2, 40, range=(0,1), density=True)[0], yerr=np.sqrt(np.histogram(y_pred_train_prob_background2, 40, range=(0,1), density=True)[0]), color='purple', lw=1, ls='none')
#            # Compute a Chi2 test of compatibility between the train and test distributions
#            chi2, p = stats.chisquare(np.histogram(y_pred_test_prob_background, 40, range=(0,1), density=True)[0], np.histogram(y_pred_train_prob_background, 40, range=(0,1), density=True)[0])
#            # Add the chi2 test result to the plot
#            plt.text(0.5, 0.9, r'$\chi^2$ = {chi2:.2f}, p = {p:.2f}'.format(chi2=chi2, p=p), horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
#            chi2, p = stats.chisquare(np.histogram(y_pred_test_prob_signal, 40, range=(0,1), density=True)[0], np.histogram(y_pred_train_prob_signal, 40, range=(0,1), density=True)[0])
#            plt.text(0.5, 0.8, r'$\chi^2$ = {chi2:.2f}, p = {p:.2f}'.format(chi2=chi2, p=p), horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
#            if self.isMultiClass:
#                chi2, p = stats.chisquare(np.histogram(y_pred_test_prob_background2, 40, range=(0,1), density=True)[0], np.histogram(y_pred_train_prob_background2, 40, range=(0,1), density=True)[0])   
#                plt.text(0.5, 0.7, r'$\chi^2$ = {chi2:.2f}, p = {p:.2f}'.format(chi2=chi2, p=p), horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
            
            plt.xlim([0,1])
            plt.xlabel('Event probability of being classified as {process}'.format(process=self.classDict[i]))
            plt.title(r"$\bf{CMS}$ $\it{Preliminary}$", fontsize=20, loc='left')
            plt.legend(loc="upper right")
            plt.grid(True)
            plt.savefig(outputPath + "outputModel"+ str(i) +".png")
            plt.savefig(outputPath + "outputModel"+ str(i) +".pdf")
            plt.close()
    
    def plotInputVariables(self, outputPath):
        # Do a loop over the variables
        for var in self.vars:
            # Plot the input variables
            plt.figure(figsize=(8,5))
            plt.rcParams.update({'font.size': 15}) 
            plt.hist(self.X_train[var][self.y_train == 1], 40, range=(self.X_train[var].min(), self.X_train[var].max()), color='orange', edgecolor='orange', lw=2, label='{process}'.format(process=self.classDict[1]), alpha=0.2, density=True)
            plt.hist(self.X_train[var][self.y_train == 0], 40, range=(self.X_train[var].min(), self.X_train[var].max()), color='red', edgecolor='red', lw=2, label='{process}'.format(process=self.classDict[0]), alpha=0.2, density=True)
            if self.isMultiClass:
                plt.hist(self.X_train[var][self.y_train == 2], 40, range=(self.X_train[var].min(), self.X_train[var].max()), color='blue', edgecolor='blue', lw=2, label='{process}'.format(process=self.classDict[2]), alpha=0.2, density=True)
            plt.xlabel(var)
            plt.title(r"$\bf{CMS}$ $\it{Preliminary}$", fontsize=20, loc='left')
            plt.legend(loc="upper left")
            plt.grid(True)
            plt.savefig(outputPath + var + ".png")
            plt.savefig(outputPath + var + ".pdf")
            plt.close()
    
    def computeCorrelationMatrix(self, outputPath, data):
        # Compute a matrix with all pairs of variables and compute the correlation coefficient
        # Select only the columns corresponding to the variables you want to include
        df_subset = self.X_train[self.vars]

        # Compute the correlation matrix
        correlation_matrix = df_subset.corr()

        # Plot the correlation matrix as a heatmap
        plt.figure(figsize=(10, 8))
        plt.rcParams.update({'font.size': 15}) 
        plt.title(r"$\bf{CMS}$ $\it{Preliminary}$", fontsize=20, loc='left')
        sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f", linewidths=0.5)
        plt.xticks(rotation=45,fontsize=9)
        plt.yticks(rotation=0,fontsize=9)
        plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
        plt.savefig(outputPath + "CorrMatrix.png")
        plt.savefig(outputPath + "CorrMatrix.pdf")
        plt.close()    

        # Now we do the same for data
        df_subset = data[self.vars]

        # Compute the correlation matrix
        correlation_matrix = df_subset.corr()

        # Plot the correlation matrix as a heatmap
        plt.figure(figsize=(10, 8))
        plt.rcParams.update({'font.size': 15}) 
        plt.title(r"$\bf{CMS}$ $\it{Preliminary}$", fontsize=20, loc='left')
        sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f", linewidths=0.5)
        plt.xticks(rotation=45,fontsize=9)
        plt.yticks(rotation=0,fontsize=9)
        plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
        plt.savefig(outputPath + "CorrMatrixData.png")
        plt.savefig(outputPath + "CorrMatrixData.pdf")
        plt.close()       
    
    def plotTrainingEvents(self, outputPath):
        # Plot a bar histogram with the number of signal and background events in the training and test sample
        plt.figure(figsize=(8,5))
        plt.rcParams.update({'font.size': 15})
        if self.isMultiClass:
            barValues = [len(self.y_train[self.y_train==1]), len(self.y_train[self.y_train==0]), len(self.y_train[self.y_train==2]), len(self.y_test[self.y_test==1]), len(self.y_test[self.y_test==0]), len(self.y_test[self.y_test==2])]
            plt.bar([0,1,2,3,4,5], barValues, color=['orange', 'red', 'blue','orange', 'red', 'blue'], alpha=0.5)
            plt.xticks([0,1,2,3,4,5], ['{process} train'.format(process=self.classDict[1]), '{process} train'.format(process=self.classDict[0]), '{process} train'.format(process=self.classDict[2]), '{process} test'.format(process=self.classDict[1]), '{process} test'.format(process=self.classDict[0]), '{process} test'.format(process=self.classDict[2])])
            # Reduce the text size of xticks
            plt.xticks(fontsize=12)
            # Logy scale
            plt.yscale('log')
            # Add the number of events in each bin on top of the bars in scientific notation with 1 decimal. Placed below right on top of its bar.
            for i in range(6):
                plt.text(x = i, y = 1.02*barValues[i], s = "{:.1e}".format(barValues[i]), ha = 'center', va = 'bottom')
                
                
        else:
            barValues = [len(self.y_train[self.y_train==1]), len(self.y_train[self.y_train==0]), len(self.y_test[self.y_test==1]), len(self.y_test[self.y_test==0])]
            plt.bar([0,1,2,3], barValues, color=['orange', 'red', 'orange', 'red'], alpha=0.5)
            plt.xticks([0,1,2,3], ['{process} train'.format(process=self.classDict[1]), '{process} train'.format(process=self.classDict[0]), '{process} test'.format(process=self.classDict[1]), '{process} test'.format(process=self.classDict[0])])
            # Add the number of events on top of the bars
            for i in range(4):
                plt.text(x = i, y = 1.02*barValues[i], s = "{:.1e}".format(barValues[i]), ha = 'center', va = 'bottom')

        plt.ylabel('Number of events')
        plt.title(r"$\bf{CMS}$ $\it{Preliminary}$", fontsize=20, loc='left')
        plt.grid(True)
        plt.savefig(outputPath + "trainingEvents.png")
        plt.savefig(outputPath + "trainingEvents.pdf")
        plt.close()






if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage = "python3 tW_MVA_trainer.py [options]", description = "Helper for training MVAs.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--model',      '-M', metavar = 'model',      dest = "model",    required = True, default = "RF")
    parser.add_argument('--region',    '-r', metavar = 'region',    dest = "region",  required = False, default = "1j1b")
    parser.add_argument('--outFolder',   '-o', metavar = 'outFolder',   dest = "outFolder", required = False, default = "RF_1j1b")
    parser.add_argument('--nEvents',  '-nE', metavar = 'nEvents',  dest = "nEvents", required = False, default = -1, type = int) # Number of events to train on
    parser.add_argument('--nthreads',  '-j', metavar = 'nthreads',  dest = "nthreads",required = False, default = 0, type = int)

    # Parse the arguments
    args     = parser.parse_args()
    modelId    = args.model
    region   = args.region
    outFolder = args.outFolder
    nEvents  = args.nEvents
    n_jobs   = args.nthreads

    isMultiClass = False
    if modelId in ["MultiRF"]:
        isMultiClass = True

    # -- Possible models -- #
    modelsId = [
                "RF",      # Random Forest
                "BDT",     # Gradient Boosted Decision Tree
                "NN",      # Neural network
                "MLP",     # Multi-layer perceptron
                "MultiRF", # Multi-class RF
                ]
    if modelId not in modelsId:
        print("ERROR: Model not recognized. Please choose one of the following:")
        print(models)
        sys.exit(1)
    
    # -- Possible regions -- #
    regions = ["nojets", "1j1b", "2j1b", "2j2b", "1j1b-mm", "1j1b-ee"]
    if region not in regions:
        print("ERROR: Region not recognized. Please choose one of the following:")
        print(regions)
        sys.exit(1)    
    # Choose variables based on the region
    if region == "nojets":
        vars = ["train_lep1_pt"]
        cuts = "(njets >= 0) and (channel == 1)"
        channel = "em"
    elif region == "1j1b":
        vars = vars_1j1b
        cuts = "(njets == 1) and (nbjets == 1) and (channel == 1)"
        channel = "em"
    elif region == "2j1b":
        vars = vars_2j1b
        cuts = "(njets == 2) and (nbjets == 1) and (channel == 1)"
        channel = "em"
    elif region == "2j2b":
        vars = ["train_lep1_pt"]    
        cuts = "(njets == 2) and (nbjets == 2) and (channel == 1)"
        channel = "em"
    elif region == "1j1b-mm":
        vars = vars_1j1b
        cuts = "(njets == 1) and (nbjets == 1) and (channel == 2) and (abs(train_lep1lep2_m - 91.2) > 15)"
        region = "1j1b"
        channel = "mm"
    elif region == "1j1b-ee":
        vars = vars_1j1b
        cuts = "(njets == 1) and (nbjets == 1) and (channel == 3) and (abs(train_lep1lep2_m - 91.2) > 15)"
        region = "1j1b"
        channel = "ee"
        
    

    # -- Training -- #
    start = 0 # Starting event
    stop = nEvents # Ending event

    # Get the data
    tw = DataHandler(inputPath, years, twFiles, 1, cuts, start, stop)
    ttbar = DataHandler(inputPath, years, ttbarFiles, 0, cuts, start, stop)
    if isMultiClass:
        if region == "1j1b":
            loosenedCuts = cuts.replace("nbjets == 1", "nbjets >= 0") ##### <--- This is a hack to get more DY events
            dy = DataHandler(inputPath, years, dyFiles, 2, loosenedCuts, start, stop)
        if region == "2j1b":
            #loosenedCuts = cuts.replace("channel == 1", "(channel == 1 or channel == 2 or channel == 3)") ##### <--- This is a hack to get more NonW/Z events
            nonworz = DataHandler(inputPath, years, nonworzFiles, 2, cuts, start, stop)
    
    # Use real data to plot correlation matrix
    dataReal = DataHandler(inputPath, years, dataFiles, 100, cuts, start, stop).dataFrame

    # Concatenate the dataframes
    if not isMultiClass:
        df = pd.concat([tw.dataFrame, ttbar.dataFrame])
    else:
        if region == "1j1b":
            df = pd.concat([tw.dataFrame, ttbar.dataFrame, dy.dataFrame])
        if region == "2j1b":
            df = pd.concat([tw.dataFrame, ttbar.dataFrame, nonworz.dataFrame])
    
    # Split the data
    X_train, X_test, y_train, y_test = train_test_split(df[vars], df["label"], test_size=test_size, random_state=0)
    
    # Remove from X_train, X_test the allweights variable but before save it in a separate variable
    #X_train_weights = X_train["allweights"]
    #X_test_weights = X_test["allweights"]
    #X_train = X_train.drop(columns=["allweights"])
    #X_test = X_test.drop(columns=["allweights"])

    if modelId in ["NN", "MLP"]:
        # Scale the data because the NN is sensitive to the scale of the input variables
        scaler = StandardScaler()
        scaler.fit(X_train)
        X_train = scaler.transform(X_train)
        X_test = scaler.transform(X_test)
    

    # Build the model
    model = modelConstructor(modelId, region, n_jobs).getModel()
    # Train the model (keras)
    if modelId in ["NN"]:
        # Keras
        lr_schedule = keras.callbacks.LearningRateScheduler(step_decay)
        histObj = model.fit(X_train, keras.utils.to_categorical(y_train), epochs=10, callbacks=[lr_schedule], batch_size=64)
    else:
        # Apply weights to the samples
        if isMultiClass:
            if region == "1j1b" and channel == "em":
                weights = {0: 1, 1: 1, 2: 0.15} # If something is very unbalanced you can use weights to balance it
            elif region == "2j1b" and channel == "em":
                weights = {0: 1, 1: 1, 2: 0.55}
            else:
                weights = {0: 1, 1: 1, 2: 1}
            sample_weight = np.ones(len(y_train))
            sample_weight[y_train == 0] = weights[0]
            sample_weight[y_train == 1] = weights[1]
            sample_weight[y_train == 2] = weights[2]
        
            model.fit(X_train, y_train, sample_weight=sample_weight)
        # Scikit-learn
        else:
            model.fit(X_train, y_train)

    # Create the output directory
    outputPath = outputPath + outFolder + "/" + region + "/" # Redifine the output path to be the model name
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)
    
    # Evaluate the model
    modelEvaluator = ModelEvaluator(model, X_train, X_test, y_train, y_test, vars, isMultiClass, region)
    modelEvaluator.evaluate(outputPath)
    modelEvaluator.plotROCforClassProb(outputPath)
    modelEvaluator.overtraining()
    modelEvaluator.plotVarImportance(outputPath)
    modelEvaluator.plotOutputModel(outputPath)
    modelEvaluator.plotInputVariables(outputPath)
    modelEvaluator.computeCorrelationMatrix(outputPath, dataReal)
    modelEvaluator.plotTrainingEvents(outputPath)
    
    
    # Save the model
    modelHandler = ModelHandler(model, modelId)
    modelHandler.save(outputPath, asJobLib=True)
    
    # Load the model
    #modelHandler = ModelHandler(model, modelId)
    #modelHandler.load(outputPath)
