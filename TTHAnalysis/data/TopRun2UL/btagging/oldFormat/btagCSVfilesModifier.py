'''
This code modifies the centrally provided csv files
to match the same desing as the ones used in the reReco campaign.
'''

import os
import sys  
import ROOT as r
import numpy as np

path = "./"
inputcsv = "wp_deepCSV_106XUL18_v2.csv"

outputPath = "./"
outputcsv = "wp_deepCSV_106XUL18_v2_mod.csv"

# Read the input csv file
dict = {}
with open(path+inputcsv, "r") as f:
    lines = f.readlines()
    for line in lines:
        # Save the header as dict key
        if line.startswith("OperatingPoint"):
            header = line.split(",")
            for i in range(len(header)):
                # remove spaces and newlines
                header[i] = header[i].replace(" ", "")
                dict[header[i]] = []
            continue
        # Save the values as dict values
        values = line.split(",")
        for i in range(len(values)):
            # Add the values to the dict using the correct header key
            dict[header[i]].append(values[i].replace(" ", "").replace("\n", ""))
            # If the header is formula enclose the formula in quotes
            if "formula" in header[i]:
                dict[header[i]][-1] = '"'+dict[header[i]][-1]+'"' + "\n"
            # If the header is jetFlavor substitute 4 by 1 and 5 by 0
            if "jetFlavor" in header[i]:
                dict[header[i]][-1] = dict[header[i]][-1].replace("4", "1")
                dict[header[i]][-1] = dict[header[i]][-1].replace("5", "0")
        #dict[values[0]] = values[1:]

# Write the modified csv file
with open(outputPath+outputcsv, "w") as f:
    # Write the header using header variable
    for i in range(len(header)):
        f.write(header[i])
        if i != len(header)-1:
            f.write(",")
    # Write the values using dict
    for lineNumber in range(len(dict[header[0]])):
        for key in header:
            f.write(dict[key][lineNumber])
            if key != header[-1]:
                f.write(",")





