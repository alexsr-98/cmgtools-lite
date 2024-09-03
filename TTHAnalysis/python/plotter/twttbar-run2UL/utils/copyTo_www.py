import os, sys

# This script takes a folder as an input an moves all the files in that folder to the www folder and puts an index file in all the subfolders

www_path = "/nfs/fanae/user/asoto/www/private/bb4lRun2UL/"
indexPath = "tw-run3/utils/index.php"

# Get the folder name from the command line
folder = sys.argv[1]

# Copy the index.php into all the subfolders with plots
for i in os.walk(folder):
    ifolder, isubfolder, ifiles = i
    for file in ifiles:
        if "png" in file or "pdf" in file:
            os.system("cp " + indexPath + " " + ifolder)
            break

# Copy everything in the folder to the www folder
os.system("cp -R " + folder + " " + www_path)

