#This script opens a root file and reads it. The it will fill a simple histogram and write it to a new root file.
#Based on tutorial from 2018 HASCO Summer School: https://indico.cern.ch/event/704163/contributions/2936719/
#Kate Evans - August 2nd 2024

import ROOT
import sys

#First, read in file names as user input. These will be used by executing "python3 practice.py <input file> <output file>". Thus, there should be three arguments read in. We will check that there are exactly 3.
if len(sys.argv) != 3:
    print ("USAGE: %s <input file> <output file>"%(sys.argv[0]))
    sys.exit(1)
inFileName = sys.argv[1] #Second user input
outFileName = sys.argv[2] #Third user input

print ("Reading from", inFileName, "and writing to", outFileName)

#Open the input file and read from it.
inFile = ROOT.TFile.Open(inFileName,"READ")

#After opening the file, we want to get the TTree from it. This is done by using "Get" on the file. Use the proper TTree name from the file, in this case "P".
tree = inFile.Get("P")

h_W = ROOT.TH1D("h_W", ";W", 150, 0.0, 3.0) #1D histogram
h_W.Sumw2() #This configures the way we want to handle statistical uncertainties in the histogram. Most of the time you want uncertainties to be calculated using the sum of the weighted squared, or Sumw2. This is also true in this case -- is it?

#We now have a tree and a histogram. The next step is to actually loop over all of the events in a tree so that we can fill the histogram. Loop over the events in the TTree as follows:
for entryNum in range(0, tree.GetEntries()):
    tree.GetEntry(entryNum)
    #Fill your histogram with the desired variable.
    h_W.Fill(getattr(tree, "e_kine_W2")) #Use the special python function "getattr" (get attribute) to retrive variables from the tree after they have been loaded into memory. This is not possible in C++.

h_W.SetDirectory(0) #Once the histogram has been filled, we want to make sure that it doesn’t disappear. By default, histograms are linked to the last file that was opened, so when you close the file the histogram may disappear.We want to change this behaviour, to say that the histogram we just created has no parent file, and thus should not be removed when any files are closed.

#Now that we have set the histogram to have a global scope, we can safely close the input file. It is good practice to close files when we are done using them to help with memory management and to avoild multiple file access problems.
inFile.Close()

#Now that we have filled the histogram and closed the input file, we want to open a new file to store the results. We will open it using the "RECREATE" mode, which means that we are opening it for writing AND if a file with the specified name alreadu exists, we are ok to overwrite it with our new file.
outHistFile = ROOT.TFile.Open(outFileName, "RECREATE")

#After opening the file, we need to switch our scope to being within that file.
outHistFile.cd()

#Now we can write our histogram to the output file.
h_W.Write()

#Finally, we can close the output file.
outHistFile.Close()
