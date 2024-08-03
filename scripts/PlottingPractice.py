#This script uses uproot to read in a root file and create a pandas dataframe from it. It then makes some simple plots to test the script.
#To execute: "python3 scripts/PlottingPractice.py <input file>"
#Kate Evans - August 2nd 2024

##Import necessary modules
import math
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

##Test the number of user inputs
if len(sys.argv) != 2:
    print ("USAGE: %s <input file> <output file>"%(sys.argv[0]))
    sys.exit(1)
inFileName = sys.argv[1] #Second user input is the root file that we will read

class PlottingPractice:

    def __init__(self):
        self.orig = pd.DataFrame() #Declare a pandas dataframe that we can call throughout the class

    def GenNumpyArray(self,filename):

        file = uproot.open(filename) #Use uproot to open the root file
        T=file["P"] #Define T as the TTree from the root file

        ##Note: my TTree is named "P", but you should edit the line above to match the Tree structure from your root file

        cutData = T.arrays(["bb_etot_over_p","bb_gr_size","bb_gr_tot","bb_gr_x","bb_gr_y","bb_ps_e","bb_ps_x","bb_ps_y","bb_sh_e","bb_sh_x","bb_sh_y","bb_tr_n","bb_tr_ph","bb_tr_th","bb_tr_x","bb_tr_y","e_kine_W2"],library="pd")
        #Create an array of values from the TTree. This allows you to choose what values to include in your script.

        self.orig = cutData #Set self.orig as the array of TTree values so that you can call this array anywhere in the script

    def DrawHistograms(self):

        ##Define a new dataframe and set it equal to self.orig. If you've run GenNumpyArray first, then df will equal cutData.
        df = pd.DataFrame()
        df = self.orig

        ##Define a figure with multiple subplot axes. Here, I am making a 2x2 figure.
        fig, axs = plt.subplots(2, 2, figsize=(20,20))

        ##Use matplotlib commands to plot 1d and 2d histograms of various variables from the TTree array
        
        axs[0,0].hist2d(df["bb_sh_e"], df["bb_gr_size"], bins=[100,10], range=[[0.0, 5.0], [0.0, 10.0]])
        axs[0,0].set_title('Shower Energy [GeV] vs GRINCH Cluster Size')

        axs[0,1].hist(df["bb_ps_e"], bins=100)
        axs[0,1].set_title('Preshower Energy [GeV]')

        axs[1,0].hist(df["bb_sh_e"], bins=100)
        axs[1,0].set_title('Shower Energy [GeV]')

        axs[1,1].hist(df["bb_etot_over_p"], bins=100)
        axs[1,1].set_title('Total BBCal Energy over Momentum')

        plt.show()

        ##Need to add line to save the plot to an output file

if __name__=='__main__':

    #inputFile = "../outfiles/GEN2_He3_run2321.root"

    ##Name your PlottingPractice() class something shorter
    cls = PlottingPractice()

    ##Run the two methods from PlottingPractice()
    cls.GenNumpyArray(str(inFileName))
    cls.DrawHistograms()
