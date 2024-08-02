#This script uses uproot to read in a root file and create a pandas dataframe from it. It then makes some simple plots to test the script.
#To execute: "python3 scripts/PlottingPractice.py <input file>"
#Kate Evans - August 2nd 2024

import math
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    print ("USAGE: %s <input file> <output file>"%(sys.argv[0]))
    sys.exit(1)
inFileName = sys.argv[1] #Second user input

class PlottingPractice:

    def __init__(self):
        self.orig = pd.DataFrame()

    def GenNumpyArray(self,filename):

        file = uproot.open(filename)
        T=file["P"]

        cutData = T.arrays(["bb_etot_over_p","bb_gr_size","bb_gr_tot","bb_gr_x","bb_gr_y","bb_ps_e","bb_ps_x","bb_ps_y","bb_sh_e","bb_sh_x","bb_sh_y","bb_tr_n","bb_tr_ph","bb_tr_th","bb_tr_x","bb_tr_y","e_kine_W2"],library="pd")

        self.orig = cutData

    def DrawHistograms(self):

        df = pd.DataFrame()
        df = self.orig

        fig, axs = plt.subplots(2, 2, figsize=(20,20))

        axs[0,0].hist2d(df["bb_sh_e"], df["bb_gr_size"], bins=[100,10], range=[[0.0, 5.0], [0.0, 10.0]])
        axs[0,0].set_title('Shower Energy [GeV] vs GRINCH Cluster Size')

        axs[0,1].hist(df["bb_ps_e"], bins=100)
        axs[0,1].set_title('Preshower Energy [GeV]')

        axs[1,0].hist(df["bb_sh_e"], bins=100)
        axs[1,0].set_title('Shower Energy [GeV]')

        axs[1,1].hist(df["bb_etot_over_p"], bins=100)
        axs[1,1].set_title('Total BBCal Energy over Momentum')

        plt.show()

if __name__=='__main__':

    #inputFile = "../outfiles/GEN2_He3_run2321.root"

    cls = PlottingPractice()

    cls.GenNumpyArray(str(inFileName))
    cls.DrawHistograms()
