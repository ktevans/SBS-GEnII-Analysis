import ROOT as r
import math
import array
import os
import sys
import random
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm

import CONFIG
import DBPARSE
from UTILITIES import *
from SIMFITS2D import DistributionFits2D
from ROOT import gStyle, TChain, TH1F, TCanvas, TLegend

class ORGANIZATION:

    def list_branches(rootfile_path, tree_name="Tout"):
      # Open the ROOT file
      root_file = r.TFile.Open(rootfile_path)

      # Check if the file was successfully opened
      if not root_file or root_file.IsZombie():
          print(f"Error: Could not open ROOT file '{rootfile_path}'.")
          return

      # Access the TTree
      tree = root_file.Get(tree_name)

      # Check if the tree exists
      if not tree:
          print(f"Error: Tree '{tree_name}' not found in file '{rootfile_path}'.")
          root_file.Close()
          return

      # Get the list of branches
      branch_list = tree.GetListOfBranches()

      # Print all branch names
      print(f"Branches in tree '{tree_name}':")
      for branch in branch_list:
          print(branch.GetName())

      # Close the ROOT file
      root_file.Close()

    def Load_HCal(config):
        rootfile = f"../outfiles/HCal_data_GEN3_sbs100p_nucleon_np_model1.root"
        #rootfile = f"../outfiles/oldhcal.root"
        rootfile2=f"../outfiles/HCal_data_GEN3_sbs100p_nucleon_p_model1.root"
        C = r.TChain("Tout")
        C.Add(rootfile)
        C.SetBranchStatus("*", 0)


        #HCAL VARS_____________________________________
        C.SetBranchStatus("cx", 1)
        C.SetBranchStatus("cy", 1)
        C.SetBranchStatus("cblktime",1)
        C.SetBranchStatus("cblkatime",1)
        C.SetBranchStatus("cblkid",1)
        C.SetBranchStatus("cblke",1)
        C.SetBranchStatus("tdc",1)
        C.SetBranchStatus("atime",1)
        C.SetBranchStatus("bbcal_time",1)
        C.SetBranchStatus("pblkid",1)
        C.SetBranchStatus("nclus",1)
        C.SetBranchStatus("nblk",1)

        bbtime=array.array('d',[0])
        tdc=array.array('d',[0])
        atime=array.array('d',[0])
        pblkid=array.array('d',[0])
        cx = array.array('d', [0]*10)
        cy = array.array('d', [0]*10)
        cblktime=array.array('d',[0]*25)
        cblkatime=array.array('d',[0]*25)
        cblkid=array.array('d',[0]*25)
        cblke=array.array('d',[0]*25)
        nclus=array.array('d',[0])
        nblk=array.array('i',[0])


        C.SetBranchAddress("cx", cx)
        C.SetBranchAddress("cy", cy)
        C.SetBranchAddress("cblktime",cblktime)
        C.SetBranchAddress("cblkatime",cblkatime)
        C.SetBranchAddress("cblkid",cblkid)
        C.SetBranchAddress("cblke",cblke)
        C.SetBranchAddress("tdc",tdc)
        C.SetBranchAddress("atime",atime)
        C.SetBranchAddress("pblkid",pblkid)
        C.SetBranchAddress("bbcal_time",bbtime)
        C.SetBranchAddress("nclus",nclus)
        C.SetBranchAddress("nblk",nblk)

        #_____________________________________


        #CUT VARS_____________________________________
        C.SetBranchStatus("W2",1)
        C.SetBranchStatus("dx",1)
        C.SetBranchStatus("dy",1)

        W2=array.array('d',[0])
        dx=array.array('d',[0])
        dy=array.array('d',[0])

        C.SetBranchAddress("W2",W2)
        C.SetBranchAddress("dx",dx)
        C.SetBranchAddress("dy",dy)

        #_____________________________________

        cx_array=[]
        cy_array=[]
        cblktime_array=[]
        cblkatime_array=[]
        cblkid_array=[]
        cblke_array=[]
        tdc_array=[]
        bbtime_array=[]
        atime_array=[]
        pblkid_array=[]
        target=[]
        nclus_array=[]
        nblk_array=[]
        #tester=[]
        C.GetEntry(0)
        Entries=C.GetEntries()
        passedcut=0
        #Entries=6000000
        chunkValue=2000000
        totalIterations=Entries//chunkValue
        for j in range(0,totalIterations):
            print(f"Beginning Chunk {j} out of {totalIterations}")

            cblktime_array=[]
            cblkatime_array=[]
            cblkid_array=[]
            cblke_array=[]
            nblk_array=[]

            #for i in range(0,Entries):
            for i in range(j*chunkValue,(j+1)*chunkValue):
                C.GetEntry(i)
                #COOLPROGRESSTRACKER___________________________________________________________________________
                if i % 10000 == 0 or i == Entries - 1:
                    progress = f'Processing entry {i + 1}/{Entries} ({(i + 1) / Entries * 100:.2f}%)\r'
                    sys.stdout.write(progress)
                    sys.stdout.flush()
                #______________________________________________________________________________________________

                #tester.append(list(cblkatime))
                #print(tester,'\n')
                #cut----------------
                #wcut=W2min<W2[0]<W2max
                #dxcut=dxmin<dx[0]<dxmax
                #dycut=dymin<dy[0]<dymax
                #cut=wcut and nblk[0]>1
                cut = nblk[0]>1
                #-------------------
                #EnergyOfCluster=np.sum(np.array(cblke[:nblk[0]]))
                #cut=EnergyOfCluster>.1
                if cut:
                    passedcut+=1
                    #cx_array.append(cx)
                    #cy_array.append(cy)
                    #tdc_array.append(tdc[0])
                    #bbtime_array.append(bbtime[0])
                    #atime_array.append(atime[0])

                    #pblkid_array.append(pblkid[0])
                    cblktime_array.append(list(cblktime))
                    cblkatime_array.append(list(cblkatime))
                    #nclus_array.append(nclus[0])
                    nblk_array.append(nblk[0])
                    cblkid_array.append(list(cblkid))
                    #cblke_array.append(list(cblke))
                    #target.append(0)

            outFile = f"../outfiles/HCalArrays/hcal_data{j}.npz"

            np.savez(f"../outfiles/HCalArrays/hcal_data{j}.npz",
                        cblktime_array=cblktime_array,
                        cblkatime_array=cblkatime_array, nblk_array=nblk_array,cblkid_array=cblkid_array)

        print(f"Processing complete. Data saved to 'HCalArrays'")
        return outFile

    def draw_grid(input_array, nblk):
        rows, cols = 24, 12
        grid = np.zeros((rows, cols))

        for idx, value in enumerate(input_array[:nblk]):
            if value >= 0 and value < rows * cols:
                row = value // cols
                col = value % cols
                if idx == 0:
                    grid[row, col] = 2  # Mark the first cell as blue
                else:
                    grid[row, col] = 1  # Mark other cells as green

        fig, ax = plt.subplots(figsize=(4, 8))
        cmap = mcolors.ListedColormap(['red', 'green', 'blue'])
        ax.imshow(grid, cmap=cmap, aspect='auto')
        ax.set_xticks(np.arange(-0.5, cols, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, rows, 1), minor=True)
        ax.grid(which="minor", color="black", linestyle='-', linewidth=2)
        ax.tick_params(which="both", bottom=False, left=False, labelbottom=False, labelleft=False)

        plt.show()

    def draw_combined_grid(cblkid_arrays, nblk_array, rows=24, cols=12):
        grid = np.zeros((rows, cols))  # Initialize the grid as all red (0)

        # Loop through all cblkid arrays and their corresponding nblk values
        for i, cblkid_array in enumerate(cblkid_arrays):
            # Ensure the index exists in nblk_array and limit the range of values accordingly
            #print(i)
            for value in cblkid_array[:nblk_array[i]]:
                #print(value)
                if value >= 0 and value < rows * cols:
                    row = value // cols
                    col = value % cols
                    grid[row, col] = 1  # Mark the block as green

        # Plot the grid
        fig, ax = plt.subplots(figsize=(4, 8))
        cmap = mcolors.ListedColormap(['red', 'green'])
        ax.imshow(grid, cmap=cmap, aspect='auto')
        ax.set_xticks(np.arange(-0.5, cols, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, rows, 1), minor=True)
        ax.grid(which="minor", color="black", linestyle='-', linewidth=2)
        ax.tick_params(which="both", bottom=False, left=False, labelbottom=False, labelleft=False)

        plt.show()

    def Load_HCalNew(config):
        rootfile = f"../outfiles/HCal_data_GEN3_sbs100p_nucleon_np_model1.root"
        #rootfile = f"../outfiles/oldhcal.root"
        rootfile2=f"../outfiles/HCal_data_GEN3_sbs100p_nucleon_p_model1.root"
        C = r.TChain("Tout")
        C.Add(rootfile)
        C.SetBranchStatus("*", 0)


        #HCAL VARS_____________________________________
        #HCAL VARS_____________________________________
        C.SetBranchStatus("cx", 1)
        C.SetBranchStatus("cy", 1)
        C.SetBranchStatus("cblktime",1)
        C.SetBranchStatus("cblkatime",1)
        C.SetBranchStatus("cblkid",1)
        C.SetBranchStatus("cblke",1)
        C.SetBranchStatus("tdc",1)
        C.SetBranchStatus("atime",1)
        C.SetBranchStatus("bbcal_time",1)
        C.SetBranchStatus("pblkid",1)
        C.SetBranchStatus("nclus",1)
        C.SetBranchStatus("nblk",1)

        bbtime=array.array('d',[0])
        tdc=array.array('d',[0])
        atime=array.array('d',[0])
        pblkid=array.array('d',[0])
        cx = array.array('d', [0]*10)
        cy = array.array('d', [0]*10)
        cblktime=array.array('d',[0]*25)
        cblkatime=array.array('d',[0]*25)
        cblkid=array.array('d',[0]*25)
        cblke=array.array('d',[0]*25)
        nclus=array.array('d',[0])
        nblk=array.array('i',[0])


        C.SetBranchAddress("cx", cx)
        C.SetBranchAddress("cy", cy)
        C.SetBranchAddress("cblktime",cblktime)
        C.SetBranchAddress("cblkatime",cblkatime)
        C.SetBranchAddress("cblkid",cblkid)
        C.SetBranchAddress("cblke",cblke)
        C.SetBranchAddress("tdc",tdc)
        C.SetBranchAddress("atime",atime)
        C.SetBranchAddress("pblkid",pblkid)
        C.SetBranchAddress("bbcal_time",bbtime)
        C.SetBranchAddress("nclus",nclus)
        C.SetBranchAddress("nblk",nblk)

        #_____________________________________


        #CUT VARS_____________________________________
        C.SetBranchStatus("W2",1)
        C.SetBranchStatus("dx",1)
        C.SetBranchStatus("dy",1)

        W2=array.array('d',[0])
        dx=array.array('d',[0])
        dy=array.array('d',[0])

        C.SetBranchAddress("W2",W2)
        C.SetBranchAddress("dx",dx)
        C.SetBranchAddress("dy",dy)

        #_____________________________________

        cx_array=[]
        cy_array=[]
        cblktime_array=[]
        cblkatime_array=[]
        cblkid_array=[]
        cblke_array=[]
        tdc_array=[]
        bbtime_array=[]
        atime_array=[]
        pblkid_array=[]
        target=[]
        weighted_time=[]

        C.GetEntry(0)
        Entries=C.GetEntries()
        passedcut=0
        for i in range(0,Entries):

            C.GetEntry(i)

            #cut----------------
            wcut=W2min<W2[0]<W2max
            dxcut=dxmin<dx[0]<dxmax
            dycut=dymin<dy[0]<dymax
            cut=wcut and dxcut and dycut
            cut = wcut and dycut
            #-------------------

            if cut:
                try:

                    # Amplitude weighted time

                            #adjust offsets of blocks affected
                    for i in range(0,len(cblkid[:nblk[0]])):
                        cblktime[i]=cblktime[i]-old_offset[int(cblkid[i])-1]+final_offset[int(cblkid[i])-1]

                    #print(cblke[:nblk[0]])

                    npENERGY=np.array(cblke[:nblk[0]])
                    npTIME=np.array(cblktime[:nblk[0]])
                    weighted_time.append(np.sum(npENERGY*npTIME)/np.sum(npENERGY))







                    # Amplitude weighted time


                    tdc_array.append(tdc[0]-old_offset[int(pblkid[0]-1)]+final_offset[int(pblkid[0]-1)])

                    #atime_array.append(atime[0]-old_offsetADC[int(pblkid[0]-1)]+final_offsetADC[int(pblkid[0]-1)])
                    bbtime_array.append(bbtime[0])
                    passedcut+=1
                    cx_array.append(cx)
                    cy_array.append(cy)

                    pblkid_array.append(pblkid[0])
                    cblktime_array.append(cblktime)
                    cblkatime_array.append(cblkatime)

                    cblkid_array.append(cblkid)
                    cblke_array.append(cblke)
                    target.append(0)
                except:
                    print("Error with event")






        HCalArrays=[cx_array,cy_array,cblktime_array,cblkid_array,cblke_array,tdc_array,pblkid_array,cblkatime_array
                    ,"atime_array",target,bbtime_array,weighted_time]
        return HCalArrays


def method():
  print("ORGANIZATION")
