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

import ORGANIZATION
import PathFinding
import PathOffset

#______________Add include directory__________________
current_dir = os.getcwd()
include_dir = os.path.join(current_dir, '../include/pyth')
sys.path.insert(0, include_dir)
#_____________________________________________________

pas = input("Which mass replay pass are you looking at?")
config = input("Which kinematic configuration are you looking at? 2, 3, or 4?")

#Check if I need to make these config files.

W2min=CONFIG.Function_JSON("W2min",f"../config/cuts{config}.cfg")
W2max=CONFIG.Function_JSON("W2max",f"../config/cuts{config}.cfg")
dxmin=CONFIG.Function_JSON("dxmin",f"../config/cuts{config}.cfg")
dxmax=CONFIG.Function_JSON("dxmax",f"../config/cuts{config}.cfg")
dymin=CONFIG.Function_JSON("dymin",f"../config/cuts{config}.cfg")
dymax=CONFIG.Function_JSON("dymax",f"../config/cuts{config}.cfg")
dybgmin=CONFIG.Function_JSON("dybgmin",f"../config/cuts{config}.cfg")
dybgmax=CONFIG.Function_JSON("dybgmax",f"../config/cuts{config}.cfg")
coinmin=CONFIG.Function_JSON("coinmin",f"../config/cuts{config}.cfg")
coinmax=CONFIG.Function_JSON("coinmax",f"../config/cuts{config}.cfg")
#dymin=-.984
#dymax=.9
#W2max=1.75
#dybgmin=-1.2
#dybgmax=1.2
#dxmin=-2
#dxmax=1
#print(f'W2min: {W2min}')
#print(f'W2max: {W2max}')
#print(f'dxmin: {dxmin}')
#print(f'dxmax: {dxmax}')
#print(f'dymin: {dymin}')
#print(f'dymax: {dymax}')
#print(f'dybgmin: {dybgmin}')
#print(f'dybgmax: {dybgmax}')
#print(f'Coin Min: {coinmin} Coin Max: {coinmax}')

def main():

    org = ORGANIZATION.ORGANIZATION()
    pathF = PathFinding.PathFinding()
    pathO = PathOffset.PathOffset()

    # Load in and declare HCal information

    outputFile = org.Load_HCal(config)

    dataToLoad = np.load(outputFile)
    cblktime=dataToLoad["cblktime_array"]
    cblkatime=dataToLoad["cblkatime_array"]
    nblk_array=dataToLoad["nblk_array"].astype(int)
    cblkid=dataToLoad["cblkid_array"].astype(int)
    nblk=nblk_array[0]

    # Create structure?? comment this out??? Check that it draws the grid right?

    eventCheck=[0,1,2,3,4,5]
    number=eventCheck[0]
    np.round(cblktime[number].astype(float),3)
    np.round(cblkatime[number].astype(float),3)
    np.round(cblke[3],3)
    nblk_array[number]
    cblkid[number][:nblk_array[number]]

    for i in range(0,len(eventCheck)):
        print(nblk_array[i])
        org.draw_grid(cblkid[i]-1,nblk_array[i])

    org.draw_combined_grid(cblkid, nblk_array)

    # Get the path stuff

    master_block = 150
    target_block = 113  # Target block for the path calculation

    path = pathF.find_path_to_block(master_block, target_block)
    print(f"Path from block {master_block} to block {target_block}: {path}")

    pathF.draw_grid_with_path(master_block, target_block, path)

    path_variables = pathF.initialize_adjacent_path_variables()

    for path_name, value in list(path_variables.items())[:10]:
        print(f"{path_name}: {value}")

    histograms = pathF.initialize_adjacent_histograms()
    pathF.plot_histogram(histograms, 'h149to148')

    pathF.populate_histograms_multiple_events(cblkid, cblktime, nblk_array, histograms)

    pathF.plot_histogram_with_gaussian(histograms, 'h151to150')

    pathF.plot_histogram_with_double_gaussian(histograms, 'h151to152')

    primary_fit, secondary_fit = pathO.get_histogram_with_gaussian_secondary(histograms, 'h152to151')

    target_block = np.arange(1,289)
    id_offset=[]
    for i in target_block:
        mean=[]
        error=[]
        if i == 150:
            continue

        path = pathF.find_path_to_block(master_block, i)  # Assume you already have the path function
        #print(path)
        #draw_grid_with_path(master_block, i, path)
        histogramPath= ['h' + p for p in path]
        for j in range(0,len(histogramPath)):
            result=pathO.get_histogram_with_gaussian_secondary(histograms, histogramPath[j])
            mean.append(result[0])
            error.append(result[1]/result[2])
        npError=np.array(error)
        offset=np.sum(mean)
        e=np.sqrt(np.sum(npError**2))
        id_offset.append([i,offset,e])

    # put these elsewhere!!

    lw=1
    plt.figure(figsize=(10,5))
    result=np.array(id_offset)
    result=np.transpose(result)
    plt.errorbar(result[0],result[1],result[2],ls='none',color='black',lw=lw*1.5)

    plt.plot(result[0],result[1],'.',color='blue',markersize=8)
    plt.plot(result[0],result[1],color='black',lw=lw/2)

    plt.xlabel('Block ID')
    plt.ylabel('Block 150 time - Block i time (ns) (ADC)')
    plt.title('Travel Calibration Offsets')
    plt.ylim(-2,12)

    # -----

    lw=1
    plt.figure(figsize=(10,5))
    result=np.array(id_offset)
    result=np.transpose(result)
    plt.errorbar(result[0],result[1],result[2],ls='none',color='black',lw=lw*1.5)

    plt.plot(result[0],result[1],'.',color='blue',markersize=8)
    plt.plot(result[0],result[1],color='black',lw=lw/2)

    plt.xlabel('Block ID')
    plt.ylabel('Block 150 time - Block i time (ns) (ADC)')
    plt.title('Travel Calibration Offsets')
    plt.ylim(-2,12)

    first_array=result[0]
    second_array=result[1]

    insert_index = np.searchsorted(first_array, 150)

    # Insert 150 into the first array
    first_array = np.insert(first_array, insert_index, 150)

    # Insert 0 into the second array
    second_array = np.insert(second_array, insert_index, 0)
    new_offset=np.nan_to_num(second_array, nan=0)
    np.save("ADCOffsetForOffset.npy",new_offset)
    reshaped_array = new_offset.reshape(24, 12)

    # Print the array in the desired format
    for row in reshaped_array:
        formatted_row = ' '.join(f"{num:.2f}" for num in row)  # Format each number to 2 decimal places
        print(formatted_row)

    old_offset = """
    -468.602 -480.413 -460.666 -463.216 -599.934 -570.151 -555.733 -580.574 -474.458 -468.996 -448.016 -492.204
    -454.851 -456.442 -458.228 -443.921 -571.557 -559.782 -584.339 -594.457 -450.178 -455.171 -482.058 -502.654
    -452.762 -491.849 -466.624 -446.454 -576.061 -569.666 -592.834 -602.355 -457.19 -457.027 -503.931 -493.858
    -435.905 -430.477 -425.544 -447.108 -554.623 -553.823 -543.9 -543.204 -439.544 -439.814 -460.386 -478.601
    -429.22 -430.361 -435.06 -433.308 -506.704 -539.935 -540.986 -555.92 -439.751 -439.88 -449.729 -448.084
    -461.392 -456.132 -445.911 -452.946 -540.176 -519.325 -560.799 -569.588 -447.04 -459.5 -441.881 -448.997
    -509.104 -494.054 -507.673 -506.049 -602.321 -607.923 -571.288 -594.528 -462.069 -467.995 -483.508 -488.063
    -471.52 -478.381 -475.626 -475.523 -571.478 -575.791 -581.95 -594.263 -489.027 -485.849 -492.143 -471.779
    -489.231 -478.469 -476.787 -467.032 -576.452 -568.905 -599.653 -590.211 -471.311 -466.66 -486.604 -480.432
    -513.739 -485.985 -479.14 -483.749 -565.799 -569.408 -579.32 -579.375 -468.253 -484.836 -490.031 -468.128
    -451.968 -491.658 -470.477 -479.242 -551.069 -567.825 -569.174 -570.887 -469.936 -468.896 -486.663 -496.796
    -470.336 -463.544 -477.125 -479.511 -567.538 -571.207 -586.871 -586.853 -463.292 -466.959 -463.123 -469.813
    -481.866 -476.62 -477.546 -468.525 -598.605 -574.593 -571.796 -566.034 -435.476 -451.238 -442.327 -444.296
    -433.508 -458.145 -411.467 -449.979 -522.575 -552.225 -533.764 -554.731 -449.031 -448.876 -417.484 -439.317
    -465.175 -452.153 -444.439 -462.377 -555.627 -536.788 -555.365 -562.627 -447.107 -499.819 -449.81 -450.057
    -430.738 -447.683 -419.6 -426.271 -552.742 -569.145 -537.887 -547.837 -446.126 -430.31 -453.853 -457.736
    -437.706 -478.238 -460.915 -456.077 -580.014 -543.548 -570.615 -598.805 -450.183 -466.92 -437.854 -482.645
    -616.248 -636.845 -580.48 -612.485 -529.121 -535.711 -698.794 -716.263 -597.231 -620.625 -421.505 -439.451
    -445.957 -457.081 -418.991 -447.909 -701.289 -693.981 -540.879 -538.506 -439.215 -440.964 -594.326 -616.598
    -434.788 -447.829 -432.764 -433.973 -707.958 -708.962 -537.945 -560.822 -443.196 -426.77 -615.626 -603.881
    -455.333 -446.375 85.6939 -443.408 -691.401 -694.388 -553.712 -560.125 -452.685 -447.345 -622.116 -607.324
    -468.932 -478.302 -439.775 -451.642 -551.364 -534.927 -521.95 -573.76 -450.34 -471.031 -455.979 -463.832
    -458.336 -419.851 -429.116 -434.81 -550.646 -536.144 -554.91 -563.769 -435.664 -444.5 -452.173 -449.867
    -603.818 -610.448 -464.781 -472.051 -556.096 -533.099 -699.243 -737.89 -449.623 -465.017 -453.213 -459.424
    """

    # Convert the string data into a list of floats
    old_offset_list = list(map(float, old_offset.split()))

    # Convert the list into a NumPy array
    old_offset = np.array(old_offset_list)

    #subtracted offsets
    final_offset=new_offset+old_offset

    final_reshaped = final_offset.reshape(24, 12)

    # Print the array in the desired format
    for row in final_reshaped:
        formatted_row = ' '.join(f"{num:.2f}" for num in row)  # Format each number to 2 decimal places
        print(formatted_row)

    filterResult = result[1][~np.isnan(result[1])]
    BadPaths=np.where(~np.isnan(result[1])==False)

    cheker=plt.hist(old_offset-final_offset,bins=30)

    org.HCalArraysNew=Load_HCalNew("3")


if __name__ == '__main__':

    main()
