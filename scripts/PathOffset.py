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

class PathOffset:

    # Function to plot a histogram and fit a Gaussian to it
    def get_histogram_with_gaussian(histograms, hist_name):
        if hist_name in histograms:
            datah = histograms[hist_name]

            # Filter the data to only include values in the range [-30, 30]
            data = datah[(datah >= -10) & (datah <= 10)]

            # Fit a Gaussian to the filtered data
            mu, std = norm.fit(data)

            # Generate the fitted Gaussian curve
            x = np.linspace(-10, 10, 100)
            p = norm.pdf(x, mu, std)

            # Plot the fitted Gaussian curve
            #plt.plot(x, p, 'k', linewidth=2, label=f'mean={mu:.2f}, std={std:.2f}')
            #plt.hist(datah, bins=60, alpha=0.75, range=(-30, 30), density=True, label='Data')
            # Add labels and title
            #plt.title(f"Histogram and Gaussian Fit for {hist_name}")
            #plt.xlabel('Time Difference (ns)')
            #plt.ylabel("Normalized")
            #plt.legend()

            # Show the plot
            return [mu,std,np.sqrt(len(data))]
        else:
            print(f"Histogram {hist_name} not found!")

    def get_histogram_with_gaussian_secondary(histograms, hist_name):
        if hist_name in histograms:
            datah = histograms[hist_name]

            # Filter the data to only include values in the range [-10, 10]
            data = datah[(datah >= -10) & (datah <= 10)]

            # First Gaussian fit to the filtered data
            mu, std = norm.fit(data)

            # Perform a secondary fit using data within one standard deviation of the first mean
            data_secondary = data[(data >= mu - std) & (data <= mu + std)]

            # Fit a Gaussian to the secondary filtered data
            mu_secondary, std_secondary = norm.fit(data_secondary)
            return [mu_secondary, std,np.sqrt(len(data))]
        else:
            print(f"Histogram {hist_name} not found!")
            return None

def method():
  print("PathOffset")
