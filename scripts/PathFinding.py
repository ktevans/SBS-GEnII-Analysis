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

class PathFinding:

    # Define the grid dimensions
    n_cols = 12
    n_rows = 24

    # Convert a block ID (1-288) to row and column position
    def block_to_position(block_id):
        row = (block_id - 1) // n_cols
        col = (block_id - 1) % n_cols
        return row, col

    # Convert a row and column position to block ID
    def position_to_block(row, col):
        return row * n_cols + col + 1

    # Function to determine path from master block (150) to any other block
    def find_path_to_block(master_block, target_block):
        # Starting and target positions
        master_row, master_col = block_to_position(master_block)
        target_row, target_col = block_to_position(target_block)

        current_block = master_block
        path = []

        # Traverse vertically first (row-wise)
        while master_row != target_row:
            if master_row < target_row:
                next_block = position_to_block(master_row + 1, master_col)  # Move down
                master_row += 1
            else:
                next_block = position_to_block(master_row - 1, master_col)  # Move up
                master_row -= 1
            path.append(f"{current_block}to{next_block}")
            current_block = next_block

        # Traverse horizontally next (column-wise)
        while master_col != target_col:
            if master_col < target_col:
                next_block = position_to_block(master_row, master_col + 1)  # Move right
                master_col += 1
            else:
                next_block = position_to_block(master_row, master_col - 1)  # Move left
                master_col -= 1
            path.append(f"{current_block}to{next_block}")
            current_block = next_block

        return path

    # Function to draw grid and plot path from master block to target block
    def draw_grid_with_path(master_block, target_block, path):
        rows, cols = 24, 12
        grid = np.zeros((rows, cols))

        # Mark the path on the grid
        for connection in path:
            start_block, end_block = map(int, connection.split('to'))

            # Get the row and column of the end block
            row, col = (end_block - 1) // cols, (end_block - 1) % cols
            grid[row, col] = 1  # Mark path cells as green

        # Mark the master block and target block separately
        master_row, master_col = (master_block - 1) // cols, (master_block - 1) % cols
        target_row, target_col = (target_block - 1) // cols, (target_block - 1) % cols

        grid[master_row, master_col] = 2  # Mark the master block as blue
        grid[target_row, target_col] = 3  # Mark the target block as yellow

        # Create a colormap
        cmap = mcolors.ListedColormap(['red', 'green', 'blue', 'yellow'])  # Red: other, Green: path, Blue: master, Yellow: target

        # Plot the grid
        fig, ax = plt.subplots(figsize=(4, 8))
        ax.imshow(grid, cmap=cmap, aspect='auto')

        # Draw the gridlines
        ax.set_xticks(np.arange(-0.5, cols, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, rows, 1), minor=True)
        ax.grid(which="minor", color="black", linestyle='-', linewidth=2)
        ax.tick_params(which="both", bottom=False, left=False, labelbottom=False, labelleft=False)

        # Show the plot
        plt.show()

    def get_adjacent_blocks(block_id):
        row, col = block_to_position(block_id)
        adjacent_blocks = []

        # Left neighbor
        if col > 0:
            adjacent_blocks.append(position_to_block(row, col - 1))

        # Right neighbor
        if col < n_cols - 1:
            adjacent_blocks.append(position_to_block(row, col + 1))

        # Above neighbor
        if row > 0:
            adjacent_blocks.append(position_to_block(row - 1, col))

        # Below neighbor
        if row < n_rows - 1:
            adjacent_blocks.append(position_to_block(row + 1, col))

        return adjacent_blocks

    def initialize_adjacent_path_variables():
        path_values = {}  # Dictionary to store path variables and their values

        # Loop through all blocks in the grid
        for block in range(1, n_cols * n_rows + 1):  # Blocks from 1 to 288
            adjacent_blocks = get_adjacent_blocks(block)  # Get adjacent blocks

            # Create path variables for each adjacent block
            for adj_block in adjacent_blocks:
                path_name = f"{block}to{adj_block}"
                reverse_path_name = f"{adj_block}to{block}"

                # Assign random values between 0 and 2
                path_values[path_name] = random.randint(0, 2)
                path_values[reverse_path_name] = random.randint(0, 2)

        return path_values

    def initialize_adjacent_histograms():
        histograms = {}  # Dictionary to store histograms

        # Loop through all blocks in the grid
        for block in range(1, n_cols * n_rows + 1):  # Blocks from 1 to 288
            adjacent_blocks = get_adjacent_blocks(block)  # Get adjacent blocks

            # Create empty histograms for each adjacent block pair
            for adj_block in adjacent_blocks:
                hist_name = f"h{block}to{adj_block}"
                reverse_hist_name = f"h{adj_block}to{block}"

                # Initialize histograms as empty arrays
                histograms[hist_name] = np.array([])  # Empty array for histogram
                histograms[reverse_hist_name] = np.array([])  # Empty array for reverse path

        return histograms

    # Function to plot a specific histogram by name
    def plot_histogram(histograms, hist_name):
        if hist_name in histograms:
            plt.hist(histograms[hist_name], bins=150,range=(-30,30), alpha=0.75,color='dodgerblue')
            plt.title(f"{hist_name}")
            plt.xlabel('Time Difference')

            plt.show()
        else:
            print(f"Histogram {hist_name} not found!")

    def populate_histograms_multiple_events(cblkid_2d, cblktime_2d, nblk_array, histograms):
        for k in range(0,3):
            print(f"Starting File {k}")
            dataToLoad=np.load(f"../outfiles/HCalArrays/hcal_data{k}.npz")
            # cblktime=dataToLoad["cblktime_array"]
            cblkatime=dataToLoad["cblkatime_array"]
            nblk_array=dataToLoad["nblk_array"].astype(int)
            cblkid=dataToLoad["cblkid_array"].astype(int)

            cblkid_2d=cblkid
            cblktime_2d=cblkatime
            nbar_array=nblk_array

            num_events = len(cblkid_2d)
            num_events = len(cblkid_2d)  # Number of events

            # Loop over each event
            for event_idx in range(num_events):

                cblkid = cblkid_2d[event_idx]# Block IDs for this event
                cblktime = cblktime_2d[event_idx]  # Block times for this event
                nblk = nblk_array[event_idx]
                if nblk>16:
                    continue
                # Number of valid blocks in this event
                #print(nblk)
                # Populate histograms for this event
                for i in range(nblk):
                    for j in range(i + 1, nblk):  # Compare with the next block in the list

                        # Get the block IDs
                        blk1 = cblkid[i]
                        blk2 = cblkid[j]

                        # Ensure blk1 and blk2 are adjacent (either horizontally or vertically)
                        row1, col1 = blk1 // n_cols, blk1 % n_cols
                        row2, col2 = blk2 // n_cols, blk2 % n_cols

                        # Check if blk1 and blk2 are adjacent (horizontally or vertically)
                        if (abs(row1 - row2) == 1 and col1 == col2) or (row1 == row2 and abs(col1 - col2) == 1):

                            # Get the time difference between the blocks
                            time_diff = cblktime[i] - cblktime[j]
                            if abs(cblktime[i])>500 or abs(cblktime[j]>500):
                                continue
                            if time_diff==0:
                                print(event_idx,cblktime[i],cblktime[j])

                            # Populate histograms with time difference
                            histograms[f"h{blk1}to{blk2}"] = np.append(histograms[f"h{blk1}to{blk2}"], time_diff)
                            histograms[f"h{blk2}to{blk1}"] = np.append(histograms[f"h{blk2}to{blk1}"], -time_diff)

    def plot_histogram_with_gaussian(histograms, hist_name):
        if hist_name in histograms:
            datah = histograms[hist_name]
            plt.figure(figsize=(8,4))
            # Filter the data to only include values in the range [-30, 30]
            data = datah[(datah >= -10) & (datah <= 10)]

            # First Gaussian fit to the initial data
            mu, std = norm.fit(data)

            # Plot the histogram of the original data with counts on the y-axis
            counts, bins, _ = plt.hist(datah, bins=120, alpha=0.75,color='dodgerblue', range=(-30, 30), label='Data')
            bin_width = bins[1] - bins[0]

            # Generate the fitted Gaussian curve for the first fit, scaled to match histogram counts
            x = np.linspace(-10, 10, 100)
            p = norm.pdf(x, mu, std) * len(data) * bin_width  # Scale to match the histogram within the [-10, 10] range

            # Plot the fitted Gaussian curve for the first fit
            plt.plot(x, p, 'k', linewidth=2, label=f'First fit: mean={mu:.2f}, std={std:.2f}')

            # Perform a secondary fit using data within one standard deviation of the first mean
            data_secondary = data[(data >= mu - std) & (data <= mu + std)]
            if len(data_secondary) > 0:  # Make sure there's enough data for a secondary fit
                mu_secondary, std_secondary = norm.fit(data_secondary)

                # Generate the fitted Gaussian curve for the secondary fit, scaled to match histogram counts in the relevant range
                x_secondary = np.linspace(mu - std, mu + std, 100)
                p_secondary = norm.pdf(x_secondary, mu_secondary, std_secondary) * len(data_secondary) * bin_width  # Scale to match histogram counts in the secondary range

                # Plot the fitted Gaussian curve for the secondary fit
                plt.plot(x_secondary, p_secondary, 'r', linewidth=2, label=f'Second fit: mean={mu_secondary:.2f}, std={std_secondary:.2f}')

            # Add labels and title
            plt.title(f"Histogram and Gaussian Fits for {hist_name}")
            plt.xlabel('Time Difference (ns)')
            plt.ylabel("Counts")
            plt.legend(loc='upper left')

            # Show the plot
            plt.show()
            return [mu, std], [mu_secondary, std_secondary] if len(data_secondary) > 0 else None
        else:
            print(f"Histogram {hist_name} not found!")

    def double_gaussian(x, mu1, sigma1, A1, mu2, sigma2, A2):
        """Double Gaussian function: sum of two Gaussian curves."""
        return (A1 * norm.pdf(x, mu1, sigma1)) + (A2 * norm.pdf(x, mu2, sigma2))

    def find_peak_of_double_gaussian(mu1, sigma1, A1, mu2, sigma2, A2):
        # Create a fine grid of x-values within the relevant range (adjust this range as needed)
        x = np.linspace(-15, 15, 1000)

        # Evaluate the double Gaussian function over this grid
        y = double_gaussian(x, mu1, sigma1, A1, mu2, sigma2, A2)

        # Find the x-value where the function reaches its maximum
        peak_x = x[np.argmax(y)]
        peak_y = np.max(y)

        return peak_x, peak_y

    def calculate_fwhm(x, y, peak_y):
        half_max = peak_y / 2  # Half of the maximum value
        # Find indices where the Gaussian curve crosses the half-max value
        indices = np.where(y >= half_max)[0]
        if len(indices) < 2:
            return None, None, None, None  # In case FWHM cannot be determined

        # The left and right boundaries for the FWHM
        left_idx = indices[0]
        right_idx = indices[-1]

        # The x-values corresponding to the FWHM
        left_fwhm = x[left_idx]
        right_fwhm = x[right_idx]

        fwhm_value = right_fwhm - left_fwhm

        return left_fwhm, right_fwhm, fwhm_value, half_max

    def plot_histogram_with_double_gaussian(histograms, hist_name):
        if hist_name in histograms:
            datah = histograms[hist_name]
            plt.figure(figsize=(10, 4))
            # Filter the data to only include values in the range [-10, 10]
            data = datah[(datah >= -10) & (datah <= 10)]
            nEntries = len(data)
            # Plot the histogram of the original data with counts on the y-axis
            counts, bins, _ = plt.hist(datah, bins=60, alpha=0.75, color='dodgerblue', range=(-15, 15), label='Data')
            bin_width = bins[1] - bins[0]

            # Use bin centers for fitting
            bin_centers = (bins[:-1] + bins[1:]) / 2

            # Initial guesses for the parameters of the double Gaussian
            mu1_initial, sigma1_initial = norm.fit(data)
            mu2_initial, sigma2_initial = mu1_initial + 1, sigma1_initial
            A1_initial, A2_initial = max(counts) * bin_width, max(counts) * bin_width / 2  # Scale amplitudes to match the histogram

            initial_guess = [mu1_initial, sigma1_initial, A1_initial, mu2_initial, sigma2_initial, A2_initial]

            # Fit the double Gaussian to the data
            popt, _ = curve_fit(double_gaussian, bin_centers, counts, p0=initial_guess)
            mu1, sigma1, A1, mu2, sigma2, A2 = popt
            peak_x, peak_y = find_peak_of_double_gaussian(mu1, sigma1, A1, mu2, sigma2, A2)

            x = np.linspace(-10, 10, 100)
            # Generate the fitted double Gaussian curve
            p_double_gaussian = double_gaussian(x, mu1, sigma1, A1, mu2, sigma2, A2)
            #fwhm addition
            left_fwhm, right_fwhm, fwhm_value, half_max = calculate_fwhm(x, p_double_gaussian, peak_y)
            print(f"Peak of the double Gaussian occurs at x = {peak_x:.2f} with height = {peak_y:.2f}")

            # PLOTTING
            plt.plot(x, p_double_gaussian, 'r', linewidth=2, label=f'Double Gaussian Fit: \nMean 1={mu1:.2f}, std1={sigma1:.2f} \nMean 2={mu2:.2f}, std2={sigma2:.2f}')

            if fwhm_value != None:
                plt.hlines(half_max, left_fwhm, right_fwhm, colors='blue', linestyles='-', lw=2, label=f"FWHM = {np.round(fwhm_value, 3)}ns")

            # Add labels and title
            plt.title(f"Double Gaussian Fit for {hist_name}")
            plt.xlabel('Time Difference (ns)')
            plt.ylabel("Counts")
            if fwhm_value != None:
                plt.text(-13, half_max, f"N Entries: {nEntries}")
                plt.text(-13, half_max / 1.5, f"FWHM/sqrt(N): {np.round(fwhm_value / np.sqrt(nEntries), 3)}")
                error = fwhm_value / np.sqrt(nEntries)
            else:
                error = 100
            plt.legend(loc='upper left')

            # MEAN1,2 and X^2-----------------------------------------

            # Step 1: Calculate the number of events under each Gaussian
            N1 = abs(A1) * np.sqrt(2 * np.pi) * sigma1  # Number of events under Gaussian 1
            N2 = abs(A2) * np.sqrt(2 * np.pi) * sigma2  # Number of events under Gaussian 2

            # Step 2: Calculate the accuracy as standard error
            accuracy1 = sigma1 / np.sqrt(N1)  # Standard error (accuracy) for Gaussian 1
            accuracy2 = sigma2 / np.sqrt(N2)  # Standard error (accuracy) for Gaussian 2

            # Step 3: Calculate the combined weighted mean (mean1-2)
            mean1_2 = (mu1 / accuracy1**2 + mu2 / accuracy2**2) / (1 / accuracy1**2 + 1 / accuracy2**2)

            # Step 4: Calculate the uncertainty in the combined mean
            sigma_mean1_2 = 1 / np.sqrt(1 / accuracy1**2 + 1 / accuracy2**2)

            # Step 5: Chi-squared calculation for the combined mean
            chi_squared = ((mu1 - mean1_2)**2 / accuracy1**2) + ((mu2 - mean1_2)**2 / accuracy2**2)

            # Plot the combined mean on the graph
            if A1 > 0 and A2 > 0:
                plt.axvline(mean1_2, color='black', lw=1, label=f"Mean1-2: {np.round(mean1_2, 3)}ns ± {np.round(sigma_mean1_2, 3)}")
            else:
                plt.axvline(peak_x, color='black', lw=1, label=f"Mean1-2: {np.round(peak_x, 3)}ns ± {np.round(error, 3)}")

            plt.show()

            # Print out the combined mean and chi-squared
            print(f"Combined Mean (mean1-2): {mean1_2:.3f}")
            print(f"Uncertainty in Combined Mean: {sigma_mean1_2:.3f}")
            print(f"Chi-Squared/DOF: {chi_squared:.3f}")
            print(f"(A1, A2): ({A1}, {A2})")

            # Optionally, you could return these values if you need to use them later
            if A1 > 0 and A2 > 0:
                return peak_x, error, mean1_2, sigma_mean1_2, chi_squared / 2
            else:
                return peak_x, error, peak_x, error, chi_squared / 2
            # ---------------------------------------------------------

        else:
            print(f"Histogram {hist_name} not found!")

def method():
  print("PathFinding")
