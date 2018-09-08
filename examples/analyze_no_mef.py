#!/usr/bin/python
"""
FlowCal Python API example, without using calibration beads data.

This script is divided in two parts. Part one processes data from five cell
samples, and generates plots of each one.

Part two exemplifies how to use the processed cell sample data with
FlowCal's plotting and statistics modules, in order to produce interesting
plots.

For details about the experiment, samples, and instrument used, please
consult readme.txt.

"""
import os
import os.path
import numpy as np
import matplotlib.pyplot as plt
import FlowCal

###
# Definition of constants
###

# Names of the FCS files containing data from cell samples
samples_filenames = ['FCFiles/Data001.fcs',
                     'FCFiles/Data002.fcs',
                     'FCFiles/Data003.fcs',
                     'FCFiles/Data004.fcs',
                     'FCFiles/Data005.fcs',
                     ]

# IPTG concentration of each cell sample, in micromolar.
iptg = np.array([0, 81, 161, 318, 1000])

# Plots will be generated after gating and transforming cell samples. These
# will be stored in the following folder.
samples_plot_dir = 'plot_samples'

if __name__ == "__main__":

    # Check that plot directory exists, create if it does not.
    if not os.path.exists(samples_plot_dir):
        os.makedirs(samples_plot_dir)

    ###
    # Part 1: Processing cell sample data
    ###
    print("\nProcessing cell samples...")

    # We will use the list ``samples`` to store processed, transformed flow
    # cytometry data of each sample.
    samples = []

    # Iterate over cell sample filenames
    for sample_id, sample_filename in enumerate(samples_filenames):

        # Load flow cytometry data from the corresponding FCS file.
        # ``FlowCal.io.FCSData(filename)`` returns an object that represents
        # flow cytometry data loaded from file ``filename``.
        print("\nLoading file \"{}\"...".format(sample_filename))
        sample = FlowCal.fcs_io.FCSData(sample_filename)
        
        # Data loaded from an FCS file is in "Channel Units", the raw numbers
        # reported from the instrument's detectors. The FCS file also contains
        # information to convert these into Relative Fluorescence Intensity
        # (RFI) values, commonly referred to as arbitrary fluorescence units
        # (a.u.). The function ``FlowCal.transform.to_rfi()`` performs this
        # conversion.
        print("Performing data transformation...")
        sample = FlowCal.transform.to_rfi(sample)

        # Gating

        # Gating is the process of removing measurements of irrelevant
        # particles, while retaining only the population of interest.
        print("Performing gating...")

        # ``FlowCal.gate.start_end()`` removes the first and last few events.
        # Transients in fluidics can make these events slightly different from
        # the rest. This may not be necessary in all instruments.
        sample_gated = FlowCal.gate.start_end(sample,
                                              num_start=250,
                                              num_end=100)

        # ``FlowCal.gate.high_low()`` removes events outside a range specified
        # by a ``low`` and a ``high`` value. If these are not specified (as
        # shown below), the function removes events outside the channel's range
        # of detection.
        # Detectors in a flow cytometer have a finite range of detection. If the
        # fluorescence of a particle is higher than the upper limit of this
        # range, the instrument will incorrectly record it with a value equal to
        # this limit. The same happens for fluorescence values lower than the
        # lower limit of detection. These saturated events should be removed,
        # otherwise statistics may be calculated incorrectly.
        # Note that this might not be necessary with newer instruments that
        # record data as floating-point numbers (and in fact it might eliminate
        # negative events). To see the data type stored in your FCS files, run
        # the following instruction: ``print sample_gated.data_type``.
        # We will remove saturated events in the forward/side scatter channels,
        # and in the fluorescence channel FL1.
        sample_gated = FlowCal.gate.high_low(sample_gated,
                                             channels=['FSC','SSC','FL1'])

        # ``FlowCal.gate.density2d()`` preserves only the densest population as
        # seen in a 2D density diagram of two channels. This helps remove
        # particle aggregations and other sparse populations that are not of
        # interest (i.e. debris).
        # We use the forward and side scatter channels, and preserve 50% of the
        # events. Finally, setting ``full_output=True`` instructs the function
        # to return two additional outputs. The last one (``gate_contour``) is
        # a curve surrounding the gated region, which we will use for plotting
        # later.
        sample_gated, __, gate_contour = FlowCal.gate.density2d(
            data=sample_gated,
            channels=['FSC','SSC'],
            gate_fraction=0.5,
            full_output=True)

        # Plot forward/side scatter 2D density plot and 1D fluorescence
        # histograms
        print("Plotting density plot and histogram...")

        # Parameters for the forward/side scatter density plot
        density_params = {}
        # We use the "scatter" mode, in which individual particles will be
        # plotted individually as in a scatter plot, but with a color
        # proportional to the particle density around.
        density_params['mode'] = 'scatter'

        # Parameters for the fluorescence histograms
        hist_params = {}
        hist_params['xlabel'] = 'FL1 Fluorescence (a.u.)'

        # Plot filename
        # The figure can be saved in any format supported by matplotlib (svg,
        # jpg, etc.) by just changing the extension.
        plot_filename = '{}/density_hist_{}.png'.format(
            samples_plot_dir,
            'S{:03}'.format(sample_id + 1))

        # Plot and save
        # The function ``FlowCal.plot.density_and_hist()`` plots a combined
        # figure with a 2D density plot at the top, and an arbitrary number of
        # 1D histograms below. In this case, we will plot the forward/side
        # scatter channels in the density plot, and a histogram of the
        # fluorescence channel FL1 below.
        # Note that we are providing data both before (``sample``) and after
        # (``sample_gated``) gating. The 1D histogram will display the ungated
        # dataset with transparency, and the gated dataset in front with a solid
        # solid color. In addition, we are providing ``gate_contour`` from the
        # density gating step, which will be displayed in the density diagram.
        # This will result in a convenient representation of the data both
        # before and after gating.
        FlowCal.plot.density_and_hist(
            sample,
            sample_gated,
            density_channels=['FSC','SSC'],
            hist_channels=['FL1'],
            gate_contour=gate_contour,
            density_params=density_params,
            hist_params=hist_params,
            savefig=plot_filename)

        # Save cell sample object
        samples.append(sample_gated)

    ###
    # Part 3: Examples on how to use processed cell sample data
    ###

    # Histogram of all samples
    # Here, we plot the fluorescence histograms of all five samples in the same
    # figure, using ``FlowCal.plot.hist1d``. Note how this function can be used
    # in the context of accessory matplotlib functions to modify the axes
    # limits and labels and add a legend, among others.
    plt.figure(figsize=(6,3.5))
    FlowCal.plot.hist1d(samples,
                        channel='FL1',
                        histtype='step',
                        bins=128)
    plt.ylim([0, 2000])
    plt.xlabel('FL1 Fluorescence (a.u.)')
    plt.legend(['{} $\mu M$ IPTG'.format(i) for i in iptg],
               loc='upper left',
               fontsize='small')
    plt.tight_layout()
    plt.savefig('histograms.png', dpi=200)
    plt.close()

    # Here we illustrate how to obtain statistics from the fluorescence of each
    # sample, and how to use them in a plot.
    # The stats module contains functions to calculate different statistics
    # such as mean, median, and standard deviation. Here, we calculate the
    # geometric mean from channel FL1 of each sample, and plot them against the
    # corresponding IPTG concentrations.
    samples_fluorescence = [FlowCal.stats.gmean(s, channels='FL1')
                            for s in samples]
    plt.figure(figsize=(5.5, 3.5))
    plt.plot(iptg,
             samples_fluorescence,
             marker='o',
             color=(0, 0.4, 0.7))
    plt.xlabel('IPTG Concentration ($\mu M$)')
    plt.ylabel('FL1 Fluorescence (a.u.)')
    plt.tight_layout()
    plt.savefig('dose_response.png', dpi=200)
    plt.close()

    print("\nDone.")
