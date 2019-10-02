#!/usr/bin/env python3

import argparse
import collections
import FlowCal
import matplotlib.pyplot as plt
import numpy as np
import sys
import seaborn as sns
from scipy import stats


def plot_expression_test(fcs_files,
                         channel,
                         density_fraction_fsc_ssc,
                         density_fraction_fsca_fsch,
                         plot_density,
                         histogram_name,
                         ):
    """takes channel, graph labels are from name of each file, one plot
    objectives: for each file, plot and store fsc vs ssc and fsc h vs a plots
    make cumulative apc-a histogram with one dist for each file
    """
    label_gated_data_dict = collections.OrderedDict()
    for fcs_file in fcs_files:
        file_base_name = fcs_file.split('.fcs')[0].split('/')[-1]

        # read in fcs file
        fcs_data_raw = FlowCal.fcs_io.FCSData(fcs_file)
        # transform from digital representation to relative fluorescence units
        fcs_data_rfi = FlowCal.transform.to_rfi(fcs_data_raw)
        # discard events that exceed high or low limit of detectors
        fcs_data_bandpass = FlowCal.gate.high_low(fcs_data_rfi, channels=['FSC-A', 'SSC-A', channel])
        # gate FSC-A and SSC-A such that f fraction of events are retained [0,1]
        fcs_data_fsc_ssc, mask_fsc_ssc, contour_fsc_ssc = FlowCal.gate.density2d(
            fcs_data_bandpass,
            channels=['FSC-A', 'SSC-A'],
            gate_fraction=density_fraction_fsc_ssc,
            full_output=True
        )
        # gate FSC-A and FSC-H
        fcs_data_fsca_fsch, mask_fsca_fsch, contour_fsca_fsch = FlowCal.gate.density2d(
            fcs_data_fsc_ssc,
            channels=['FSC-A', 'FSC-H'],
            gate_fraction=density_fraction_fsca_fsch,
            full_output=True
        )

        if plot_density:
            # plot ungated and gated data and save to file
            plt.subplot(2, 1, 1)
            FlowCal.plot.density2d(
                fcs_data_bandpass,
                channels=['FSC-A', 'SSC-A'],
                mode='scatter',
            )
            for g in contour_fsc_ssc:
                plt.plot(g[:, 0], g[:, 1], color='k', linewidth=1.25)
            plt.subplot(2, 1, 2)
            FlowCal.plot.density2d(
                fcs_data_fsc_ssc,
                channels=['FSC-A', 'FSC-H'],
                mode='scatter'
            )
            for g in contour_fsca_fsch:
                plt.plot(g[:, 0], g[:, 1], color='k', linewidth=1.25)
            plt.savefig('density_gates_{0}.png'.format(file_base_name), dpi=300)
            plt.close()
        label_gated_data_dict[file_base_name] = fcs_data_fsca_fsch[:, channel]

    # FlowCal.plot.hist1d(list(label_gated_data_dict.values()), normed_height=True)
    # plt.legend(label_gated_data_dict.keys(), loc='upper left')
    # plt.savefig('APC_histogram.png', dpi=300)

    # trying seaborn to get a quick and better looking histogram
    sns.set(color_codes=True)
    sns.set_style(style='white')
    test_list = list(label_gated_data_dict.values())[1]
    for base_name, gated_data in label_gated_data_dict.items():
        distplot = sns.kdeplot(np.log10(list(gated_data)), shade=True, label=base_name)
    # distplot = sns.kdeplot(np.log10(list(label_gated_data_dict.values())[1]), shade=True, label='GPR85')
    # distplot = sns.kdeplot(np.log10(list(label_gated_data_dict.values())[2]), shade=True, label='GPR85 miniGs')
    # # distplot = sns.kdeplot(np.log10(list(label_gated_data_dict.values())[3]), shade=True, label='positive control')
    distplot.set_xlim(0,)
    distplot.set_xlabel('Log10 Fluorescence (A.U.)')
    sns.despine(top=True, left=True, right=True)
    fig = distplot.get_figure()
    if histogram_name:
        fig.savefig(histogram_name, format='png', dpi=500)
    else:
        fig.savefig('histogram_{0}.png', format='png', dpi=500)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="script to plot data from fcs files. density gating of FSC-A/SSC-A "
                    "and of FSC-A/FSC-H is performed using density fraction input arguments. "
                    "after density gating, a histogram is generated for the specified channel")
    required = parser.add_argument_group("required")
    required.add_argument('-f', '--fcs_files', nargs='*', required=True, help='fcs files from cytometer')
    parser.add_argument('-c', '--channel', default='FL5-A',
                        help='channel to use when plotting final histogram')
    parser.add_argument('-d1', '--density_gate_fsc_ssc', type=float, default=0.6,
                        help='fraction of population to keep when gating FSC-A and SSC-A')
    parser.add_argument('-d2', '--density_gate_fsca_fsch', type=float, default=0.9,
                        help='fraction of population to keep when gating FSC-A and FSC-H')
    parser.add_argument('-p', '--plot_density', action='store_true',
                        help='plot 2d dot plots of density gates for each fcs file')
    # parser.add_argument('-dpn', '--density_plot_name', default='density_gates.png',
    #                     help='filename for density plot')
    parser.add_argument('-hn', '--histogram_name',
                        help='name of histogram. default is name of channel')
    args = parser.parse_args()
    plot_expression_test(args.fcs_files,
                         args.channel,
                         args.density_gate_fsc_ssc,
                         args.density_gate_fsca_fsch,
                         args.plot_density,
                         args.histogram_name)
