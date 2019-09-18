#!/usr/bin/env python3

import collections
import FlowCal
import matplotlib.pyplot as plt
import sys


def plot_expression_test(input_fcs_files, channel='APC'):
    """takes channel, graph labels are from name of each file, one plot
    objectives: for each file, plot and store fsc vs ssc and fsc h vs a plots
    make cumulative apc-a histogram with one dist for each file
    """
    label_gated_data_dict = collections.OrderedDict()
    for fcs_file in input_fcs_files:
        file_base_name = fcs_file.split('.fcs')[0].split('/')[-1]

        # read in fcs file
        fcs_data_raw = FlowCal.fcs_io.FCSData(fcs_file)
        # transform from digital representation to relative fluorescence units
        fcs_data_rfi = FlowCal.transform.to_rfi(fcs_data_raw)
        # discard events that exceed high or low limit of detectors
        fcs_data_bandpass = FlowCal.gate.high_low(fcs_data_rfi, channels=['FSC-A', 'SSC-A'])
        # gate FSC-A and SSC-A such that f fraction of events are retained [0,1]
        fcs_data_fsc_ssc, mask_fsc_ssc, contour_fsc_ssc = FlowCal.gate.density2d(
            fcs_data_bandpass,
            channels=['FSC-A', 'SSC-A'],
            gate_fraction=0.9,
            full_output=True
        )
        # gate FSC-A and FSC-H
        fcs_data_fsca_fsch, mask_fsca_fsch, contour_fsca_fsch = FlowCal.gate.density2d(
            fcs_data_fsc_ssc,
            channels=['FSC-A', 'FSC-H'],
            gate_fraction=0.9,
            full_output=True
        )

        # plot ungated and gated data and save to file
        plt.subplot(2, 1, 1)
        FlowCal.plot.density_and_hist(
            fcs_data_bandpass,
            gate_contour=contour_fsc_ssc,
            density_channels=['FSC-A', 'SSC-A'],
            density_params={'mode': 'scatter'},
            hist_channels=[],
        )
        plt.subplot(2, 1, 2)
        FlowCal.plot.density_and_hist(
            fcs_data_fsc_ssc,
            gate_contour=contour_fsca_fsch,
            density_channels=['FSC-A', 'FSC-H'],
            density_params={'mode': 'scatter'},
            hist_channels=[],
        )
        plt.savefig('{0}.png'.format(file_base_name), dpi=300)
        label_gated_data_dict[file_base_name] = fcs_data_fsca_fsch[:, 'FL5-A']
    FlowCal.plot.hist1d(label_gated_data_dict.values())
    plt.legend(label_gated_data_dict.keys(), loc='upper left')
    plt.savefig('APC_histogram.png', dpi=300)


if __name__ == '__main__':
    input_fcs_files = sys.argv[1:]
    plot_expression_test(input_fcs_files)
