import re
import numpy as np
import matplotlib.pyplot as plt
from .helpers import determine_tag_value, figs_assert, initiate_figs, plot_helper_settings


def get_lobster(input_file='COHPCAR.lobster', ISPIN=None, plot=False, type='COHP', bond=0, xlim=None, ylim=None,
                on_figs=None):
    """
    Plot the -COHP, with consideration of spin-polarization.

    Parameters
    ----------
    input_file: string
        input file name, default to 'COHPCAR.lobster'
    ISPIN: int
        user specified ISPIN
        If not given, infer from OUTCAR/INCAR.
    plot: bool
        whether to plot the data, default to False
    type: str
        COHP or COOP
    bond: int
        the bond number in COHPCAR.lobster/COOPCAR.lobster, counting from 1. Default to 0, meaning the average
    xlim: list
        the range of x-axis, 2 values in a list
    ylim: list
        the range of y-axis, 2 values in a list(, of the spin-combined plot if ISPIN == 2)
    on_figs: list/int
        the current figure numbers to plot to, default to new figures

    Returns
    -------
    a dict, containing
        'data': a dict that has 2D array of data,
            easily to Pandas DataFrame by pd.DataFrame(**returned_dict['data'])
        'ax': the axes reference
    """
    # get data
    with open(input_file, 'r') as f:
        LOBSTERCAR = f.readlines()

    for line_num, line in enumerate(LOBSTERCAR):
        if re.match(r'No\.\d*:.*\(.*\)', line):
            break
    N_headerlines = line_num

    for line_num, line in enumerate(LOBSTERCAR[N_headerlines:]):
        if not re.match(r'No\.\d*:.*\(.*\)', line):
            break
    N_bonds = line_num

    data_start_line = N_headerlines + N_bonds

    for i in range(len(LOBSTERCAR)):
        LOBSTERCAR[i] = LOBSTERCAR[i].split()

    NEDOS = int(LOBSTERCAR[1][2])
    if ISPIN:
        print("Using user specified ISPIN.")
    else:
        ISPIN = determine_tag_value('ISPIN')

    data = np.array(LOBSTERCAR[data_start_line:data_start_line + NEDOS], dtype=float)

    # confluence and data organizing
    if ISPIN == 1:
        col_names = ['E', 'avg', 'avg_integrated']
        for n_bond in range(1, N_bonds + 1):
            col_names.extend(['bond_{0}'.format(n_bond), 'bond_{0}_integrated'.format(n_bond)])
        return_dict = {'data': {'columns': col_names, 'data': data}}
    elif ISPIN == 2:
        data1 = data[:, :(3 + 2 * N_bonds)]
        data2 = data[:, (3 + 2 * N_bonds):(5 + 4 * N_bonds)]
        data2 = np.column_stack((data[:, 0], data2))
        col_names1 = ['E', 'avg_up', 'avg_integrated_up']
        col_names2 = ['E', 'avg_down', 'avg_integrated_down']
        for n_bond in range(1, N_bonds + 1):
            col_names1.extend(['bond_{0}_up'.format(n_bond), 'bond_{0}_integrated_up'.format(n_bond)])
            col_names2.extend(['bond_{0}_down'.format(n_bond), 'bond_{0}_integrated_down'.format(n_bond)])

        return_dict = {'data_spin_up': {'columns': col_names1, 'data': data1},
                       'data_spin_down': {'columns': col_names2, 'data': data2},
                       }

    if plot:
        # start plotting
        figs_assert(on_figs, ISPIN, 'lobster')

        if ISPIN == 1:
            col_num = bond * 2 + 1
            initiate_figs(on_figs)
            if type == 'COHP':
                plt.plot(data[:, 0], -data[:, col_num])
            elif type == 'COOP':
                plt.plot(data[:, 0], data[:, col_num])
            ax = plt.gca()
            plot_helper_settings((xlim, ylim), type)
            return_dict.update({'ax': ax})

        elif ISPIN == 2:
            col_bond = bond * 2 + 1
            # Plot the combined COHP/COOP
            initiate_figs(on_figs)
            if type == 'COHP':
                plt.plot(data1[:, 0], -data1[:, col_bond] - data2[:, col_bond], label='spin up + down')
            elif type == 'COOP':
                plt.plot(data1[:, 0], data1[:, col_bond] + data2[:, col_bond], label='spin up + down')
            ax1 = plt.gca()
            plot_helper_settings((xlim, ylim), type)
            # Plot the separated COHP/COOP
            initiate_figs(on_figs)
            if type == 'COHP':
                plt.plot(data1[:, 0], -data1[:, col_bond], label='spin up')
                plt.plot(data2[:, 0], -data2[:, col_bond], label='spin down')
            elif type == 'COOP':
                plt.plot(data1[:, 0], data1[:, col_bond], label='spin up')
                plt.plot(data2[:, 0], data2[:, col_bond], label='spin down')
            ax2 = plt.gca()
            ylim_sp = None
            if ylim:
                ylim_sp = ylim[:]
                ylim_sp[0] /= 2.
                ylim_sp[1] /= 2.
            plot_helper_settings((xlim, ylim_sp), type)
            return_dict.update({'ax_spin_combined': ax1, 'ax_spin_separated': ax2})

    return return_dict
