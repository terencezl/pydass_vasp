import os
import re
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .helpers import determine_tag_value, figs_assert, initiate_figs, plot_helper_settings


def get_lobster(bond=0, filepath='COHPCAR.lobster', ISPIN=None, plot=False, xlim=None, ylim=None,
                on_figs=None):
    """
    Get the COHPCAR or COOPCAR, with consideration of spin-polarization.

    Parameters
    ----------
    bond: int
        the bond number in 'COHPCAR.lobster'/'COOPCAR.lobster', counting from 1. Default to 0, meaning the average
    filepath: string
        filepath, default to 'COHPCAR.lobster'
    ISPIN: int
        user specified ISPIN
        If not given, infer from 'OUTCAR'/'INCAR'.
    plot: bool
        whether to plot the data, default to False
    xlim: list
        the range of x-axis, 2 values in a list
    ylim: list
        the range of y-axis, 2 values in a list(, of the spin-combined plot if ISPIN == 2)
    on_figs: list/int
        the current figure numbers to plot to, default to new figures

    Returns
    -------
    a dict, containing
        'data': a pandas dataframe
        'ax': the axes reference
    """
    # get data
    basename = os.path.basename(filepath)
    datatype = 'COHP' if 'COHP' in basename.upper() else 'COOP'
    with open(filepath, 'r') as f:
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
        ISPIN = determine_tag_value('ISPIN', filepath)

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

        return_dict = {'data_spin_up': pd.DataFrame(**{'columns': col_names1, 'data': data1}),
                       'data_spin_down': pd.DataFrame(**{'columns': col_names2, 'data': data2}),
                       }

    if plot:
        # start plotting
        figs_assert(on_figs, ISPIN, 'lobster')

        if ISPIN == 1:
            col_num = bond * 2 + 1
            initiate_figs(on_figs)
            if datatype == 'COHP':
                plt.plot(data[:, 0], -data[:, col_num])
            elif datatype == 'COOP':
                plt.plot(data[:, 0], data[:, col_num])
            ax = plt.gca()
            plot_helper_settings((xlim, ylim), datatype)
            return_dict.update({'ax': ax})

        elif ISPIN == 2:
            col_bond = bond * 2 + 1
            # Plot the combined COHP/COOP
            initiate_figs(on_figs)
            if datatype == 'COHP':
                plt.plot(data1[:, 0], -data1[:, col_bond] - data2[:, col_bond], label='spin up + down')
            elif datatype == 'COOP':
                plt.plot(data1[:, 0], data1[:, col_bond] + data2[:, col_bond], label='spin up + down')
            ax1 = plt.gca()
            plot_helper_settings((xlim, ylim), datatype)
            # Plot the separated COHP/COOP
            initiate_figs(on_figs)
            if datatype == 'COHP':
                plt.plot(data1[:, 0], -data1[:, col_bond], label='spin up')
                plt.plot(data2[:, 0], -data2[:, col_bond], label='spin down')
            elif datatype == 'COOP':
                plt.plot(data1[:, 0], data1[:, col_bond], label='spin up')
                plt.plot(data2[:, 0], data2[:, col_bond], label='spin down')
            ax2 = plt.gca()
            ylim_sp = None
            if ylim:
                ylim_sp = ylim[:]
                ylim_sp[0] /= 2.
                ylim_sp[1] /= 2.
            plot_helper_settings((xlim, ylim_sp), datatype)
            return_dict.update({'ax_spin_combined': ax1, 'ax_spin_separated': ax2})

    return return_dict


def get_integrated_lobster(filepath='ICOHPLIST.lobster', return_total=True):
    """
    Get the ICOHPLIST or ICOOPLIST, with consideration of spin-polarization.

    Parameters
    ----------
    filepath: string
        filepath, default to 'ICOHPLIST.lobster'
    return_total: bool
        whether return the spin-up and spin-down summed dataframe or a dict with
        key 1 meaning spin-up, key 2 meaning spin-down, default to True

    Returns
    ----------
    a pandas dataframe if return_total is True, or a dict if return_total is False
    """
    if 'COHP' in os.path.basename(filepath):
        filetype = 'COHP'
    elif 'COOP' in os.path.basename(filepath):
        filetype = 'COOP'

    linenum_list = [int(i.split(':')[0]) for i in
                    subprocess.getoutput(' '.join(['grep -n', filetype, filepath])).split('\n')]
    ILOBSTERLIST_dict = {}
    if len(linenum_list) == 1:
        is_mag = False
        ILOBSTERLIST_dict[1] = pd.read_table(filepath, sep='\s+', index_col=filetype + '#', usecols=range(5))
    else:
        is_mag = True
        n_interactions = linenum_list[1] - linenum_list[0] - 1
        ILOBSTERLIST_dict[1] = pd.read_table(filepath, nrows=n_interactions, sep='\s+', index_col=filetype + '#',
                                             usecols=range(5))
        ILOBSTERLIST_dict[-1] = pd.read_table(filepath, skiprows=n_interactions + 1, nrows=n_interactions, sep='\s+',
                                              index_col=filetype + '#', usecols=range(5))

    ILOBSTER = ILOBSTERLIST_dict[1] + ILOBSTERLIST_dict[-1] if is_mag else ILOBSTERLIST_dict[1]
    if return_total:
        return ILOBSTER
    else:
        return ILOBSTERLIST_dict


def filter_lobster_by_elements(e1, e2, ILOBSTER):
    """
    Get the filtered ICOHPLIST or ICOOPLIST, by the constituting elements.

    Parameters
    ----------
    e1, e2: string
        element symbols. Order doesn't matter
    ILOBSTER: pandas dataframe of ICOHPLIST or ICOOPLIST

    Returns
    ----------
    a filtered pandas dataframe by elements
    """
    e1_regex = '^' + e1 + '\d+$'
    e2_regex = '^' + e2 + '\d+$'
    mask = ((ILOBSTER['atomMU'].str.contains(e1_regex) & ILOBSTER['atomNU'].str.contains(e2_regex)) | \
            (ILOBSTER['atomNU'].str.contains(e1_regex) & ILOBSTER['atomMU'].str.contains(e2_regex)))
    return ILOBSTER[mask]
