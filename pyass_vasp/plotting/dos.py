import re
import numpy as np
import matplotlib.pyplot as plt
from helpers import determine_tag_value, plot_helper_figs_assert, plot_helper_figs, plot_helper_settings, plot_helper_post

def plot_tdos(axis_range=None, DOSCAR='DOSCAR', ISPIN=None, display=True,
    on_figs=None, save_figs=False, save_data=False, output_prefix='TDOS', return_states_at_Ef=False):
    """
    Plot the total density of states, with consideration of spin-polarization.

    Parameters
    ----------
    axis_range: list
        the range of axes x and y, 4 values in a list
    DOSCAR: string
        DOSCAR file name, default to 'DOSCAR'
    ISPIN: int
        user specified ISPIN. If not given, infer from OUTCAR, or INCAR
    display: bool
        display figures or not
    on_figs: list
        the current figure numbers to plot to, default to new figures
    save_figs: bool
        save figures or not
    save_data: bool
        save data or not
    output_prefix: string
        prefix string before the output files, default to 'TDOS'
    return_states_at_Ef: bool
        calculate the TDOS at Ef or not

    Returns
    -------
    a dict, containing 'data', a 2D numpy array of data from DOSCAR,
        and 'columns', a list of column labels
    """

    with open(DOSCAR, 'r') as f:
        DOSCAR = f.readlines()
    for i in range(len(DOSCAR)):
        DOSCAR[i] = DOSCAR[i].split()

    N_steps = int(DOSCAR[5][2])
    Ef = float(DOSCAR[5][3])

    if ISPIN:
        print "Using user specified ISPIN."
    else:
        ISPIN = determine_tag_value('ISPIN')

    plot_helper_figs_assert(on_figs, ISPIN, 'tdos')

    if ISPIN == 2:
        col_names = ['E', 'total_up', 'total_down', 'integrated_up', 'integrated_down']
        DOS_data = np.array(DOSCAR[6:6+N_steps], dtype=float)
        DOS_data[:, 0] -= Ef

        # Plot the separated TDOS
        plot_helper_figs(on_figs)
        plt.plot(DOS_data[:, 0], DOS_data[:, 1], label='spin up')
        plt.plot(DOS_data[:, 0], DOS_data[:, 2], label='spin down')
        plot_helper_settings(axis_range, 'tdos')
        if axis_range:
            plt.axis([axis_range[0], axis_range[1], axis_range[2], axis_range[3] / 2.])
        if save_figs:
            plt.savefig(output_prefix + '-spin-separated.pdf')
        plot_helper_post(display)

        # Plot the combined TDOS
        plot_helper_figs(on_figs)
        plt.plot(DOS_data[:, 0], DOS_data[:, 1] + DOS_data[:, 2], label='spin up + down')
        plot_helper_settings(axis_range, 'tdos')
        if save_figs:
            plt.savefig(output_prefix + '-spin-combined.pdf')
        plot_helper_post(display)

        if return_states_at_Ef:
            half_window = 0.2
            TDOS_at_Ef = (DOS_data[abs(DOS_data[:, 0] - half_window).argmin(), 3] - \
                    DOS_data[abs(DOS_data[:, 0] + half_window).argmin(), 3] + \
                    DOS_data[abs(DOS_data[:, 0] - half_window).argmin(), 4] - \
                    DOS_data[abs(DOS_data[:, 0] + half_window).argmin(), 4]) / half_window / 2.


    elif ISPIN == 1:
        col_names = ['E', 'total', 'integrated']
        DOS_data = np.array(DOSCAR[6:6+N_steps], dtype=float)
        DOS_data[:, 0] -= Ef

        plot_helper_figs(on_figs)
        plt.plot(DOS_data[:, 0], DOS_data[:, 1])
        plot_helper_settings(axis_range, 'tdos')
        if save_figs:
            plt.savefig(output_prefix + '.pdf')
        plot_helper_post(display)

        if return_states_at_Ef:
            half_window = 0.2
            TDOS_at_Ef = (DOS_data[abs(DOS_data[:, 0] - 0.2).argmin(), 2] - \
                       DOS_data[abs(DOS_data[:, 0] + 0.2).argmin(), 2]) / half_window / 2.

    if save_data:
        np.savetxt(output_prefix + '.txt', DOS_data, '%15.6E', header=' '.join(col_names))

    if return_states_at_Ef:
        return {'dos_data': {'columns': col_names, 'data': DOS_data}, 'TDOS_at_Ef': TDOS_at_Ef}
    else:
        return {'dos_data':{'columns': col_names, 'data': DOS_data}}


def plot_ldos(atom, axis_range=None, DOSCAR='DOSCAR', ISPIN=None, LORBIT=None, display=True,
    on_figs=None, save_figs=False, save_data=False, output_prefix='LDOS'):
    """
    Plot the local projected density of states, with consideration of spin-polarization.

    Parameters
    ----------
    atom: int
        the atom number in DOSCAR/POSCAR interested, counting from 1
    axis_range: list
        the range of axes x and y, 4 values in a list
    DOSCAR: string
        DOSCAR file name, default to 'DOSCAR'
    ISPIN: int
        user specified ISPIN. If not given, infer from OUTCAR, or INCAR
    LORBIT: int
        user specified LORBIT. If not given, infer from OUTCAR, or INCAR
    display: bool
        display figures or not
    on_figs: list
        the current figure numbers to plot to, default to new figures
    save_figs: bool
        save figures or not
    save_data: bool
        save data or not
    output_prefix: string
        prefix string before the output files, default to 'LDOS'

    Returns
    -------
    a dict, containing 'data', a 2D numpy array of data from DOSCAR,
        and 'columns', a list of column labels
    """

    with open(DOSCAR, 'r') as f:
        DOSCAR = f.readlines()
    for i in range(len(DOSCAR)):
        DOSCAR[i] = DOSCAR[i].split()

    N_steps = int(DOSCAR[5][2])
    Ef = float(DOSCAR[5][3])

    if ISPIN:
        print "Using user specified ISPIN."
    else:
        ISPIN = determine_tag_value('ISPIN')

    if LORBIT:
        print "Using user specified LORBIT."
    else:
        LORBIT = determine_tag_value('LORBIT')

    plot_helper_figs_assert(on_figs, ISPIN, 'ldos')

    if ISPIN == 2 and (LORBIT == 11 or LORBIT == 1):
        col_names = ['E', 's_up', 's_down', 'p_y_up', 'p_y_down', 'p_z_up', 'p_z_down', 'p_x_up', 'p_x_down',
                     'd_xy_up', 'd_xy_down', 'd_yz_up', 'd_yz_down', 'd_z2_up', 'd_z2_down',
                     'd_xz_up', 'd_xz_down', 'd_x2y2_up', 'd_x2y2_down']
        DOS_data = np.array(DOSCAR[(6 + (N_steps + 1) * atom):(6 + (N_steps + 1) * atom + N_steps)], dtype=float)
        DOS_data[:, 0] -= Ef

        # Spin up
        plot_helper_figs(on_figs)
        for i in range(1, 18, 2):
            plt.plot(DOS_data[:, 0], DOS_data[:, i], label=col_names[i])

        plot_helper_settings(axis_range, 'ldos')
        if axis_range:
            plt.axis([axis_range[0], axis_range[1], axis_range[2]/2., axis_range[3]/2.])
        if save_figs:
            plt.savefig(output_prefix + '-spin-up.pdf')
        plot_helper_post(display)

        # Spin down
        plot_helper_figs(on_figs)
        for i in range(1, 18, 2):
            plt.plot(DOS_data[:, 0], DOS_data[:, i + 1], label=col_names[i + 1])

        plot_helper_settings(axis_range, 'ldos')
        if axis_range:
            plt.axis([axis_range[0], axis_range[1], axis_range[2]/2., axis_range[3]/2.])
        if save_figs:
            plt.savefig(output_prefix + '-spin-down.pdf')
        plot_helper_post(display)

        # Spin up + down
        plot_helper_figs(on_figs)
        for i in range(1, 18, 2):
            plt.plot(DOS_data[:, 0], DOS_data[:, i] + DOS_data[:, i + 1], label=col_names[i] + '+' + col_names[i + 1])

        plot_helper_settings(axis_range, 'ldos')
        if save_figs:
            plt.savefig(output_prefix + '-spin-combined.pdf')
        plot_helper_post(display)

    elif ISPIN == 2 and (LORBIT == 10 or LORBIT == 0):
        col_names = ['E', 's_up', 's_down', 'p_up', 'p_down', 'd_up', 'd_down']
        DOS_data = np.array(DOSCAR[(6 + (N_steps + 1) * atom):(6 + (N_steps + 1) * atom + N_steps)], dtype=float)
        DOS_data[:, 0] -= Ef

        # Spin up
        plot_helper_figs(on_figs)
        for i in range(1, 6, 2):
            plt.plot(DOS_data[:, 0], DOS_data[:, i], label=col_names[i])

        plot_helper_settings(axis_range, 'ldos')
        if axis_range:
            plt.axis([axis_range[0], axis_range[1], axis_range[2]/2., axis_range[3]/2.])
        if save_figs:
            plt.savefig(output_prefix + '-spin-up.pdf')
        plot_helper_post(display)

        # Spin down
        plot_helper_figs(on_figs)
        for i in range(1, 6, 2):
            plt.plot(DOS_data[:, 0], DOS_data[:, i + 1], label=col_names[i + 1])

        plot_helper_settings(axis_range, 'ldos')
        if axis_range:
            plt.axis([axis_range[0], axis_range[1], axis_range[2]/2., axis_range[3]/2.])
        if save_figs:
            plt.savefig(output_prefix + '-spin-down.pdf')
        plot_helper_post(display)

        # Spin up + down
        plot_helper_figs(on_figs)
        for i in range(1, 6, 2):
            plt.plot(DOS_data[:, 0], DOS_data[:, i] + DOS_data[:, i + 1], label=col_names[i] + '+' + col_names[i + 1])

        plot_helper_settings(axis_range, 'ldos')
        if save_figs:
            plt.savefig(output_prefix + '-spin-combined.pdf')
        plot_helper_post(display)

    elif ISPIN == 1 and (LORBIT == 11 or LORBIT == 1):
        col_names = ['E', 's', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z2', 'd_xz', 'd_x2y2']
        DOS_data = np.array(DOSCAR[(6 + (N_steps + 1) * atom):(6 + (N_steps + 1) * atom + N_steps)], dtype=float)
        DOS_data[:, 0] -= Ef

        plot_helper_figs(on_figs)
        for i in range(1, 10):
            plt.plot(DOS_data[:, 0], DOS_data[:, i], label=col_names[i])

        plot_helper_settings(axis_range, 'ldos')
        if save_figs:
            plt.savefig(output_prefix + '.pdf')
        plot_helper_post(display)

    elif ISPIN == 1 and (LORBIT == 10 or LORBIT == 0):
        col_names = ['E', 's', 'p', 'd']
        DOS_data = np.array(DOSCAR[(6 + (N_steps + 1) * atom):(6 + (N_steps + 1) * atom + N_steps)], dtype=float)
        DOS_data[:, 0] -= Ef

        plot_helper_figs(on_figs)
        for i in range(1, 4):
            plt.plot(DOS_data[:, 0], DOS_data[:, i], label=col_names[i])

        plot_helper_settings(axis_range, 'ldos')
        if save_figs:
            plt.savefig(output_prefix + '.pdf')
        plot_helper_post(display)

    if save_data:
        np.savetxt(output_prefix + '.txt', DOS_data, '%15.6E', header=' '.join(col_names))

    return {'dos_data': {'columns': col_names, 'data': DOS_data}}


def plot_cohp(bond_to_plot, axis_range=None, COHPCAR='COHPCAR.lobster', ISPIN=None, display=True,
    on_figs=None, save_figs=False, save_data=False, output_prefix='COHP'):
    """
    Plot the -COHP, with consideration of spin-polarization.

    Parameters
    ----------
    bond_to_plot: int
        the bond number in COHPCAR.lobster interested, counting from 1
    axis_range: list
        the range of axes x and y, 4 values in a list
    COHPCAR: string
        COHPCAR file name, default to 'COHPCAR.lobster'
    ISPIN: int
        user specified ISPIN. If not given, infer from OUTCAR, or INCAR
    display: bool
        display figures or not
    on_figs: list
        the current figure numbers to plot to, default to new figures
    save_figs: bool
        save figures or not
    save_data: bool
        save data or not
    output_prefix: string
        prefix string before the output files, default to 'COHP'

    Returns
    -------
    a dict, containing 'data', a 2D numpy array of data from DOSCAR,
        and 'columns', a list of column labels
    """

    if ISPIN:
        print "Using user specified ISPIN."
    else:
        ISPIN = determine_tag_value('ISPIN')

    plot_helper_figs_assert(on_figs, ISPIN, 'cohp')

    with open(COHPCAR, 'r') as f:
        COHPCAR = f.readlines()

    for line_num, line in enumerate(COHPCAR):
        if re.match(r'No\.\d*:.*\(.*\)', line):
            break
    N_headerlines = line_num

    for line_num, line in enumerate(COHPCAR[line_num:]):
        if not re.match(r'No\.\d*:.*\(.*\)', line):
            break
    N_bonds = line_num

    data_start_line = N_headerlines + N_bonds

    for i in range(len(COHPCAR)):
        COHPCAR[i] = COHPCAR[i].split()

    N_steps = int(COHPCAR[1][2])

    if ISPIN == 2:
        col_names = ['E', 'avg_up', 'avg_integrated_up']
        for n_bond in range(1, N_bonds + 1):
            col_names.extend(['No.{0}_up'.format(n_bond), 'No.{0}_integrated_up'.format(n_bond)])
        col_names.extend(['avg_down', 'avg_integrated_down'])
        for n_bond in range(1, N_bonds + 1):
            col_names.extend(['No.{0}_down'.format(n_bond), 'No.{0}_integrated_down'.format(n_bond)])
        COHP_data = np.array(COHPCAR[data_start_line:data_start_line + N_steps], dtype=float)

        col_up_to_plot = bond_to_plot * 2 + 1
        col_down_to_plot = (bond_to_plot + N_bonds + 1) * 2 + 1
        # Plot the separated COHP
        plot_helper_figs(on_figs)
        plt.plot(COHP_data[:, 0], -COHP_data[:, col_up_to_plot], label=col_names[col_up_to_plot])
        plt.plot(COHP_data[:, 0], -COHP_data[:, col_down_to_plot], label=col_names[col_down_to_plot])
        plot_helper_settings(axis_range, 'cohp')
        if axis_range:
            plt.axis([axis_range[0], axis_range[1], axis_range[2]/2., axis_range[3]/2.])
        if save_figs:
            plt.savefig(output_prefix + '-spin-separated.pdf')
        plot_helper_post(display)

        # Plot the combined COHP
        plot_helper_figs(on_figs)
        plt.plot(COHP_data[:, 0], -COHP_data[:, col_up_to_plot] - COHP_data[:, col_down_to_plot])
        plot_helper_settings(axis_range, 'cohp')
        if save_figs:
            plt.savefig(output_prefix + '-spin-combined.pdf')
        plot_helper_post(display)

    elif ISPIN == 1:
        col_names = ['E', 'avg', 'avg_integrated']
        for n_bond in range(1, N_bonds + 1):
            col_names.extend(['No.{0}'.format(n_bond), 'No.{0}_integrated'.format(n_bond)])

        COHP_data = np.array(COHPCAR[data_start_line:data_start_line + N_steps], dtype=float)

        col_to_plot = bond_to_plot * 2 + 1
        plot_helper_figs(on_figs)
        plt.plot(COHP_data[:, 0], -COHP_data[:, col_to_plot], label=col_names[col_to_plot])
        plot_helper_settings(axis_range, 'cohp')
        if save_figs:
            plt.savefig(output_prefix + '.pdf')
        plot_helper_post(display)

    if save_data:
        np.savetxt(output_prefix + '.txt', COHP_data, '%15.6E', header=' '.join(col_names))

    return {'coph_data': {'columns': col_names, 'data': COHP_data}}