import re
import warnings
import numpy as np
import matplotlib.pyplot as plt
from helpers import determine_tag_value, figs_assert, initiate_figs, display_or_close_figs
from ..xml_utils import parse


# internal
def plot_helper_settings(axis_range, data_type, save_figs, output):
    plt.axhline(y=0, c='k')
    plt.axvline(x=0, ls='--', c='k', alpha=0.5)
    if axis_range:
        plt.axis([axis_range[0], axis_range[1], axis_range[2], axis_range[3]])
    plt.xlabel('Energy (eV)')
    if data_type == 'tdos':
        plt.ylabel('TDOS (States / Unit Cell / eV)')
    elif data_type == 'ldos':
        plt.ylabel('LDOS (States / Unit Cell / eV)')
    elif data_type == 'cohp':
        plt.ylabel('-pCOHP (Arbitrary Unit / Unit Cell / eV)')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plt.legend(loc=0, fontsize='small')
    ax = plt.gca()
    ax.figure.set_tight_layout(True)
    plt.draw()
    if save_figs:
        plt.savefig(output)


def plot_tdos(axis_range=None, ISPIN=None, input_file='DOSCAR', display=True,
              on_figs=None, return_refs=False, save_figs=False, save_data=False, output_prefix='TDOS',
              return_states_at_Ef=False):
    """
    Plot the total density of states, with consideration of spin-polarization.
    Accepts input file 'DOSCAR', or 'vasprun.xml'.

    Parameters
    ----------
    axis_range: list
        the range of axes x and y, 4 values in a list
    ISPIN: int
        user specified ISPIN
        If not given, for DOSCAR-type input, infer from OUTCAR/INCAR.
        For vasprun.xml-type input, infer from 'vasprun.xml'.
    input_file: string
        input file name, default to 'DOSCAR'
        For DOSCAR-type, can be any string containing 'DOSCAR'.
        For vasprun.xml-type input, can be any string ending with '.xml'.
    display: bool
        Display figures or not. Default to True.
    on_figs: list/int
        the current figure numbers to plot to, default to new figures
    return_refs: bool
        Return the axes reference(s) drawing or not. Default to False.
    save_figs: bool
        Save figures or not. Default to False. 
    save_data: bool
        Save data or not. Default to False.
    output_prefix: string
        prefix string before the output files, default to 'TDOS'
    return_states_at_Ef: bool
        Calculate the TDOS at Ef with a 0.4 eV window of integration or not. Default to False.

    Returns
    -------
    a dict, containing
        'data': a dict that has 2D array of data,
            easily to Pandas DataFrame by pd.DataFrame(**returned_dict['data'])
        'ax': the axes reference, if return_refs == True
    """
    # get data
    if re.match(r".*\.xml", input_file):
        root = parse(input_file)

        NEDOS = int(root.find("./parameters/separator[@name='dos']/i[@name='NEDOS']").text)
        Ef = float(root.find("./calculation/dos/i[@name='efermi']").text)
        if ISPIN:
            print "Using user specified ISPIN."
        else:
            ISPIN = int(root.find(
                "./parameters/separator[@name='electronic']/separator[@name='electronic spin']/i[@name='ISPIN']").text)

        if ISPIN == 1:
            data = np.zeros((NEDOS, 3))
            for n_step, elem in enumerate(root.findall(
                    "./calculation/dos/total/array/set/set[@comment='spin 1']/r")):
                data[n_step] = elem.text.split()

        elif ISPIN == 2:
            data1 = np.zeros((NEDOS, 3))
            for n_step, elem in enumerate(root.findall(
                    "./calculation/dos/total/array/set/set[@comment='spin 1']/r")):
                data1[n_step] = elem.text.split()
            data2 = np.zeros((NEDOS, 3))
            for n_step, elem in enumerate(root.findall(
                    "./calculation/dos/total/array/set/set[@comment='spin 2']/r")):
                data2[n_step] = elem.text.split()

    elif re.match(r".*DOSCAR.*", input_file):
        with open(input_file, 'r') as f:
            DOSCAR = f.readlines()
        for i in range(len(DOSCAR)):
            DOSCAR[i] = DOSCAR[i].split()

        NEDOS = int(DOSCAR[5][2])
        Ef = float(DOSCAR[5][3])
        if ISPIN:
            print "Using user specified ISPIN."
        else:
            ISPIN = determine_tag_value('ISPIN')

        data = np.array(DOSCAR[6:6 + NEDOS], dtype=float)
        if ISPIN == 2:
            data1 = data[:, [0, 1, 3]]
            data2 = data[:, [0, 2, 4]]

    # start plotting
    figs_assert(on_figs, ISPIN, 'tdos')

    if ISPIN == 1:
        col_names = ['E', 'total', 'integrated']
        data[:, 0] -= Ef
        initiate_figs(on_figs)
        plt.plot(data[:, 0], data[:, 1])
        ax = plt.gca()
        plot_helper_settings(axis_range, 'tdos', save_figs, output=output_prefix + '.pdf')
        axes = {'ax': ax}
        return_dict = {'data': {'columns': col_names, 'data': data}}
        if save_data:
            np.savetxt(output_prefix + '.txt', data, '%15.6E', header=' '.join(col_names))

        if return_states_at_Ef:
            half_window = 0.2
            TDOS_at_Ef = (data[abs(data[:, 0] - 0.2).argmin(), 2] -
                          data[abs(data[:, 0] + 0.2).argmin(), 2]) / half_window / 2.

    elif ISPIN == 2:
        col_names1 = ['E', 'total_up', 'integrated_up']
        col_names2 = ['E', 'total_down', 'integrated_down']
        data1[:, 0] -= Ef
        data2[:, 0] -= Ef
        # Plot the combined TDOS
        initiate_figs(on_figs)
        plt.plot(data1[:, 0], data1[:, 1] + data2[:, 1], label='spin up + down')
        ax1 = plt.gca()
        plot_helper_settings(axis_range, 'tdos', save_figs, output=output_prefix + '-spin-combined.pdf')
        # Plot the overlapping TDOS
        initiate_figs(on_figs)
        plt.plot(data1[:, 0], data1[:, 1], label='spin up')
        plt.plot(data2[:, 0], data2[:, 1], label='spin down')
        ax2 = plt.gca()
        axis_range_copy = None
        if axis_range:
            axis_range_copy = axis_range[:]
            axis_range_copy[3] /= 2.
        plot_helper_settings(axis_range_copy, 'tdos', save_figs, output=output_prefix + '-spin-overlapping.pdf')
        axes = {'ax_spin_combined': ax1, 'ax_spin_overlapping': ax2}
        return_dict = {'data_spin_up': {'columns': col_names1, 'data': data1},
                       'data_spin_down': {'columns': col_names2, 'data': data2}}
        if save_data:
            np.savetxt(output_prefix + '_spin_up.txt', data1, '%15.6E', header=' '.join(col_names1))
            np.savetxt(output_prefix + '_spin_down.txt', data2, '%15.6E', header=' '.join(col_names2))

        if return_states_at_Ef:
            half_window = 0.2
            TDOS_at_Ef = (data1[abs(data1[:, 0] - half_window).argmin(), 2] -
                          data1[abs(data1[:, 0] + half_window).argmin(), 2] +
                          data2[abs(data2[:, 0] - half_window).argmin(), 2] -
                          data2[abs(data2[:, 0] + half_window).argmin(), 2]) / half_window / 2.

    return_dict = display_or_close_figs(display, return_refs, return_dict, axes)
    if return_states_at_Ef:
        return_dict['TDOS_at_Ef'] = TDOS_at_Ef
    return return_dict


def plot_ldos(atom, axis_range=None, ISPIN=None, LORBIT=None, input_file='DOSCAR', display=True,
              on_figs=None, return_refs=False, save_figs=False, save_data=False, output_prefix='LDOS'):
    """
    Plot the local projected density of states, with consideration of spin-polarization.
    Accepts input file 'DOSCAR', or 'vasprun.xml'.

    Parameters
    ----------
    atom: int
        the atom number in DOSCAR/POSCAR interested, counting from 1
    axis_range: list
        the range of axes x and y, 4 values in a list
    ISPIN: int
        user specified ISPIN
        If not given, for DOSCAR-type input, infer from OUTCAR/INCAR.
        For vasprun.xml-type input, infer from 'vasprun.xml'.
    LORBIT: int
        user specified LORBIT
        If not given, for DOSCAR-type input, infer from OUTCAR/INCAR.
        For vasprun.xml-type input, infer from 'vasprun.xml'.
    input_file: string
        input file name, default to 'DOSCAR'
        For DOSCAR-type, can be any string containing 'DOSCAR'.
        For vasprun.xml-type input, can be any string ending with '.xml'.
    display: bool
        Display figures or not. Default to True.
    on_figs: list/int
        the current figure numbers to plot to, default to new figures
    return_refs: bool
        Return the axes reference(s) drawing or not. Default to False.
    save_figs: bool
        Save figures or not. Default to False. 
    save_data: bool
        Save data or not. Default to False.
    output_prefix: string
        prefix string before the output files, default to 'LDOS'

    Returns
    -------
    a dict, containing
        'data': a dict that has 2D array of data,
            easily to Pandas DataFrame by pd.DataFrame(**returned_dict['data'])
        'ax': the axes reference, if return_refs == True
    """
    # get data
    if re.match(r".*\.xml", input_file):
        root = parse(input_file)

        NEDOS = int(root.find("./parameters/separator[@name='dos']/i[@name='NEDOS']").text)
        Ef = float(root.find("./calculation/dos/i[@name='efermi']").text)
        if ISPIN:
            print "Using user specified ISPIN."
        else:
            ISPIN = int(root.find(
                "./parameters/separator[@name='electronic']/separator[@name='electronic spin']/i[@name='ISPIN']").text)
        # vasprun.xml's LORBIT is not correct
        if LORBIT:
            print "Using user specified LORBIT."
        else:
            LORBIT = determine_tag_value('LORBIT')

        if ISPIN == 1:
            if LORBIT == 10 or LORBIT == 0:
                data = np.zeros((NEDOS, 4))
            elif LORBIT == 11 or LORBIT == 1:
                data = np.zeros((NEDOS, 10))
            for n_step, elem in enumerate(root.findall(
                                    "./calculation/dos/partial/array/set/set[@comment='ion " + str(
                                    atom) + "']/set[@comment='spin 1']/r")):
                data[n_step] = elem.text.split()

        elif ISPIN == 2:
            if LORBIT == 10 or LORBIT == 0:
                data1 = np.zeros((NEDOS, 4))
                data2 = np.zeros((NEDOS, 4))
            elif LORBIT == 11 or LORBIT == 1:
                data1 = np.zeros((NEDOS, 10))
                data2 = np.zeros((NEDOS, 10))

            for n_step, elem in enumerate(root.findall(
                                    "./calculation/dos/partial/array/set/set[@comment='ion " + str(
                                    atom) + "']/set[@comment='spin 1']/r")):
                data1[n_step] = elem.text.split()
            for n_step, elem in enumerate(root.findall(
                                    "./calculation/dos/partial/array/set/set[@comment='ion " + str(
                                    atom) + "']/set[@comment='spin 2']/r")):
                data2[n_step] = elem.text.split()

    elif re.match(r".*DOSCAR.*", input_file):
        with open(input_file, 'r') as f:
            DOSCAR = f.readlines()
        for i in range(len(DOSCAR)):
            DOSCAR[i] = DOSCAR[i].split()

        NEDOS = int(DOSCAR[5][2])
        Ef = float(DOSCAR[5][3])
        if ISPIN:
            print "Using user specified ISPIN."
        else:
            ISPIN = determine_tag_value('ISPIN')
        if LORBIT:
            print "Using user specified LORBIT."
        else:
            LORBIT = determine_tag_value('LORBIT')

        data = np.array(DOSCAR[(6 + (NEDOS + 1) * atom):(6 + (NEDOS + 1) * atom + NEDOS)], dtype=float)
        if ISPIN == 2:
            if LORBIT == 10 or LORBIT == 0:
                data1 = data[:, [0, 1, 3, 5]]
                data2 = data[:, [0, 2, 4, 6]]
            elif LORBIT == 11 or LORBIT == 1:
                data1 = data[:, [0, 1, 3, 5, 7, 9, 11, 13, 15, 17]]
                data2 = data[:, [0, 2, 4, 6, 8, 10, 12, 14, 16, 18]]

    # start plotting
    figs_assert(on_figs, ISPIN, 'ldos')

    if ISPIN == 1:
        data[:, 0] -= Ef
        initiate_figs(on_figs)
        if LORBIT == 10 or LORBIT == 0:
            col_names = ['E', 's', 'p', 'd']
            for i in range(1, 4):
                plt.plot(data[:, 0], data[:, i], label=col_names[i])
        elif LORBIT == 11 or LORBIT == 1:
            col_names = ['E', 's', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z2', 'd_xz', 'd_x2y2']
            for i in range(1, 10):
                plt.plot(data[:, 0], data[:, i], label=col_names[i])
        ax = plt.gca()
        plot_helper_settings(axis_range, 'ldos', save_figs, output=output_prefix + '.pdf')
        axes = {'ax': ax}
        return_dict = {'data': {'columns': col_names, 'data': data}}
        if save_data:
            np.savetxt(output_prefix + '.txt', data, '%15.6E', header=' '.join(col_names))

    elif ISPIN == 2:
        data1[:, 0] -= Ef
        data2[:, 0] -= Ef
        # plot spin combined
        initiate_figs(on_figs)
        if LORBIT == 10 or LORBIT == 0:
            col_names1 = ['E', 's_up', 'p_up', 'd_up']
            col_names2 = ['E', 's_down', 'p_down', 'd_down']
            for i in range(1, 4):
                plt.plot(data1[:, 0], data1[:, i] + data2[:, i], label=col_names1[i] + ' + ' + col_names2[i])
        elif LORBIT == 11 or LORBIT == 1:
            col_names1 = ['E', 's_up', 'p_y_up', 'p_z_up', 'p_x_up', 'd_xy_up', 'd_yz_up', 'd_z2_up', 'd_xz_up',
                          'd_x2y2_up']
            col_names2 = ['E', 's_down', 'p_y_down', 'p_z_down', 'p_x_down', 'd_xy_down', 'd_yz_down', 'd_z2_down',
                          'd_xz_down', 'd_x2y2_down']
            for i in range(1, 10):
                plt.plot(data1[:, 0], data1[:, i] + data2[:, i], label=col_names1[i] + ' + ' + col_names2[i])
        ax1 = plt.gca()
        plot_helper_settings(axis_range, 'ldos', save_figs, output=output_prefix + '-spin-combined.pdf')
        # plot spin overlapping
        initiate_figs(on_figs)
        if LORBIT == 10 or LORBIT == 0:
            for i in range(1, 4):
                plt.plot(data1[:, 0], data1[:, i], label=col_names1[i])
                plt.plot(data2[:, 0], -data2[:, i], label=col_names2[i])
        elif LORBIT == 11 or LORBIT == 1:
            for i in range(1, 10):
                plt.plot(data1[:, 0], data1[:, i], label=col_names1[i])
                plt.plot(data2[:, 0], -data2[:, i], label=col_names2[i])
        ax2 = plt.gca()
        axis_range_copy = None
        if axis_range:
            axis_range_copy = axis_range[:]
            axis_range_copy[3] /= 2.
            axis_range_copy[2] /= -axis_range_copy[3]
        plot_helper_settings(axis_range_copy, 'ldos', save_figs, output=output_prefix + '-spin-overlapping.pdf')
        axes = {'ax_spin_combined': ax1, 'ax_spin_overlapping': ax2}
        return_dict = {'data_spin_up': {'columns': col_names1, 'data': data1},
                       'data_spin_down': {'columns': col_names2, 'data': data2}}
        if save_data:
            np.savetxt(output_prefix + '_spin_up.txt', data1, '%15.6E', header=' '.join(col_names1))
            np.savetxt(output_prefix + '_spin_down.txt', data2, '%15.6E', header=' '.join(col_names2))

    return_dict = display_or_close_figs(display, return_refs, return_dict, axes)
    return return_dict


def plot_cohp(bond, axis_range=None, ISPIN=None, input_file='COHPCAR.lobster', display=True,
              on_figs=None, return_refs=False, save_figs=False, save_data=False, output_prefix='COHP'):
    """
    Plot the -COHP, with consideration of spin-polarization.

    Parameters
    ----------
    bond: int
        the bond number in COHPCAR.lobster interested, counting from 1
    axis_range: list
        the range of axes x and y, 4 values in a list
    ISPIN: int
        user specified ISPIN
        If not given, infer from OUTCAR/INCAR.
    input_file: string
        input file name, default to 'COHPCAR.lobster'
    display: bool
        Display figures or not. Default to True.
    on_figs: list/int
        the current figure numbers to plot to, default to new figures
    return_refs: bool
        Return the axes reference(s) drawing or not. Default to False.
    save_figs: bool
        Save figures or not. Default to False. 
    save_data: bool
        Save data or not. Default to False.
    output_prefix: string
        prefix string before the output files, default to 'COHP'

    Returns
    -------
    a dict, containing
        'data': a dict that has 2D array of data,
            easily to Pandas DataFrame by pd.DataFrame(**returned_dict['data'])
        'ax': the axes reference, if return_refs == True
    """
    # get data
    with open(input_file, 'r') as f:
        COHPCAR = f.readlines()

    for line_num, line in enumerate(COHPCAR):
        if re.match(r'No\.\d*:.*\(.*\)', line):
            break
    N_headerlines = line_num

    for line_num, line in enumerate(COHPCAR[N_headerlines:]):
        if not re.match(r'No\.\d*:.*\(.*\)', line):
            break
    N_bonds = line_num

    data_start_line = N_headerlines + N_bonds

    for i in range(len(COHPCAR)):
        COHPCAR[i] = COHPCAR[i].split()

    NEDOS = int(COHPCAR[1][2])
    if ISPIN:
        print "Using user specified ISPIN."
    else:
        ISPIN = determine_tag_value('ISPIN')

    data = np.array(COHPCAR[data_start_line:data_start_line + NEDOS], dtype=float)

    # start plotting
    figs_assert(on_figs, ISPIN, 'cohp')

    if ISPIN == 1:
        col_names = ['E', 'avg', 'avg_integrated']
        for n_bond in range(1, N_bonds + 1):
            col_names.extend(['bond_{0}'.format(n_bond), 'bond_{0}_integrated'.format(n_bond)])

        col_num = bond * 2 + 1
        initiate_figs(on_figs)
        plt.plot(data[:, 0], -data[:, col_num])
        ax = plt.gca()
        plot_helper_settings(axis_range, 'cohp', save_figs, output=output_prefix + '.pdf')
        axes = {'ax': ax}
        return_dict = {'data': {'columns': col_names, 'data': data}}
        if save_data:
            np.savetxt(output_prefix + '.txt', data, '%15.6E', header=' '.join(col_names))

    elif ISPIN == 2:
        col_names1 = ['E', 'avg_up', 'avg_integrated_up']
        col_names2= ['E', 'avg_down', 'avg_integrated_down']
        for n_bond in range(1, N_bonds + 1):
            col_names1.extend(['bond_{0}_up'.format(n_bond), 'bond_{0}_integrated_up'.format(n_bond)])
            col_names2.extend(['bond_{0}_down'.format(n_bond), 'bond_{0}_integrated_down'.format(n_bond)])
        data1 = data[:, :(3+2*N_bonds)]
        data2 = data[:, (3+2*N_bonds):(5+4*N_bonds)]
        data2 = np.column_stack((data[:, 0], data2))

        col_bond = bond * 2 + 1
        # Plot the combined COHP
        initiate_figs(on_figs)
        plt.plot(data1[:, 0], -data1[:, col_bond] - data2[:, col_bond], label='spin up + down')
        ax1 = plt.gca()
        plot_helper_settings(axis_range, 'cohp', save_figs, output=output_prefix + '-spin-combined.pdf')
        # Plot the overlapping COHP
        initiate_figs(on_figs)
        plt.plot(data1[:, 0], -data1[:, col_bond], label='spin up')
        plt.plot(data2[:, 0], -data2[:, col_bond], label='spin down')
        ax2 = plt.gca()
        axis_range_copy = None
        if axis_range:
            axis_range_copy = axis_range[:]
            axis_range_copy[3] /= 2.
        plot_helper_settings(axis_range_copy, 'cohp', save_figs, output=output_prefix + '-spin-overlapping.pdf')
        axes = {'ax_spin_combined': ax1, 'ax_spin_overlapping': ax2}
        return_dict = {'data_spin_up': {'columns': col_names1, 'data': data1},
                       'data_spin_down': {'columns': col_names2, 'data': data2}}
        if save_data:
            np.savetxt(output_prefix + '_spin_up.txt', data, '%15.6E', header=' '.join(col_names1))
            np.savetxt(output_prefix + '_spin_down.txt', data, '%15.6E', header=' '.join(col_names2))

    return_dict = display_or_close_figs(display, return_refs, return_dict, axes)
    return return_dict
