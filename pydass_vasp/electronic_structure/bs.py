import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .helpers import determine_tag_value, figs_assert, initiate_figs
from ..xml_utils import parse


# Effective mass calculation funcitons.
def find_band_edges(kp_edge, prec_range, eigenvalues):
    """
    Given the k-point index number on the x-axis, search for the indices of bands
    right below or above Ef.

    Parameters
    ----------
    kp_edge : int
        The k-point index number that corresponds to the band edges (VBM and CBM).
    prec_range : float
        The searching range in eV.
    eigenvalues : 2D numpy array
        The 2D array that contains eigenvalues. Each row denotes a k-point, and each column a band.
    """
    # Examine valence band edge.
    print('The possible valence bands are', \
        np.where(np.logical_and(eigenvalues[kp_edge] > -prec_range, eigenvalues[kp_edge] < 0))[1])
    # Examine conduction band edge.
    print('The possible conduction bands are', \
        np.where(np.logical_and(eigenvalues[kp_edge] < prec_range, eigenvalues[kp_edge] > 0))[1])


def get_effective_mass(band, kp_start, kp_end, kps_linearized, eigenvalues):
    """
    Given the band index number, k-point start and end indices, fit the included curve
    to a 2nd-order polynomial, and obtain the effective mass of the carrier electron or hole.

    Parameters
    ----------
    band : int
        The band index number of interest.
    kp_start : int
        The index number of the starting k-point.
    kp_end : int
        The index number of the ending k-point.
    kps_linearized : 1D numpy array
        The full x-axis of k-points.
    eigenvalues : 2D numpy array
        The 2D array that contains eigenvalues. Each row denotes a k-point, and each column a band.

    Returns
    -------
    the effective mass
    """
    h_bar = 1.054571726e-34
    e = 1.6021176462e-19
    m_e = 9.10938291e-31
    scaling_const = 6.3743775177e-10

    # Decide on the fitting range, characterized by indices.
    selected_kp_array = kps_linearized[kp_start:kp_end + 1]
    selected_energy_array = eigenvalues[kp_start:kp_end + 1, band]
    p = np.poly1d(np.polyfit(selected_kp_array, selected_energy_array, 2))
    axis_fitted = -p[1]/2/p[2]
    axis_actual = selected_kp_array[selected_energy_array.argmin() if p[2] > 0 else selected_energy_array.argmax()]
    print("The fitted x coord at energy extrema is {0}, and the actual is {1}.".format(axis_fitted, axis_actual))
    k_fit = np.linspace(kps_linearized[kp_start], kps_linearized[kp_end], 200)
    plt.plot(k_fit, p(k_fit), lw=2)

    d2E_dk2 = e * p[2] / (2 * np.pi / scaling_const) ** 2
    effective_mass = h_bar ** 2 / d2E_dk2 / m_e
    return effective_mass


# internal
def plot_helper_settings(ylim, reciprocal_point_locations, reciprocal_point_labels):
    plt.xlim(reciprocal_point_locations[0], reciprocal_point_locations[-1])
    ax = plt.gca()
    ax.xaxis.set_ticks(reciprocal_point_locations)
    ax.xaxis.set_ticklabels(reciprocal_point_labels)
    plt.axhline(0, ls='--', c='k', alpha=0.5)
    if ylim:
        plt.ylim(ylim)
    for kp_end_point in range(len(reciprocal_point_locations)):
        plt.axvline(reciprocal_point_locations[kp_end_point], ls='--', c='k', alpha=0.5)
    plt.ylabel('Energy (eV)')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plt.legend(loc=0, fontsize='small')
    fig = plt.gcf()
    fig.set_tight_layout(True)
    plt.draw()


def get_bs(filepath='EIGENVAL', ISPIN=None, N_kps_per_section=None, reciprocal_point_labels=None, Ef=None, plot=False,
            ylim=None, on_figs=None):
    """
    Plot the band structure, with consideration of spin-polarization.
    Accepts file type 'EIGENVAL', or 'vasprun.xml'.

    Parameters
    ----------
    filepath: string
        file path, default to 'EIGENVAL'
        For EIGENVAL-type file, can be any string containing 'EIGENVAL'.
        For vasprun.xml-type file, can be any string ending with '.xml'.
    ISPIN: int
        user specified ISPIN
        If not given, for EIGENVAL-type file, infer from 'OUTCAR'/'INCAR'.
    N_kps_per_section: int
        user specified number of k-points per line section
        If not given, for EIGENVAL-type file, infer from 'KPOINTS'.
    reciprocal_point_labels: list
        user specified list of reciprocal point string symbols, like ['G','X','A']
        Its length has to be the number of line sections + 1
        If not given, infer from 'OUTCAR'/'KPOINTS'.
    Ef: float
        user specified Ef
        If not given, for EIGENVAL-type file, infer from 'OUTCAR'/'DOSCAR'
    plot: bool
        whether to plot the data, default to False
    ylim: list
        the range of y-axis, 2 values in a list
    on_figs: list/int
        the current figure numbers to plot to, default to new figures

    Returns
    -------
    a dict, containing
        'reciprocal_point_labels', the symbols of those points
        'reciprocal_point_locations', their locations on the x-axis
        'data', a pandas dataframe
        'ax': the axes reference, if return_refs == True
    """
    dirname = os.path.dirname(filepath)
    if re.match(r".*\.xml", filepath):
        root = parse(filepath)

        if ISPIN:
            print("Using user specified ISPIN.")
        else:
            ISPIN = int(root.find(
                "./parameters/separator[@name='electronic']/separator[@name='electronic spin']/i[@name='ISPIN']").text)
        if Ef:
            print("Using user specified Ef.")
        else:
            Ef = float(root.find("./calculation/dos/i[@name='efermi']").text)

        N_bands = int(root.find("./parameters/separator[@name='electronic']/i[@name='NBANDS']").text)
        N_kps_per_section = int(root.find("./kpoints/generation[@param='listgenerated']/i[@name='divisions']").text)
        N_kps = len(root.findall("./kpoints/varray[@name='kpointlist']/v"))
        N_kp_sections = int(N_kps / N_kps_per_section)
        N_kp_sections_alt = len(root.findall("./kpoints/generation[@param='listgenerated']/v")) - 1
        assert N_kp_sections == N_kp_sections_alt, \
            "The product of N of kpoints per section and N of sections does not match N of total kpoints. Strange."

        # get the full k-points list
        kps = np.zeros((N_kps, 3))
        for kp, elem in enumerate(root.findall("./kpoints/varray[@name='kpointlist']/v")):
            kps[kp] = elem.text.split()[:3]

        # get eigenvalues data
        if ISPIN == 1:
            data = np.zeros((N_kps, N_bands))
            for kp in range(N_kps):
                for band, elem in enumerate(root.findall(
                                        "./calculation/eigenvalues/array/set/set[@comment='spin 1']/set[@comment='kpoint "
                                        + str(kp + 1) + "']/r")):
                    data[kp, band] = elem.text.split()[0]

        if ISPIN == 2:
            data1 = np.zeros((N_kps, N_bands))
            for kp in range(N_kps):
                for band, elem in enumerate(root.findall(
                                        "./calculation/eigenvalues/array/set/set[@comment='spin 1']/set[@comment='kpoint "
                                        + str(kp + 1) + "']/r")):
                    data1[kp, band] = elem.text.split()[0]
            data2 = np.zeros((N_kps, N_bands))
            for kp in range(N_kps):
                for band, elem in enumerate(root.findall(
                                        "./calculation/eigenvalues/array/set/set[@comment='spin 2']/set[@comment='kpoint "
                                        + str(kp + 1) + "']/r")):
                    data2[kp, band] = elem.text.split()[0]

    elif re.match(r".*EIGENVAL.*", filepath):
        # get ISPIN
        if ISPIN:
            print("Using user specified ISPIN.")
        else:
            ISPIN = determine_tag_value('ISPIN', filepath)
        # get Ef
        if Ef:
            print("Using user specified Ef.")
        else:
            try:
                with open(os.path.join(dirname, 'OUTCAR'), 'r') as f:
                    for line in f:
                        if re.match(r"\s*E-fermi :.*", line):
                            Ef = float(line.split()[2])
            except IOError:
                try:
                    with open(os.path.join(dirname, 'DOSCAR'), 'r') as f:
                        for i in range(6):
                            line = f.readline()
                    # Fermi energy. Found in DOSCAR, 6th line, 4th number.
                    Ef = float(line.split()[3])
                except IOError:
                    raise IOError("Can't determine Ef! Either manually specify it, or provide OUTCAR/DOSCAR.")

        # read the main file
        with open(filepath, 'r') as f:
            EIGENVAL = f.readlines()
        for i in range(len(EIGENVAL)):
            EIGENVAL[i] = EIGENVAL[i].split()
        # How many bands are to be drawn? 6th line, 3rd number.
        N_bands = int(EIGENVAL[5][2])
        # How many KPs in total? Can be found in EIGENVAL, 6th line, 2nd number.
        N_kps = int(EIGENVAL[5][1])

        # get nkp per sections
        if N_kps_per_section:
            print("Using user specified number of k-points per line section.")
        else:
            try:
                with open(os.path.join(dirname, 'KPOINTS'), 'r') as f:
                    KPOINTS = f.readlines()
                N_kps_per_section = int(KPOINTS[1])
            except IOError:
                raise IOError(
                    "Can't determine number of k-points per line section! Either manually specify it, or provide KPOINTS.")

        N_kp_sections = int(N_kps / N_kps_per_section)

        # get the full k-points list
        kps = np.zeros((N_kps, 3))
        for kp in range(N_kps):
            kps[kp] = [float(i) for i in EIGENVAL[7 + (N_bands + 2) * kp][:3]]

        # get eigenvalues data
        if ISPIN == 1:
            data = np.zeros((N_kps, N_bands))
            for kp in range(N_kps):
                EIGENVAL_at_kp = EIGENVAL[8 + (N_bands + 2) * kp: 8 + N_bands + (N_bands + 2) * kp]
                data[kp] = [float(i[1]) for i in EIGENVAL_at_kp]

        if ISPIN == 2:
            data1 = np.zeros((N_kps, N_bands))
            data2 = np.zeros((N_kps, N_bands))
            for kp in range(N_kps):
                EIGENVAL_at_kp = EIGENVAL[8 + (N_bands + 2) * kp: 8 + N_bands + (N_bands + 2) * kp]
                data1[kp] = [float(i[1]) for i in EIGENVAL_at_kp]
                data2[kp] = [float(i[2]) for i in EIGENVAL_at_kp]

    # confluence of two processing approaches
    # generate the kp_section_pairs
    kp_section_pairs = np.zeros((N_kp_sections, 2, 3))
    for section in range(N_kp_sections):
        kp_section_pairs[section] = [kps[N_kps_per_section * section],
                                     kps[N_kps_per_section * (section + 1) - 1]]

    # generate the linearized kps_linearized as x-axis.
    reciprocal_point_locations = np.zeros(N_kp_sections + 1)
    kps_linearized_sectioned = np.zeros((N_kp_sections, N_kps_per_section))
    for section, section_pair in enumerate(kp_section_pairs):
        reciprocal_point_locations[section + 1] = reciprocal_point_locations[section] + np.linalg.norm(
            section_pair[1] - section_pair[0])
        kps_linearized_sectioned[section] = np.linspace(reciprocal_point_locations[section],
                                                        reciprocal_point_locations[section + 1], N_kps_per_section)
    kps_linearized = kps_linearized_sectioned.flatten()

    # get reciprocal point symbols
    if reciprocal_point_labels:
        print("Using user specified reciprocal point symbols.")
    else:
        try:
            with open(os.path.join(dirname, 'OUTCAR'), 'r') as f:
                for line in f:
                    if re.match(r".*k-points in units of 2pi/SCALE and weight:.*", line):
                        reciprocal_point_labels = line.replace(
                            'k-points in units of 2pi/SCALE and weight:', '').strip().split('-')
                        break
        except IOError:
            try:
                with open(os.path.join(dirname, 'KPOINTS'), 'r') as f:
                    KPOINTS = f.readlines()
                reciprocal_point_labels = KPOINTS[0].strip().split('-')
            except IOError:
                print(
                    "Can't determine reciprocal point symbols! Either manually specify it, or provide OUTCAR/KPOINTS.")

    # make Greek letters real
    reciprocal_point_labels = ['$' + i + '$' if '\\' in i else i for i in reciprocal_point_labels]

    return_dict = {
        'reciprocal_point_locations': reciprocal_point_locations,
        'reciprocal_point_labels': reciprocal_point_labels,
    }

    if ISPIN == 1:
        col_names = ['band_' + str(i) for i in range(N_bands + 1)]
        col_names[0] = 'k_points'
        data -= Ef
        data = np.column_stack((kps_linearized, data))
        return_dict.update({'data': {'columns': col_names, 'data': data}})
    elif ISPIN == 2:
        col_names1 = ['band_' + str(i) + '_up' for i in range(N_bands + 1)]
        col_names1[0] = 'k_points'
        col_names2 = ['band_' + str(i) + '_down' for i in range(N_bands + 1)]
        col_names2[0] = 'k_points'
        data1 -= Ef
        data2 -= Ef
        data1 = np.column_stack((kps_linearized, data1))
        data2 = np.column_stack((kps_linearized, data2))
        return_dict.update({'data_spin_up': pd.DataFrame(**{'columns': col_names1, 'data': data1}),
                            'data_spin_down': pd.DataFrame(**{'columns': col_names2, 'data': data2}),
                            })

    if plot:
        # plot
        # figs_assert(on_figs, ISPIN, 'bs')
        color_cycle = plt.rcParams['axes.color_cycle']
        if ISPIN == 1:
            initiate_figs(on_figs)
            ax = plt.gca()
            for band in range(1, N_bands + 1):
                plt.plot(kps_linearized, data[:, band], color_cycle[0])
            plot_helper_settings(ylim, reciprocal_point_locations, reciprocal_point_labels)
        elif ISPIN == 2:
            initiate_figs(on_figs)
            ax = plt.gca()
            for band in range(1, N_bands + 1):
                plt.plot(kps_linearized, data1[:, band], color_cycle[0], label='spin up')
                plt.plot(kps_linearized, data2[:, band], color_cycle[1], label='spin down')
            plot_helper_settings(ylim, reciprocal_point_locations, reciprocal_point_labels)

        return_dict.update({'ax': ax})

    return return_dict
