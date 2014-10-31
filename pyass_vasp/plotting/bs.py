import re
import numpy as np
import matplotlib.pyplot as plt
from helpers import determine_tag_value, plot_helper_figs_assert, plot_helper_figs, plot_helper_post

###### integrate plotting function #####


# Effective mass calculation funcitons.
def find_band_edges(kp_edge, prec_range, E):
    """
    Given the k-point index number on the x-axis, search for the indices of bands
    right below or above Ef.
    Best used in an interactive Python interpreter and having pyplot.ion().

    Parameters
    ----------
    kp_edge : int
        The k-point index number that corresponds to the band edges (VBM and CBM).
    prec_range : float
        The searching range in eV.
    E : 2D numpy array
        The 2D array that contains eigenvalues. Each row denotes a band, and each column a k-point.
    """
    # Examine valence band edge.
    print 'The possible valence bands are', \
        np.where(np.logical_and(E[:, kp_edge] > -prec_range, E[:, kp_edge] < 0))[0]
    # Examine conduction band edge.
    print 'The possible conduction bands are', \
        np.where(np.logical_and(E[:, kp_edge] < prec_range, E[:, kp_edge] > 0))[0]


def get_reduced_effective_mass(band, kp_start, kp_end, kp_linearized_array, E):
    """
    Given the band index number, k-point start and end indices, fit the included curve
    to a 2nd-order polynomial, and obtain the effective mass of the carrier electron or hole.
    Best used in an interactive Python interpreter and having pyplot.ion().

    Parameters
    ----------
    band : int
        The band index number of interest.
    kp_start : int
        The index number of the starting k-point.
    kp_end : int
        The index number of the ending k-point.
    kp_linearized_array : 1D numpy array
        The full x-axis of k-points.
    E : 2D numpy array
        The 2D array that contains eigenvalues. Each row denotes a band, and each column a k-point.

    Returns
    -------
    effective_mass_reduced : float
        The reduced effective mass.
    """
    h_bar = 1.054571726e-34
    e = 1.6021176462e-19
    m_e = 9.10938291e-31
    scaling_const = 6.3743775177e-10

    # Decide on the fitting range, characterized by indices.
    selected_kp_array = kp_linearized_array[kp_start:kp_end + 1]
    selected_energy_array = E[band, kp_start:kp_end + 1]
    p = np.poly1d(np.polyfit(selected_kp_array, selected_energy_array, 2))
    axis_fitted = -p[1]/2/p[2]
    axis_actual = selected_kp_array[selected_energy_array.argmin() if p[2] > 0 else selected_energy_array.argmax()]
    print "The fitted x coord at energy extrema is {0}, and the actual is {1}.".format(axis_fitted, axis_actual)
    k_fit = np.linspace(kp_linearized_array[kp_start], kp_linearized_array[kp_end], 200)
    plt.plot(k_fit, p(k_fit), lw=2)

    d2E_dk2 = e * p[2] / (2 * np.pi / scaling_const) ** 2
    effective_mass_reduced = h_bar ** 2 / d2E_dk2 / m_e
    return effective_mass_reduced


def plot_bs(axis_range=None, ISPIN=None, Ef=None, EIGENVAL='EIGENVAL', display=True,
    on_figs=None, save_figs=False, save_data=False, output_prefix='BS'):
    """
    Plot the band structure, with consideration of spin-polarization.

    Parameters
    ----------
    axis_range: list
        the range of axes x and y, 4 values in a list
    ISPIN: int
        user specified ISPIN. If not given, infer from OUTCAR, or INCAR
    Ef: float
        user specified Ef. If not given, infer from DOSCAR, or OUTCAR
    EIGENVAL: string
        EIGENVAL file name, default to 'EIGENVAL'
    display: bool
        display figures or not
    on_figs: list
        the current figure numbers to plot to, default to new figures
    save_figs: bool
        save figures or not
    save_data: bool
        save data or not
    output_prefix: string
        prefix string before the output files, default to 'BS'

    Returns
    -------
    a dict, containing 'reciprocal_points', the symbols of those points,
    'reciprocal_point_locations_on_axis', their locations on the x-axis,
    'E_matrix', energy eigenvalue matrix, with the shape [N_bands, N_kps],
    'bs_data', which is another dict, containing
        'data', a 2D numpy array of data from EIGENVAL, and
        'columns', a list of column labels
    """

    # from vasprun.xml
    from ..xml_utils import parse
    root = parse('vasprun.xml')

    if Ef:
        print "Using user specified Ef."
    else:
        Ef = float(root.find("./calculation/dos/i[@name='efermi']").text)

    if ISPIN:
        print "Using user specified ISPIN."
    else:
        ISPIN = int(root.find(
        "./parameters/separator[@name='electronic']/separator[@name='electronic spin']/i[@name='ISPIN']").text)

    plot_helper_figs_assert(on_figs, ISPIN, 'bs')

    # with open(EIGENVAL, 'r') as f:
    #     EIGENVAL = f.readlines()
    # for i in range(len(EIGENVAL)):
    #     EIGENVAL[i] = EIGENVAL[i].split()

    # # How many KPs in total? Can be found in EIGENVAL, 6th line, 2nd number.
    # N_kps = int(EIGENVAL[5][1])
    # # How many bands are to be drawn? 6th line, 3rd number.
    # N_bands = int(EIGENVAL[5][2])

    # with open('KPOINTS', 'r') as f:
    #     KPOINTS = f.readlines()
    # N_kps_per_section = int(KPOINTS[1])
    # N_sections = N_kps / N_kps_per_section

    N_kps_per_section = int(root.find("./kpoints/generation[@param='listgenerated']/i[@name='divisions']").text)
    N_sections = len(root.findall("./kpoints/generation[@param='listgenerated']/v")) - 1
    N_kps = len(root.findall("./kpoints/varray[@name='kpointlist']/v"))
    assert N_kps_per_section * N_sections == N_kps, \
        "The product of N of kpoints per section and N of sections does not match N of total kpoints. Strange."
    N_bands = int(root.find("./parameters/separator[@name='electronic']/i[@name='NBANDS']").text)


    ##### first try vasprun.xml. make fallback to OUTCAR #####

    # Get the start and end point coordinate of each section. From OUTCAR.
    kp_list = np.zeros((N_kps, 3))
    with open('OUTCAR', 'r') as f:
        for line in f:
            if re.match(r".*k-points in units of 2pi/SCALE and weight:.*", line):
                kp_end_symbol_list = line.replace(
                    'k-points in units of 2pi/SCALE and weight:', '').strip().split('-')
                break
        for kp in range(N_kps):
            kp_list[kp] = f.next().split()[:3]

    kp_section_start_end_pair_array = np.zeros((N_sections, 2, 3))
    for section in range(N_sections):
        kp_section_start_end_pair_array[section] = [kp_list[N_kps_per_section * section],
                                                    kp_list[N_kps_per_section * (section + 1) - 1]]

    # Gerenate the linearized kp_linearized_array as x-axis.
    kp_end_point = 0
    kp_end_point_list = [0] * (N_sections + 1)
    kp_section_linearized_array = np.zeros((N_sections, N_kps_per_section))
    for section, section_coord_pair in enumerate(kp_section_start_end_pair_array):
        kp_end_point_next = kp_end_point + np.linalg.norm(
            section_coord_pair[1] - section_coord_pair[0])
        kp_end_point_list[section + 1] = kp_end_point_next
        kp_section_linearized_array[section] = np.linspace(kp_end_point,
                                                           kp_end_point_next, N_kps_per_section)
        kp_end_point = kp_end_point_next

    kp_linearized_array = kp_section_linearized_array.flatten()

    # # Get energy for each band, each kpoint step.
    # E = np.zeros((N_bands, N_kps))
    # for n_b in range(0, N_bands):
    #     for n_s in range(0, N_kps):
    #         E[n_b, n_s] = float(EIGENVAL[8 + n_b + (N_bands + 2) * n_s][1])
    # E -= Ef

    # for the returned data column names
    col_names = [str(i) for i in range(N_bands + 1)]
    col_names[0] = 'k_points'

    if ISPIN == 2:
        E_spin1 = np.zeros((N_kps, N_bands))
        for n_kp in range(N_kps):
            for n_band, elem in enumerate(root.findall(
                "./calculation/eigenvalues/array/set/set[@comment='spin 1']/set[@comment='kpoint "
                + str(n_kp) + "']/r")):
                E_spin1[n_kp, n_band] = elem.text.split()[0]
        E_spin1 -= Ef

        E_spin2 = np.zeros((N_kps, N_bands))
        for n_kp in range(N_kps):
            for n_band, elem in enumerate(root.findall(
                "./calculation/eigenvalues/array/set/set[@comment='spin 2']/set[@comment='kpoint "
                + str(n_kp) + "']/r")):
                E_spin2[n_kp, n_band] = elem.text.split()[0]
        E_spin2 -= Ef

        # Plot the bands of spin up.
        plot_helper_figs(on_figs)
        ax1 = plt.subplot(111)
        for band in range(N_bands):
            plt.plot(kp_linearized_array, E_spin1[:, band])

        plt.xlim(kp_end_point_list[0], kp_end_point_list[-1])
        ax1.xaxis.set_ticks(kp_end_point_list)
        ax1.xaxis.set_ticklabels(kp_end_symbol_list)
        plt.axhline(0, ls='--', c='k', alpha=0.5)
        if axis_range:
            plt.ylim(axis_range[0], axis_range[1])
        for kp_end_point in range(len(kp_end_point_list)):
            plt.axvline(kp_end_point_list[kp_end_point], ls='--', c='k', alpha=0.5)
        plt.ylabel('Energy (eV)')
        try:
            plt.tight_layout()
        except RuntimeError:
            print "Tight layout failed... Not a big deal though."
        plt.savefig(output_prefix + '-spin-up.pdf')
        plot_helper_post()

        # Plot the bands of spin down.
        plot_helper_figs(on_figs)
        ax2 = plt.subplot(111)
        for band in range(N_bands):
            plt.plot(kp_linearized_array, E_spin2[:, band])

        plt.xlim(kp_end_point_list[0], kp_end_point_list[-1])
        ax2.xaxis.set_ticks(kp_end_point_list)
        ax2.xaxis.set_ticklabels(kp_end_symbol_list)
        plt.axhline(0, ls='--', c='k', alpha=0.5)
        if axis_range:
            plt.ylim(axis_range[0], axis_range[1])
        for kp_end_point in range(len(kp_end_point_list)):
            plt.axvline(kp_end_point_list[kp_end_point], ls='--', c='k', alpha=0.5)
        plt.ylabel('Energy (eV)')
        try:
            plt.tight_layout()
        except RuntimeError:
            print "Tight layout failed... Not a big deal though."
        plt.savefig(output_prefix + '-spin-down.pdf')
        plot_helper_post()

        # Plot the bands of up and down overlapped.
        plot_helper_figs(on_figs)
        ax3 = plt.subplot(111)
        for band in range(N_bands):
            plt.plot(kp_linearized_array, E_spin1[:, band], 'k')
            plt.plot(kp_linearized_array, E_spin2[:, band], 'b')

        plt.xlim(kp_end_point_list[0], kp_end_point_list[-1])
        ax3.xaxis.set_ticks(kp_end_point_list)
        ax3.xaxis.set_ticklabels(kp_end_symbol_list)
        plt.axhline(0, ls='--', c='k', alpha=0.5)
        if axis_range:
            plt.ylim(axis_range[0], axis_range[1])
        for kp_end_point in range(len(kp_end_point_list)):
            plt.axvline(kp_end_point_list[kp_end_point], ls='--', c='k', alpha=0.5)
        plt.ylabel('Energy (eV)')
        try:
            plt.tight_layout()
        except RuntimeError:
            print "Tight layout failed... Not a big deal though."
        plt.savefig(output_prefix + '-spin-overlapped.pdf')
        plot_helper_post()

        return {
            'reciprocal_points': kp_end_symbol_list,
            'reciprocal_point_locations_on_axis': kp_end_point_list,
            'bs_spin_up_data': {
                'columns': col_names,
                'data': np.columns_stack((kp_linearized_array, E_spin1))
            },
            'bs_spin_down_data': {
                'columns': col_names,
                'data': np.columns_stack((kp_linearized_array, E_spin2))
            },
            'axes': {'ax_spin_up': ax1, 'ax_spin_down': ax2, 'ax_spin_separated': ax3},
        }

    if ISPIN == 1:
        E = np.zeros((N_kps, N_bands))
        for n_kp in range(N_kps):
            for n_band, elem in enumerate(root.findall(
                "./calculation/eigenvalues/array/set/set[@comment='spin 1']/set[@comment='kpoint "
                + str(n_kp) + "']/r")):
                E[n_kp, n_band] = elem.text.split()[0]
        E -= Ef

        # Plot the bands.
        plot_helper_figs(on_figs)
        ax = plt.subplot(111)
        for band in range(N_bands):
            plt.plot(kp_linearized_array, E[:, band])

        plt.xlim(kp_end_point_list[0], kp_end_point_list[-1])
        ax.xaxis.set_ticks(kp_end_point_list)
        ax.xaxis.set_ticklabels(kp_end_symbol_list)
        plt.axhline(0, ls='--', c='k', alpha=0.5)
        if axis_range:
            plt.ylim(axis_range[0], axis_range[1])
        for kp_end_point in range(len(kp_end_point_list)):
            plt.axvline(kp_end_point_list[kp_end_point], ls='--', c='k', alpha=0.5)
        plt.ylabel('Energy (eV)')
        try:
            plt.tight_layout()
        except RuntimeError:
            print "Tight layout failed... Not a big deal though."
        plt.savefig(output_prefix + '.pdf')
        plot_helper_post()

        return {
            'reciprocal_points': kp_end_symbol_list,
            'reciprocal_point_locations_on_axis': kp_end_point_list,
            'bs_data': {
                'columns': col_names,
                'data': np.columns_stack((kp_linearized_array, E_spin1))
            },
            'axes': ax1
        }
