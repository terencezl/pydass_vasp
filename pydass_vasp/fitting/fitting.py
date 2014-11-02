import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# internal use
def get_r_squared(Y, Y_fit_eqlen):
    ss_resid = np.sum((Y_fit_eqlen - Y) ** 2)
    Y_avg = np.sum(Y) / len(Y)
    ss_total = np.sum((Y - Y_avg) ** 2)
    r_squared = 1 - ss_resid / ss_total
    return r_squared


# internal use
def B_M_eqn(V, V0, B0, B0_prime, E0):
    return (E0 + 9 * V0 * B0 / 16. * (((V0 / V) ** (2 / 3.) - 1) ** 3 * B0_prime +
        ((V0 / V) ** (2 / 3.) - 1) ** 2 * (6 - 4 * (V0 / V) ** (2 / 3.))))

# internal use
def B_M_eqn_fixB0prime(B0_prime):
    def B_M_eqn(V, V0, B0, E0):
        return (E0 + 9 * V0 * B0 / 16. * (((V0 / V) ** (2 / 3.) - 1) ** 3 * B0_prime +
            ((V0 / V) ** (2 / 3.) - 1) ** 2 * (6 - 4 * (V0 / V) ** (2 / 3.))))
    return B_M_eqn

# internal use
def B_M_eqn_pv(V, V0, B0, B0_prime):
    return (3*B0/2. * ((V0/V)**(7/3.) - (V0/V)**(5/3.)) *
        (1 + 3/4. * (B0_prime - 4) * ((V0/V)**(2/3.) - 1)))

# internal use
def B_M_eqn_pv_fixB0prime(B0_prime):
    def B_M_eqn_pv(V, V0, B0):
        return (3*B0/2. * ((V0/V)**(7/3.) - (V0/V)**(5/3.)) *
            (1 + 3/4. * (B0_prime - 4) * ((V0/V)**(2/3.) - 1)))
    return B_M_eqn_pv


def eos_fit(V, Y, p_v=False, fix_B0_prime=False, plot=False, new_fig=False, save_fig=False):
    """
    Fit the volume and total energy, or pressure to the Birch-Murnaghan equation of state.

    Parameters
    ----------
    V: arrays of volume (Angstrom^3)
    Y: arrays of total energy (eV) or pressure (GPa)

    Optional

    fix_B0_prime: Default to False, or give it a value, like 4
    plot: if plot the fitted arrays
    new_fig: if plot on new figures, rather than existing ones
    save_fig: if save the figure to a file. When passed a string, it will be the filename.

    Returns
    -------
    a dict, containing parameters, r-squared value and fitted_arrays
    """
    V = np.array(V)
    Y = np.array(Y)
    if not p_v:
        if not fix_B0_prime:
            initial_parameters = [V.mean(), 2.5, 4, Y.mean()]
            fit_eqn = B_M_eqn
        else:
            initial_parameters = [V.mean(), 2.5, Y.mean()]
            fit_eqn = B_M_eqn_fixB0prime(fix_B0_prime)
    else:
        if not fix_B0_prime:
            initial_parameters = [V.mean(), 2.5, 4]
            fit_eqn = B_M_eqn_pv
        else:
            initial_parameters = [V.mean(), 2.5]
            fit_eqn = B_M_eqn_pv_fixB0prime(fix_B0_prime)

    coeffs, pcov = curve_fit(fit_eqn, V, Y, initial_parameters)
    V_fit = np.linspace(sorted(V)[0], sorted(V)[-1], 1000)
    Y_fit = fit_eqn(V_fit, *coeffs)
    Y_fit_eqlen = fit_eqn(V, *coeffs)
    r_squared = get_r_squared(Y, Y_fit_eqlen)

    parameters = {'V0': coeffs[0], 'B0': coeffs[1]}
    if not fix_B0_prime:
        parameters['B0_prime'] = coeffs[2]
        if not p_v:
            parameters['E0'] = coeffs[3]
    else:
        parameters['E0'] = coeffs[2]

    if plot:
        if new_fig:
            plt.figure()
        plt.plot(V, Y, 'o')
        plt.plot(V_fit, Y_fit, '-')
        plt.xlabel(r'V ($\AA^{3}$)')
        if not p_v:
            plt.ylabel('E (eV)')
        else:
            plt.ylabel('P (GPa)')
        plt.tight_layout()
        if save_figs:
            if isinstance(save_figs, bool):
                plt.savefig('eos_fit.pdf')
            else:
                plt.savefig(save_figs)

    columns = ['V_fit', 'P_fit']
    if not p_v:
        columns[1] = 'E_fit'

    return {'parameters': parameters, 'r_squared': r_squared,
            'fitted_arrays': {'columns': columns, 'data': np.column_stack((V_fit, Y_fit))}}


def polyfit(X, Y, order, plot=False, new_fig=False, save_fig=False):
    """
    Fit the two equal length arrays to a polynomial with a specified order.

    Parameters
    ----------
    X, Y: arrays
    order: the polynomial order

    Optional

    plot: if plot the fitted arrays
    new_fig: if plot on new figures, rather than existing ones
    save_fig: if save the figure to a file. When passed a string, it will be the filename.

    Returns
    -------
    a dict, containing parameters, r-squared value and fitted_arrays
    """
    popts = np.polyfit(X, Y, order)
    p = np.poly1d(popts)
    X_fit = np.linspace(sorted(X)[0], sorted(X)[-1], 1000)
    Y_fit = p(X_fit)
    Y_fit_eqlen = p(X)
    r_squared = get_r_squared(Y, Y_fit_eqlen)

    if plot:
        if new_fig:
            plt.figure()
        plt.plot(X, Y, 'o')
        plt.plot(X_fit, Y_fit, '-')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.tight_layout()
        if save_figs:
            if isinstance(save_figs, bool):
                plt.savefig('polyfit.pdf')
            else:
                plt.savefig(save_figs)

    return {'coeffs': p, 'r_squared': r_squared,
        'fitted_arrays': {'columns': ['X_fit', 'Y_fit'], 'data': np.column_stack((X_fit, Y_fit))}}