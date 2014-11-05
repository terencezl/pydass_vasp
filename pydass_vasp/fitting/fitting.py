import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from ..plotting.helpers import initiate_figs, display_or_close_figs


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


def eos_fit(V, Y, p_v=False, fix_B0_prime=None, display=False, on_figs=None, return_refs=False,
            save_figs=False, save_data=False, output_prefix='eos_fit'):
    """
    Fit the volume and total energy, or pressure to the Birch-Murnaghan equation of state.

    Parameters
    ----------
    V: arrays
        volume (Angstrom^3)
    Y: arrays
        total energy (eV), or pressure (GPa) if p_v == True
    fix_B0_prime: float
        Keep B0_prime fixed to a given value or not. Default to None.
    display: bool
        Display figures or not. Default to False.
    on_figs: list/int
        the current figure numbers to plot to, default to new figures
    return_refs: bool
        Return the axes reference(s) drawing or not. Default to False.
    save_figs: bool
        Save figures or not. Default to False.
    save_data: bool
        Save data or not. Default to False.
    output_prefix: string
        prefix string before the output files, default to 'eos_fit'

    Returns
    -------
    a dict, containing
        'parameters': fitting parameters
        'r-squared': value for evaluating error
        'fitted_data': a dict that has 2D array of fitted data
            easily to Pandas DataFrame by pd.DataFrame(**returned_dict['fitted_data'])
        'ax': the axes reference, if return_refs == True
    """
    V = np.array(V)
    Y = np.array(Y)

    if not fix_B0_prime:
        if not p_v:
            initial_parameters = [V.mean(), 2.5, 4, Y.mean()]
            fit_eqn = B_M_eqn
        else:
            initial_parameters = [V.mean(), 2.5, 4]
            fit_eqn = B_M_eqn_pv
    else:
        if not p_v:
            initial_parameters = [V.mean(), 2.5, Y.mean()]
            fit_eqn = B_M_eqn_fixB0prime(fix_B0_prime)
        else:
            initial_parameters = [V.mean(), 2.5]
            fit_eqn = B_M_eqn_pv_fixB0prime(fix_B0_prime)

    coeffs, pcov = curve_fit(fit_eqn, V, Y, initial_parameters)
    V_fit = np.linspace(sorted(V)[0], sorted(V)[-1], 1000)
    Y_fit = fit_eqn(V_fit, *coeffs)
    Y_fit_eqlen = fit_eqn(V, *coeffs)
    r_squared = get_r_squared(Y, Y_fit_eqlen)
    data = np.column_stack((V_fit, Y_fit))
    if not p_v:
        col_names = ['V_fit', 'E_fit']
    else:
        col_names = ['V_fit', 'P_fit']

    parameters = {'V0': coeffs[0], 'B0': coeffs[1]}
    if not fix_B0_prime:
        parameters['B0_prime'] = coeffs[2]
        if not p_v:
            parameters['E0'] = coeffs[3]
    else:
        if not p_v:
            parameters['E0'] = coeffs[2]

    return_dict = {'parameters': parameters, 'r_squared': r_squared,
                'fitted_data': {'columns': col_names, 'data': data}}

    if display or save_figs or return_refs:
        initiate_figs(on_figs)
        plt.plot(V, Y, 'o')
        plt.plot(V_fit, Y_fit, '-')
        axes = {'ax': plt.gca()}
        plt.xlabel(r'V ($\AA^{3}$)')
        if not p_v:
            plt.ylabel('E (eV)')
        else:
            plt.ylabel('P (GPa)')
        plt.tight_layout()
        if save_figs:
            plt.savefig(output_prefix + '.pdf')
        return_dict = display_or_close_figs(display, return_refs, return_dict, axes)

    if save_data:
        np.savetxt(output_prefix + '.txt', data, '%15.6E', header=' '.join(col_names))

    return return_dict


def polyfit(X, Y, order, display=False, on_figs=None, return_refs=False,
            save_figs=False, save_data=False, output_prefix='polyfit'):
    """
    Fit the two equal length arrays to a polynomial with a specified order.

    Parameters
    ----------
    X: arrays
    Y: arrays
    order: int
        the polynomial order
    display: bool
        Display figures or not. Default to False.
    on_figs: list/int
        the current figure numbers to plot to, default to new figures
    return_refs: bool
        Return the axes reference(s) drawing or not. Default to False.
    save_figs: bool
        Save figures or not. Default to False.
    save_data: bool
        Save data or not. Default to False.
    output_prefix: string
        prefix string before the output files, default to 'eos_fit'

    Returns
    -------
    a dict, containing
        'coeffs': fitting coefficients
        'r-squared': value for evaluating error
        'fitted_data': a dict that has 2D array of fitted data
            easily to Pandas DataFrame by pd.DataFrame(**returned_dict['fitted_data'])
        'ax': the axes reference, if return_refs == True
    """
    popts = np.polyfit(X, Y, order)
    p = np.poly1d(popts)
    X_fit = np.linspace(sorted(X)[0], sorted(X)[-1], 1000)
    Y_fit = p(X_fit)
    Y_fit_eqlen = p(X)
    r_squared = get_r_squared(Y, Y_fit_eqlen)
    data = np.column_stack((X_fit, Y_fit))
    col_names = ['X_fit', 'Y_fit']

    return_dict = {'coeffs': p, 'r_squared': r_squared,
        'fitted_data': {'columns': col_names, 'data': data}}

    if display or save_figs or return_refs:
        initiate_figs(on_figs)
        plt.plot(X, Y, 'o')
        plt.plot(X_fit, Y_fit, '-')
        axes = {'ax': plt.gca()}
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.tight_layout()
        if save_figs:
            plt.savefig(output_prefix + '.pdf')
        return_dict = display_or_close_figs(display, return_refs, return_dict, axes)

    if save_data:
        np.savetxt(output_prefix + '.txt', data, '%15.6E', header=' '.join(col_names))

    return return_dict
