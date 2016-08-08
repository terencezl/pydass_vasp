import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cf
from ..electronic_structure.helpers import initiate_figs


# internal use
def get_r_squared(Y, Y_fit_eqlen):
    ss_resid = np.sum((Y_fit_eqlen - Y) ** 2)
    Y_avg = np.sum(Y) / len(Y)
    ss_total = np.sum((Y - Y_avg) ** 2)
    r_squared = 1 - ss_resid / ss_total
    return r_squared


def birch_murnaghan(V, V0, B0, B0_prime, E0):
    """
    3rd order Birch-Murnaghan equation of state, in the energy-volume form
    """
    V = np.array(V)
    return E0 + 9 * V0 * B0 / 16. * (
        ((V0 / V) ** (2 / 3.) - 1) ** 3 * B0_prime +
        ((V0 / V) ** (2 / 3.) - 1) ** 2 * (6 - 4 * (V0 / V) ** (2 / 3.)))


def birch_murnaghan_p(V, V0, B0, B0_prime):
    """
    3rd order Birch-Murnaghan equation of state, in the pressure-volume form
    The unit of pressure depends on the unit of the input B0.
    """
    V = np.array(V)
    return 3 * B0 / 2. * ((V0 / V) ** (7 / 3.) - (V0 / V) ** (5 / 3.)) * \
           (1 + 3 / 4. * (B0_prime - 4) * ((V0 / V) ** (2 / 3.) - 1))


def vinet(V, V0, B0, B0_prime, E0):
    """
    Vinet equation of state in the energy-volume form
    """
    V = np.array(V)
    x = (V / V0) ** (1.0 / 3)
    xi = 3.0 / 2 * (B0_prime - 1)
    return E0 + (9 * B0 * V0 / (xi ** 2) *
                 (1 + (xi * (1 - x) - 1) * np.exp(xi * (1 - x))))


def vinet_p(V, V0, B0, B0_prime):
    """
    Vinet equation of state in the pressure-volume form
    """
    V = np.array(V)
    x = (V / V0) ** (1.0 / 3)
    xi = 3.0 / 2 * (B0_prime - 1)
    return 3 * B0 / x ** 2 * (1 - x) * np.exp(xi * (1 - x))


# a decorator
def fix_B0_prime(f, B0_prime):
    if '_p' in f.__name__:
        def g(V, V0, B0, E0):
            return f(V, V0, B0, B0_prime, E0)
    elif '_p' not in f.__name__:
        def g(V, V0, B0):
            return f(V, V0, B0, B0_prime)
    else:
        raise Exception('Wrong usage!')
    return g


def eos_fit(V, Y, eos='birch_murnaghan', B0_prime=None, plot=False, on_figs=None):
    """
    Fit the volume and total energy, or pressure to the Birch-Murnaghan equation of state.

    Note: bulk modulus B0 will be returned in the unit of GPa.

    Parameters
    ----------
    V: array
        volume (Angstrom^3)
    Y: array
        total energy (eV), or pressure (GPa) if eos has '_p'
    eos: string
        chosen from ['birch_murnaghan', 'vinet']. Default to 'birch_murnaghan'
    B0_prime: float
        Keep B0_prime fixed to a given value or not. Default to None
    plot: bool
        whether to plot the data, default to False.
    on_figs: list/int
        the current figure numbers to plot to, default to new figures

    Returns
    -------
    a dict, containing
        'params': fitting parameters
        'r_squared': value for evaluating error
        'fitted_data': a dict that has 2D array of fitted data
            easily to Pandas DataFrame by pd.DataFrame(**returned_dict['fitted_data'])
    """
    V = np.array(V)
    Y = np.array(Y)

    if not B0_prime:
        if '_p' not in eos:
            initial_parameters = [V.mean(), 2.5, 4, Y.mean()]
            fit_eqn = eval(eos)
        else:
            initial_parameters = [V.mean(), 2.5, 4]
            fit_eqn = eval(eos)
    else:
        if '_p' not in eos:
            initial_parameters = [V.mean(), 2.5, Y.mean()]
            fit_eqn = fix_B0_prime(eval(eos), B0_prime)
        else:
            initial_parameters = [V.mean(), 2.5]
            fit_eqn = fix_B0_prime(eval(eos), B0_prime)

    popt, pcov = cf(fit_eqn, V, Y, initial_parameters)
    V_fit = np.linspace(sorted(V)[0], sorted(V)[-1], 1000)
    Y_fit = fit_eqn(V_fit, *popt)
    Y_fit_eqlen = fit_eqn(V, *popt)
    r_squared = get_r_squared(Y, Y_fit_eqlen)
    data = np.column_stack((V_fit, Y_fit))
    if '_p' not in eos:
        col_names = ['V_fit', 'E_fit']
    else:
        col_names = ['V_fit', 'p_fit']

    params = {'V0': popt[0], 'B0': popt[1] * 160.2}
    if not B0_prime:
        params['B0_prime'] = popt[2]
        if '_p' not in eos:
            params['E0'] = popt[3]
    else:
        if '_p' not in eos:
            params['E0'] = popt[2]

    return_dict = {'params': params, 'r_squared': r_squared,
                'fitted_data': {'columns': col_names, 'data': data}}

    if plot:
        initiate_figs(on_figs)
        plt.plot(V, Y, 'o')
        plt.plot(V_fit, Y_fit, '-')
        axes = {'ax': plt.gca()}
        plt.xlabel(r'V ($\AA^{3}$)')
        if '_p' not in eos:
            plt.ylabel('E (eV)')
        else:
            plt.ylabel('P (GPa)')
        plt.tight_layout()
        return_dict.update(axes)

    return return_dict


def polyfit(X, Y, order, plot=False, on_figs=None):
    """
    Fit the two equal length arrays to a polynomial with a specified order.

    Parameters
    ----------
    X: array
    Y: array
    order: int
        the polynomial order
    plot: bool
        plot figures or not. Default to False
    on_figs: list/int
        the current figure numbers to plot to, default to new figures

    Returns
    -------
    a dict, containing
        'coeffs': fitting coefficients
        'r_squared': value for evaluating error
        'fitted_data': a dict that has 2D array of fitted data
            easily to Pandas DataFrame by pd.DataFrame(**returned_dict['fitted_data'])
        'ax': the axes reference, if return_refs == True
    """
    X = np.array(X)
    Y = np.array(Y)

    p = np.polyfit(X, Y, order)
    p = np.poly1d(p)
    X_fit = np.linspace(sorted(X)[0], sorted(X)[-1], 1000)
    Y_fit = p(X_fit)
    Y_fit_eqlen = p(X)
    r_squared = get_r_squared(Y, Y_fit_eqlen)
    data = np.column_stack((X_fit, Y_fit))
    col_names = ['X_fit', 'Y_fit']

    return_dict = {'coeffs': p, 'r_squared': r_squared,
        'fitted_data': {'columns': col_names, 'data': data}}

    if plot:
        initiate_figs(on_figs)
        plt.plot(X, Y, 'o')
        plt.plot(X_fit, Y_fit, '-')
        axes = {'ax': plt.gca()}
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.tight_layout()
        return_dict.update(axes)

    return return_dict


def curve_fit(fit_eqn, X, Y, p0=None, sigma=None, absolute_sigma=False, plot=False, on_figs=None, **kw):
    """
    Fit the two equal length arrays to a certain function Y = f(X).

    Parameters
    ----------
    X: array
    Y: array
    plot: bool
        plot figures or not. Default to False.
    on_figs: list/int
        the current figure numbers to plot to, default to new figures

    Returns
    -------
    a dict, containing
        'params': fitting parameters
        'r_squared': value for evaluating error
        'fitted_data': a dict that has 2D array of fitted data
            easily to Pandas DataFrame by pd.DataFrame(**returned_dict['fitted_data'])
        'ax': the axes reference, if return_refs == True
    """
    X = np.array(X)
    Y = np.array(Y)
    popt, pcov = cf(fit_eqn, X, Y, p0, sigma, absolute_sigma, **kw)
    X_fit = np.linspace(sorted(X)[0], sorted(X)[-1], 1000)
    Y_fit = fit_eqn(X_fit, *popt)
    Y_fit_eqlen = fit_eqn(X, *popt)
    r_squared = get_r_squared(Y, Y_fit_eqlen)
    data = np.column_stack((X_fit, Y_fit))
    col_names = ['X_fit', 'Y_fit']

    return_dict = {'params': popt, 'r_squared': r_squared,
        'fitted_data': {'columns': col_names, 'data': data}}

    if plot:
        initiate_figs(on_figs)
        plt.plot(X, Y, 'o')
        plt.plot(X_fit, Y_fit, '-')
        axes = {'ax': plt.gca()}
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.tight_layout()
        return_dict.update(axes)

    return return_dict
