import numpy as np
import matplotlib
if not matplotlib.is_interactive():
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
try:
    plt.style.use('ggplot')
except AttributeError:
    print "If you upgrade to matplotlib 1.4 and I will change the style to ggplot, just prettier."
from scipy.optimize import curve_fit


def get_r_squared(Y, Y_fit_eqlen):
    ss_resid = np.sum((Y_fit_eqlen - Y) ** 2)
    Y_avg = np.sum(Y) / len(Y)
    ss_total = np.sum((Y - Y_avg) ** 2)
    r_squared = 1 - ss_resid / ss_total
    return r_squared

################################# EOS ################################

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


def eos_fit(V, Y, p_v=False, fix_B0_prime=False, plot=False, new_fig=False):
    """
    Takes volume (Angstrom^3) and total energy (eV) or pressure (GPa) in arrays.
    Optional: fix_B0_prime. Default to False, or give it a value, like 4.
    Returns dictionary: parameters, r-squared value and fitted_arrays.
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

    return {'parameters': parameters, 'r_squared': r_squared,
            'fitted_arrays': {'V': V_fit, 'Y': Y_fit}}


################################# Polyfit #################################

def polyfit(X, Y, order, plot=False, new_fig=False):
    """
    Takes X and Y arrays, and the polynomial order.
    Returns dictionary: parameters, r-squared value and fitted_arrays.
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

    return {'coeffs': p, 'r_squared': r_squared, 'fitted_arrays': {'X': X_fit, 'Y': Y_fit}}
