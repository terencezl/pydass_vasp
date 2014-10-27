#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
if not matplotlib.is_interactive():
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
try:
    plt.style.use('ggplot')
except AttributeError:
    print "If you upgrade to matplotlib 1.4 and I will change the style to ggplot, just prettie    r."
from scipy.optimize import curve_fit

# internal use
def B_M_eqn(V, V0, B0, B0_prime):
    return (3*B0/2. * ((V0/V)**(7/3.) - (V0/V)**(5/3.)) * (1 + 3/4. * (B0_prime - 4) * ((V0/V)**(2/3.) - 1)))

# internal use
def B_M_fit(V, P):
    coeffs, pcov = curve_fit(B_M_eqn, V, P, [[V.mean(), 2.5, 4]])
    x_fit = np.linspace(sorted(V)[0], sorted(V)[-1], 1000)
    y_fit = B_M_eqn(x_fit, *coeffs)
    y_fit_eqlen = B_M_eqn(V, *coeffs)
    ss_resid = np.sum((y_fit_eqlen - P)**2)
    y_avg = np.sum(P)/len(P)
    ss_total = np.sum((P-y_avg)**2)
    r_squared = 1 - ss_resid/ss_total
    return (coeffs, r_squared, x_fit, y_fit)

# internal use
def B_M_eqn_fixB0prime(B0_prime):
    def B_M_eqn(V, V0, B0):
        return (3*B0/2. * ((V0/V)**(7/3.) - (V0/V)**(5/3.)) * (1 + 3/4. * (B0_prime - 4) * ((V0/V)**(2/3.) - 1)))
    return B_M_eqn

# internal use
def B_M_fit_fixB0prime(V, P, B0_prime):
    coeffs, pcov = curve_fit(B_M_eqn_fixB0prime(B0_prime), V, P, [[V.mean(), 2.5]])
    x_fit = np.linspace(sorted(V)[0], sorted(V)[-1], 1000)
    y_fit = B_M_eqn_fixB0prime(B0_prime)(x_fit, *coeffs)
    y_fit_eqlen = B_M_eqn_fixB0prime(B0_prime)(V, *coeffs)
    ss_resid = np.sum((y_fit_eqlen - P)**2)
    y_avg = np.sum(P)/len(P)
    ss_total = np.sum((P-y_avg)**2)
    r_squared = 1 - ss_resid/ss_total
    return (coeffs, r_squared, x_fit, y_fit)

def eos_pv(V, P, fix_B0_prime=False):
    """
    Takes volume (Angstrom^3) and pressure (GPa) in arrays.
    Optional: fix_B0_prime. Default to False, or give it a value, like 4.
    Returns dictionary: parameters, r-squared value and fitted_arrays.
    """
    V = np.array(V)
    P = np.array(P)
    if fix_B0_prime:
        (coeffs, r_squared, V_fit, P_fit) = B_M_fit_fixB0prime(V, P, fix_B0_prime)
    else:
        (coeffs, r_squared, V_fit, P_fit) = B_M_fit(V, P)

    parameters = {'V0': coeffs[0], 'B0': coeffs[1]}
    if not fix_B0_prime:
        parameters['B0_prime'] = coeffs[2]

    plt.plot(V, P, 'o')
    plt.plot(V_fit, P_fit, '-')
    plt.xlabel(r'V ($\AA^{3}$)')
    plt.ylabel('P (GPa)')
    plt.legend()
    plt.tight_layout()
    return {'parameters': parameters, 'r_squared': r_squared,
            'fitted_arrays': np.column_stack((V_fit, P_fit))}


def polyfit(x, y, order):
    """
    Takes x and y arrays, and the polynomial order.
    Returns dictionary: parameters, r-squared value and fitted_arrays.
    """
    popts = np.polyfit(x, y, order, full=True)
    p = np.poly1d(popts[0])
    x_fit = np.linspace(sorted(x)[0], sorted(x)[-1], 1000)
    y_fit = p(x_fit)
    y_avg = np.sum(y) / len(y)
    ss_total = np.sum((y - y_avg) ** 2)
    r_squared = 1 - popts[1][0] / ss_total
    return {'coeffs': p, 'r_squared': r_squared, 'fitted_arrays': np.column_stack((x_fit, y_fit))}
