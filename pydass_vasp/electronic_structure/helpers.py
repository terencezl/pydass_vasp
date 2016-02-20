import os
import re
import warnings
import matplotlib.pyplot as plt


# internal
def determine_tag_value(tag, filepath):
    dirname = os.path.dirname(filepath)
    try:
        with open(os.path.join(dirname, 'OUTCAR'), 'r') as f:
            for line in f:
                if tag in line:
                    tag_value = int(line.split()[2])
                    break
    except IOError:
        try:
            with open(os.path.join(dirname, 'INCAR'), 'r') as f:
                for line in f:
                    m = re.match(r'\s*' + tag + r'\s*=\s*(\d)\s*.*', line)
                    if m:
                        tag_value = int(m.group(1))
                        break
        except IOError:
            raise IOError("Can't determine " + tag +
                "! Either manually specify it, or provide OUTCAR or INCAR.")
    return tag_value


# internal
def figs_assert(on_figs, ISPIN, data_type):
    if ISPIN == 2:
        if data_type == 'tdos' or data_type == 'ldos' or data_type == 'lobster':
            assert on_figs is None or \
                   (isinstance(on_figs, list) and len(on_figs) == 2), \
                'The number of figures should be 2!'
    elif ISPIN == 1:
        assert on_figs is None or \
               isinstance(on_figs, int) or \
               (isinstance(on_figs, list) and len(on_figs) == 1), \
            'The number of figures should be 1!'


# internal
def initiate_figs(on_figs):
    if on_figs is None:
        plt.figure()
    elif isinstance(on_figs, int):
        plt.figure(on_figs)
    else:
        plt.figure(on_figs.pop(0))


# internal
def plot_helper_settings(axis_range, data_type):
    plt.axhline(y=0, c='k')
    plt.axvline(x=0, ls='--', c='k', alpha=0.5)
    xlim, ylim = axis_range
    if xlim and (xlim[0] != None) and (xlim[1] != None):
        plt.xlim(xlim)
    if ylim and (ylim[0] != None) and (ylim[1] != None):
        plt.ylim(ylim)
    plt.xlabel('Energy (eV)')
    if data_type == 'tdos':
        plt.ylabel('TDOS (States / Unit cell / eV)')
    elif data_type == 'ldos':
        plt.ylabel('PDOS (States / Unit cell / eV)')
    elif data_type == 'COHP':
        plt.ylabel('-pCOHP per bond')
    elif data_type == 'COOP':
        plt.ylabel('pCOOP per bond')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plt.legend(loc=0, fontsize='small')
    fig = plt.gcf()
    fig.set_tight_layout(True)
    plt.draw()
