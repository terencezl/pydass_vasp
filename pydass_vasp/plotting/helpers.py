import re
import matplotlib.pyplot as plt


# internal
def determine_tag_value(tag):
    try:
        with open('OUTCAR', 'r') as f:
            for line in f:
                if tag in line:
                    tag_value = int(line.split()[2])
    except IOError:
        try:
            with open('INCAR', 'r') as f:
                for line in f:
                    m = re.match(r'\s*' + tag + r'\s*=\s*(\d)\s*.*', line)
                    if m:
                        tag_value = int(m.group(1))
        except IOError:
            raise IOError("Can't determine " + tag +
                "! Either manually specify it, or provide OUTCAR or INCAR.")
    return tag_value


# internal
def plot_helper_figs_assert(on_figs, ISPIN, data_type):
    if ISPIN == 2:
        if data_type == 'tdos' or data_type == 'ldos' or data_type == 'cohp':
            assert on_figs is None or (isinstance(on_figs, list) and len(on_figs) == 2), \
                'The number of figures should be 2!'
        # elif data_type == 'bs':
        #     assert on_figs is None or (isinstance(on_figs, list) and len(on_figs) == 3), \
        #         'The number of figures should be 3!'
    elif ISPIN == 1:
        assert on_figs is None or (isinstance(on_figs, str)) or \
        (isinstance(on_figs, list) and len(on_figs) == 1), 'The number of figures should be 1!'


# internal
def plot_helper_figs(on_figs):
    if on_figs is None:
        plt.figure()
    else:
        plt.figure(on_figs.pop(0))