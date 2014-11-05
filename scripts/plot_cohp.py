#!/usr/bin/env python
import argparse
import matplotlib

parser = argparse.ArgumentParser(description="""Plot the -COHP, with consideration of spin-polarization.""")
parser.add_argument('bond', type=int, help='bond number to plot')
parser.add_argument('-a', '--axis-range', type=eval, help='''the x and y range of axis in the form of
            '[Xmin,Xmax,Ymin,Ymax]'. If ISPIN=2, this option specifies the combined spin.''')
parser.add_argument('--ISPIN', type=int, help="manually override ISPIN detection")
parser.add_argument('-p', '--display', action="store_true", help="display figure but not save it.")
parser.add_argument('--save-data', action="store_true", help="save extracted COHPCAR data or not.")
parser.add_argument('-i', '--input', metavar='COHPCAR', default='COHPCAR.lobster', help="the input COHPCAR file name")
parser.add_argument('-o', '--output-prefix', default='COHP', help="the output files' prefix")
args = parser.parse_args()
if args.display:
    display = True
    save_figs = False
else:
    display = False
    save_figs = True

if not args.display:
    matplotlib.use('Agg')
else:
    if matplotlib.get_backend() == 'MacOSX':
        try:
            matplotlib.use('TkAgg')
            print "Switched from MacOSX backend to TkAgg."
        except ValueError:
            print "Can't switched to TkAgg. Keep using MacOSX backend."
import matplotlib.pyplot as plt
try:
    plt.style.use('ggplot')
    print "Changed style to ggplot, just prettier."
except AttributeError:
    print "If you upgrade to matplotlib 1.4 and I will change the style to ggplot, just prettier."

from pydass_vasp.plotting import plot_cohp
plot_cohp(bond=args.bond, axis_range=args.axis_range, ISPIN=args.ISPIN, input_file=args.input, display=display, return_refs=False,
    save_figs=save_figs, save_data=args.save_data, output_prefix=args.output_prefix)