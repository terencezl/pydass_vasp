#!/usr/bin/env python
import argparse
import matplotlib

parser = argparse.ArgumentParser(description="""Plot the local projected density of states, with
            consideration of spin-polarization. Accepts input file 'DOSCAR', or 'vasprun.xml'.""")
parser.add_argument('atom', type=int, help='atom to plot')
parser.add_argument('-a', '--axis-range', type=eval, help='''the x and y range of axis in the form of
            '[Xmin,Xmax,Ymin,Ymax]'. If ISPIN=2, this option specifies the combined spin.''')
parser.add_argument('--ISPIN', type=int, help="manually override ISPIN detection")
parser.add_argument('--LORBIT', type=int, help="manually override LORBIT detection")
parser.add_argument('-p', '--display', action="store_true", help="display figure but not save it.")
parser.add_argument('-s', '--save-data', action="store_true", help="save extracted EIGENVAL data or not.")
parser.add_argument('-t', '--style', help="use a plotting style for matplotlib >= 1.4")
parser.add_argument('-i', '--input', metavar='DOSCAR', default='DOSCAR', help="the input DOSCAR file name")
parser.add_argument('-o', '--output-prefix', default='LDOS', help="the output files' prefix")
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
            print("Switched from MacOSX backend to TkAgg.")
        except ValueError:
            print("Can't switched to TkAgg. Keep using MacOSX backend.")
import matplotlib.pyplot as plt
if args.style:
    try:
        plt.style.use(args.style)
    except AttributeError:
        print("Style is only supported for matplotlib >= 1.4.")

from pydass_vasp.plotting import plot_ldos
plot_ldos(atom=args.atom, axis_range=args.axis_range, ISPIN=args.ISPIN, LORBIT=args.LORBIT, input_file=args.input, display=display, return_refs=False,
    save_figs=save_figs, save_data=args.save_data, output_prefix=args.output_prefix)