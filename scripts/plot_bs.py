#!/usr/bin/env python
import argparse
import matplotlib

parser = argparse.ArgumentParser(description="""Plot the band structure, with
            consideration of spin-polarization. Accepts input file 'EIGENVAL', or 'vasprun.xml'.""")
parser.add_argument('-i', '--input', metavar='EIGENVAL', default='EIGENVAL', help="the input EIGENVAL file name")
parser.add_argument('--ISPIN', type=int, help="manually override ISPIN detection")
parser.add_argument('--Ef', type=float, help="manually override Ef detection")
parser.add_argument('--ylim', type=eval, help="the range of y-axis in the form of '[min,max]'.")
parser.add_argument('-p', '--display', action="store_true", help="display figure but not save it.")
parser.add_argument('-s', '--save-data', action="store_true", help="save extracted EIGENVAL data or not.")
parser.add_argument('-t', '--style', help="use a plotting style for matplotlib >= 1.4")
parser.add_argument('-o', '--output-prefix', default='BS', help="the output files' prefix")
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

from pydass_vasp.plotting import plot_bs
plot_bs(input_file=args.input, ISPIN=args.ISPIN, Ef=args.Ef, ylim=args.ylim, display=display, return_refs=False,
    save_figs=save_figs, save_data=args.save_data, output_prefix=args.output_prefix)
