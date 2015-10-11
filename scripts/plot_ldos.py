#!/usr/bin/env python
import argparse
import re
import matplotlib as mpl

parser = argparse.ArgumentParser(description="""Plot the local projected density of states, with
            consideration of spin-polarization. Accepts input file 'DOSCAR', or 'vasprun.xml'.""")
parser.add_argument('atom', type=int, help='atom to plot')
parser.add_argument('-i', metavar='input', default='DOSCAR', help="the input file name, default to 'DOSCAR'")
parser.add_argument('--ISPIN', type=int, help="manually override ISPIN detection")
parser.add_argument('--LORBIT', type=int, help="manually override LORBIT detection")
parser.add_argument('--xlim', nargs=2, metavar=('from', 'to'), type=float, help="the range of x-axis")
parser.add_argument('--ylim_upper', type=float, help="the upper limit of y-axis (of the spin-combined plot if ISPIN == 2)")
parser.add_argument('-s', metavar='savefig_name', help="if specified, save to file but not display")
parser.add_argument('-t', help="use a plotting style for matplotlib >= 1.4")
args = parser.parse_args()
if args.s:
    mpl.use('Agg')
import matplotlib.pyplot as plt

if args.t:
    try:
        plt.style.use(args.t)
    except AttributeError:
        print("Style is only supported for matplotlib >= 1.4.")

from pydass_vasp.plotting import plot_ldos
plot_ldos(atom=args.atom, input_file=args.i, ISPIN=args.ISPIN, LORBIT=args.LORBIT, xlim=args.xlim, ylim_upper=args.ylim_upper)

if args.s:
    figname_sp = re.match('(.*)(\..*)', args.s).groups()
    plt.figure(1)
    plt.savefig(figname_sp[0] + '-spin-combined' + figname_sp[1])
    plt.figure(2)
    plt.savefig(figname_sp[0] + '-spin-separated' + figname_sp[1])
else:
    plt.show()
