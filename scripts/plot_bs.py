#!/usr/bin/env python
import argparse
import re
import matplotlib as mpl

parser = argparse.ArgumentParser(description="""Plot the band structure, with
            consideration of spin-polarization. Accepts input file 'EIGENVAL', or 'vasprun.xml'.""")
parser.add_argument('-i', metavar='input', default='EIGENVAL', help="the input file name, default to 'EIGENVAL'")
parser.add_argument('--ISPIN', type=int, help="manually override ISPIN detection")
parser.add_argument('--Ef', type=float, help="manually override Ef detection")
parser.add_argument('--ylim', nargs=2, metavar=('from', 'to'), type=float, help="the range of y-axis")
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

from pydass_vasp.plotting import plot_bs
plot_bs(input_file=args.i, ISPIN=args.ISPIN, Ef=args.Ef, ylim=args.ylim)

if args.s:
    figname_sp = re.match('(.*)(\..*)', args.s).groups()
    plt.figure(1)
    plt.savefig(figname_sp[0] + '-spin-combined' + figname_sp[1])
    plt.figure(2)
    plt.savefig(figname_sp[0] + '-spin-separated' + figname_sp[1])
else:
    plt.show()
