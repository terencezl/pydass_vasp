#!/usr/bin/env python
import argparse
import re
import matplotlib as mpl

parser = argparse.ArgumentParser(description="""Plot the band structure, with
            consideration of spin-polarization. Accepts input file 'EIGENVAL', or 'vasprun.xml'.""")
parser.add_argument('-f', metavar='filepath', default='EIGENVAL', help="the input filepath, default to 'EIGENVAL'")
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

from pydass_vasp.electronic_structure import get_bs

get_bs(filepath=args.f, ISPIN=args.ISPIN, Ef=args.Ef, plot=True, ylim=args.ylim)

if args.s:
    plt.savefig(args.s)
else:
    plt.show()
