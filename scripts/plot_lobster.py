#!/usr/bin/env python
import argparse
import re
import matplotlib as mpl

parser = argparse.ArgumentParser(description="""Plot the LOBSTER -COHP or COOP, with consideration of spin-polarization.""")
parser.add_argument('bond', type=int, help='bond number to plot')
parser.add_argument('-f', metavar='filepath', default='COHPCAR.lobster',
                    help="the filepath, default to 'COHPCAR.lobster'")
parser.add_argument('--ISPIN', type=int, help="manually override ISPIN detection")
parser.add_argument('-y', metavar='type', default='COHP', help="the input file type, default to 'COHP'")
parser.add_argument('--xlim', nargs=2, metavar=('from', 'to'), type=float, help="the range of x-axis")
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

from pydass_vasp.electronic_structure import get_lobster

get_lobster(bond=args.bond, type=args.y, filepath=args.f, ISPIN=args.ISPIN, plot=True, xlim=args.xlim, ylim=args.ylim)

if args.s:
    figname_sp = re.match('(.*)(\..*)', args.s).groups()
    plt.figure(1)
    plt.savefig(figname_sp[0] + '-spin-combined' + figname_sp[1])
    plt.figure(2)
    plt.savefig(figname_sp[0] + '-spin-separated' + figname_sp[1])
else:
    plt.show()
