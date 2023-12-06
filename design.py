import os
import argparse
import dataclasses
import string
import numpy as np
import pylab
import vehicle
from graph import AvlInput, plot_all

def write_avl(input_avl_fname, output_avl_fname, wing, Cref):
    Bref = wing.sections[-1].xyzLE[1] * 2
    Sref = Bref * Cref
    iSec = 0
    with open(input_avl_fname) as fin, open(output_avl_fname, 'wt') as fout:
        for line in fin:
            fout.write(line)
            if line.strip() == '!Sref     Cref     Bref':
                line = ' '.join([f'{v:8.3f}' for v in [Sref, Cref, Bref]])
                fout.write(line + '\n')
                fin.readline()
            elif line.strip() == '!     Xle      Yle      Zle    Chord     Ainc':
                sec = wing.sections[iSec]
                sec_vars = (*sec.xyzLE, sec.chord, np.rad2deg(sec.aInc))
                line = ' '.join([f'{v:8.3f}' for v in sec_vars])
                fout.write(line + '\n')
                fin.readline()
                iSec += 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='substitute',
            description='substitute new planform into AVL file')
    parser.add_argument('input_avl_fname')
    parser.add_argument('cylinders_fname')
    parser.add_argument('output_avl_fname')
    parser.add_argument('--x-scale', type=float, default=1.0)
    parser.add_argument('--yz-scale', type=float, default=1.0)
    parser.add_argument('--thick-scale', type=float, default=1.0)
    parser.add_argument('--sweep-angle', type=float, default=30.58915328)
    root_spar_x_over_c = 0.4
    tip_spar_x_over_c = 0.4

    args = parser.parse_args()
    avlInput = AvlInput(args.input_avl_fname)
    cylinders = np.loadtxt(args.cylinders_fname, skiprows=1)

    wing = avlInput.surfaces[0]
    root, tip = wing.sections[0], wing.sections[-1]
    xSparRoot = root.xyzLE[0] + root.chord * root_spar_x_over_c
    xSparTip = tip.xyzLE[0] + tip.chord * tip_spar_x_over_c

    dx = wing.sections[0].chord * root_spar_x_over_c * (args.x_scale - args.yz_scale)
    for sec in wing.sections:
        xSpar0 = xSparRoot + (xSparTip - xSparRoot) * sec.xyzLE[1] / tip.xyzLE[1]
        xLE0, ySpar0, zSpar0 = sec.xyzLE
        ySpar, zSpar = np.array([ySpar0, zSpar0]) * args.yz_scale
        xSpar = xSparRoot + np.tan(np.deg2rad(args.sweep_angle)) * ySpar
        xLE = xSpar + (xLE0 - xSpar0) * args.x_scale + dx
        sec.xyzLE = (xLE, ySpar, zSpar)
        sec.chord = sec.chord * args.x_scale
        sec.airfoil[:,1] *= args.thick_scale

    Cref = avlInput.Cref * args.x_scale
    write_avl(args.input_avl_fname, args.output_avl_fname, wing, Cref)
    plot_all(wing, cylinders, 0.4, 0.4)
    pylab.show()
