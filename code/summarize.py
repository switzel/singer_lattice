#!/usr/bin/python
from base import *
from pretty import *
import pickle
import sys
import argparse

what_s_what = { 'diff_mat' : Matrix,
                'abelianization' : Homology,
                'moufang' : List(Space(Maybe, 1)),
                'split' : List(Space(Str, 1.5), raise_long=1),
                'aut_dm' : AutDiffMat,
                'K_1^2' : List(Space(Order, 1.9), raise_long = 0),
                'K_1^3' : List(Space(Order, 2.4), raise_long = 1),
                'N' : List(Space(Order, 1))}

# Running as a script

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Summarize the data collected by other scripts in the lattices file.')
    parser.add_argument('q', type = int, help = 'the thickness parameter q')
    parser.add_argument('--raw', action = 'store_true', help = 'print the raw dictionary (overrides most other options)')
    parser.add_argument('--latex', action = 'store_true', help = 'print latex code')
    parser.add_argument('--properties', type = str, nargs = '+', default = ['diff_mat', 'based_diff_mat', 'aut_dm', 'abelianization', 'commutator_abelianization', 'moufang_2', 'split_1_2', 'K_1^2', 'N/C^2(1)'], help = 'properties to print')
    parser.add_argument('--instances', type = int, nargs = '*', help = 'the indices of the exponent matrices to process; if none is given, all instances are processed')
    parser.add_argument('--number', action = 'store_true', help = 'number the lattices')
    parser.add_argument('--number_format', type = str, default = '(%d, %d)', help = 'format string for numbering')
    parser.add_argument('--matrices_filename', default = 'diff_mats_%d.pickled', help = 'file from which to read the difference matrices, %%d will be replaced by q')
    parser.add_argument('--lattice_filename', default = 'lattices_%d.pickled', help = 'filename from which to read the lattices')
    args = parser.parse_args()
    relss = unpickl(filename(args.matrices_filename, args.q))
    lattices = unpickl(filename(args.lattice_filename, args.q))
    rows = []
    n = 0
    if not args.instances:
        args.instances = xrange(len(relss))
    for i in args.instances:
        rels = relss[i]
        lattices[rels]['instance'] = i
        if args.raw:
            print(lattices[rels])
        else:
            if not args.number:
                result = []
            else:
                n += 1
                result = [args.number_format % (args.q,n)]
            for prop in args.properties:
                for name, typ in what_s_what.items():
                    if name in prop:
                        result.append(typ.str(lattices[rels].get(prop, typ.default), latex = args.latex))
                        break
                else:
                    result.append(str(lattices[rels][prop]))
            if args.latex:
                rows.append('& '.join(result))
            else:
                rows.append('\t '.join(result))
    if args.latex:
        print('\\\\\n'.join(rows))
    else:
        print('\n'.join(rows))
