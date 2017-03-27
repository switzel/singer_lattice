#!/usr/bin/python
from pretty import *
import pickle
import sys
import argparse
from collections import Counter
from summarize import what_s_what
from base import *

triple_keys = {'moufang_2', 'split_1_2', 'K_1^2', 'N/C^2(1)'}

def dict2tuple(keys, dic):
    def get(key):
        try:
            result = dic[key]
            if key == 'aut_dm': # somewhat dirty hack
                return ((tuple(),) * len(result[0]), result[1])
            else:
                return dic[key]
        except KeyError:
            for prop, typ in what_s_what.items():
                if prop in key:
                    return typ.default
        return None
    def tuplify(lst):
        if isinstance(lst, list):
            return tuple(lst)
        return lst
    tuple_entries = zip(*[tuplify(get(key)) for key in keys if key in triple_keys])
    perm = sorted(range(3), key = lambda i: tuple_entries[i])
    def permute(entry):
        return tuple(entry[perm[i]] for i in range(3))
    return tuple(permute(tuplify(get(key))) if key in triple_keys else tuplify(get(key)) for key in keys)

def tuple2dict(keys, tup):
    return { key: val for key, val in zip(keys, tup) }

# Running as a script

parser = argparse.ArgumentParser(description = 'Condense the lattice file to a digest.')
parser.add_argument('q', type = int, help = 'the thickness parameter q')
parser.add_argument('--properties', type = str, nargs = '+', default = ['abelianization', 'commutator_abelianization', 'aut_dm', 'moufang_2', 'split_1_2', 'K_1^2', 'N/C^2(1)'], help = 'properties to take into acocunt')
parser.add_argument('--digest_name', default = 'digest_%d.pickled', help = 'file to which to write the digest, %%d will be replaced by q')
parser.add_argument('--matrices_name', default = 'digest_mats_%d.pickled', help = 'file to which to write the list of difference matrices')

if __name__ == '__main__':
    args = parser.parse_args()
    try:
        digest_name = args.digest_name % args.q
    except TypeError:
        digest_name = args.digest_name
    try:
        matrices_name = args.matrices_name % args.q
    except TypeError:
        matrices_name = args.matrices_name
    relss = unpickl(filename('diff_mats_%d.pickled', args.q))
    lattices = unpickl(filename('lattices_%d.pickled', args.q))
    count = Counter(dict2tuple(args.properties, lattice) for lattice in lattices.values())
    digest = {}
    digest_mats = []
    for val, num in count.items():
        for rels in relss:
            if dict2tuple(args.properties, lattices[rels]) == val:
                digest[rels] = tuple2dict(args.properties, val)
                digest[rels]['diff_mat'] = lattices[rels]['diff_mat']
                digest[rels]['based_diff_mat'] = lattices[rels]['based_diff_mat']
                digest[rels]['num'] = num
                digest_mats.append(rels)
                break
    digest_mats.sort()
    pickl(digest_name, digest)
    pickl(matrices_name, digest_mats)
