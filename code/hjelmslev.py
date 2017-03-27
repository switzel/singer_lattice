#!/usr/bin/python
import pickle
import argparse
import itertools
import sys
import os
import collections
import re
import pretty
from base import pickl, unpickl, factorize
from difference_matrices import DiffSetFactory

Adjacencies = collections.namedtuple('Adjacencies', 'below above adjacent')

STABILIZER_GAP = 'stabilizers.gap'
HJELMSLEV_GAP = 'hjelmslev.gap'
LIFTS_GAP = 'lifts.gap'

Maybe = None

def diff_sets(diff_mat, typ = set):
    d = [typ(), typ(), typ()]
    for row in diff_mat:
        for i, entry in enumerate(row):
            if typ == set:
                d[i].add(entry)
            elif typ == dict:
                d[i][entry] = row
            else:
                raise ValueError('The only allowed types for diff_set are dict and set, got %s.' % str(typ))
    return d

def rotate(diff_mat, i):
    return tuple( row[i:] + row[:i] for row in diff_mat)

class Hjelmslev:
    def __init__(self):
        self.points = []
        self.lines = []
        self.spanning_points = []

    def init(self, q, diff_mat, vi, maxlevel = 2, prime = False):
        def empty_adjacencies():
            return Adjacencies(below = [tuple() for i in range(maxlevel+1)],
                               above = [set() for i in range(maxlevel+1)],
                               adjacent = set())
        self.points = [collections.defaultdict(empty_adjacencies) for i in range(maxlevel + 1)]
        self.lines = [collections.defaultdict(empty_adjacencies) for i in range(maxlevel + 1)]
        if vi != 0:
            diff_mat = rotate(diff_mat, vi)
        self.q = q
        self.delta = q**2 + q + 1
        self.diff_mat = diff_mat
        self.prime = prime
        d = diff_sets(diff_mat, typ = dict)
        nd = [set((-x) % self.delta for x in ds.keys()) for ds in d]
        flags = [set() for i in range(maxlevel + 1)]
        for a0 in range(self.delta):
            for e in diff_mat:
                b0 = (a0 + e[0]) % self.delta
                flags[1].add(((b0,), (a0,)))
                if maxlevel <= 1:
                    continue
                for a2 in set(range(self.delta)) - nd[2]:
                    for f in diff_mat:
                        try:
                            h = d[2][(f[2] - e[2] - a2) % self.delta]
                        except KeyError:
                            continue
                        for g in diff_mat:
                            b1 = (e[1] - f[1] + g[1]) % self.delta
                            if not b1 in set(range(self.delta)) - set(d[1].keys()):
                                continue
                            flags[2].add(((b0,b1),(a0,a2)))
                            if maxlevel <= 2:
                                continue
                            for a1 in set(range(self.delta)) - nd[1]:
                                for k in diff_mat:
                                    for l in diff_mat:
                                        if not (l[1] - h[1] - a1) % self.delta in d[1].keys():
                                            continue
                                        if not (-g[0] + k[0] + f[0] - h[0] + l[0]) % self.delta in d[0].keys():
                                            continue
                                        for m in d[2].keys():
                                            b2 = (g[2] - k[2] + m) % self.delta
                                            if b2 in set(range(self.delta)) - set(d[2].keys()):
                                                flags[3].add(((b0,b1,b2),(a0,a2,a1)))
        for level, level_flags in enumerate(flags):
            for point, line in level_flags:
                self.points[level][point].adjacent.add(line)
                self.lines[level][line].adjacent.add(point)
        for pt in self.points[maxlevel]:
            for level in range(1, maxlevel+1):
                for sublevel in range(1,level+1):
                    self.points[level][pt[:level]].below[sublevel]=pt[:sublevel]
                    self.points[sublevel][pt[:sublevel]].above[level].add(pt[:level])
        for ln in self.lines[maxlevel]:
            for level in range(1, maxlevel+1):
                for sublevel in range(1,level+1):
                    self.lines[level][ln[:level]].below[sublevel]=ln[:sublevel]
                    self.lines[sublevel][ln[:sublevel]].above[level].add(ln[:level])

    def init_spanning_points(self, maxlevel = 1):
        if maxlevel > 1:
            raise NotImplemented, 'Spanning points for higher level'
        if len(self.spanning_points) >= maxlevel + 1:
            return
        self.spanning_points = [[] for i in range(maxlevel+1)]
        chosen = []
        span_pts, span_lns = dict(), dict()
        possible = set(self.points[maxlevel].keys())
        preferred = set(self.points[maxlevel].keys())
        while (possible):
            if preferred:
                pt = next(iter(preferred))
            else:
                pt = next(iter(possible))
            chosen.append(pt)
            span_pts[pt] = None
            span_pts, span_lns = self.span(span_pts, span_lns, maxlevel)
            possible -= set(span_pts.keys())
            preferred &= possible
            for ln in span_lns:
                preferred -= self.lines[maxlevel][ln].adjacent
        self.spanning_points[maxlevel] = chosen

    def graph(self, level):
        result = []
        for pt, adjs in self.points[level].items():
            for ln in adjs.adjacent:
                result.append((pt, ln))
        return result

    def gap_graph(self, level, lower_level = 1, out_file = sys.stdout, variable_postfix = ''):
        pt_lst = list(self.points[level].keys())
        ln_lst = list(self.lines[level].keys())
        pt_dict = { pt: pti+1 for pti, pt in enumerate(pt_lst) }
        ln_dict = { ln: lni+len(pt_lst)+1 for lni, ln in enumerate(ln_lst) }
        def pt_str(pt, ptln):
            return '["%s", %s]' % (ptln, ', '.join(str(n) for n in pt))
        out_file.write('q := %d;\n' % self.q);
        out_file.write('vertices%s := [%s, %s];\n' %
                           (variable_postfix,
                            ', '.join(pt_str(pt, 'p') for pt in self.points[level].keys() if pt[0] == 0),
                            ', '.join(pt_str(ln, 'l') for ln in self.lines[level].keys() if ln[0] == 0)))
        out_file.write('edges%s := [%s];\n' %
                           (variable_postfix,
                            ', '.join('[%s, %s]' % (pt_str(pt, 'p'), pt_str(ln, 'l'))
                                      for pt, adj in self.points[level].items() if pt[0] == 0
                                      for ln in adj.adjacent)))

    def metapost(self, level):
        assert level == 2
        i = 0
        for pt in self.points[1][(0,)].above[level]:
            for ln in self.points[level][pt].adjacent:
                print('edges[%d] := (%s, %s);' % (i, ','.join(str(d) for d in pt[1:]), ','.join(str(d) for d in ln)))
                i += 1

    def incidence_geometry(self, level):
        lns = { line: adj.adjacent for line, adj in hj.lines[level].items() }
        return lns

    def span(self, pts, lns, level, image_level = None):
        if image_level == None:
            image = lambda ptln, pt1, pt2 : None
        else:
            def image(ptln, obj1, obj2):
                if ptln == 'pt':
                    lines = self.lines[image_level]
                else:
                    lines = self.points[image_level]
                joint = lines[obj1].adjacent & lines[obj2].adjacent
                assert len(joint) == 1, '(Dual) points not on one (dual) line.'
                return next(iter(joint))
        span_pts, span_lns = pts, lns
        new_pts, new_lns = dict(pts), dict(lns)
        while new_pts or new_lns:
            for ptln, points, lines, span_points, span_lines, new_points, new_lines in \
                (('ln', self.points[level], self.lines[level], span_pts, span_lns, new_pts, new_lns),
                 ('pt', self.lines[level], self.points[level], span_lns, span_pts, new_lns, new_pts)):
                for (pt1, im1), (pt2, im2) in itertools.product(new_points.items(), span_points.items()):
                    if pt1 in new_points and pt2 in new_points and pt1 > pt2:
                        continue
                    if points[pt1].below[1] == points[pt2].below[1]:
                        continue
                    joint = points[pt1].adjacent & points[pt2].adjacent
                    assert len(joint) == 1, '(Dual) points not on one (dual) line'
                    new_ln = next(iter(joint))
                    ln_img = image(ptln, im1, im2)
                    if new_ln in span_lines:
                        if span_lines[new_ln] != ln_img:
                            return dict(), dict()
                    else:
                        new_lines[new_ln] = ln_img
                        span_lines[new_ln] = ln_img
                new_points.clear()
        return span_pts, span_lns


    def splits(self, level, lower_level = 1, equivariant_only = False):
        global projected_pts, projected_lns
        result = []
        self.init_spanning_points(lower_level)
        for upper_pts in itertools.product(*[self.points[lower_level][pt].above[level]
                                             for pt in self.spanning_points[lower_level]]):
            pts = {pt : upper_pt for pt, upper_pt in
                        zip(self.spanning_points[lower_level], upper_pts)}
            span_pts, span_lns = self.span(pts, dict(), lower_level, level)
            if len(span_pts) == len(self.points[lower_level]):
                result.append((span_pts, span_lns))
        return result

    def chambers_of_some_apartment(self, level):
        self.init_spanning_points(1)
        lifts = [next(iter(self.points[1][pt].above[level])) for pt in self.spanning_points[1][:3]]
        for pta, ptb in itertools.permutations(lifts, 2):
            yield (pta, next(iter(self.points[level][pta].adjacent & self.points[level][ptb].adjacent)))
        return

    def enough_roots(self, level):
        for pt, ln in self.chambers_of_some_apartment(level = level):
            for bdry_ln in self.points[level][pt].adjacent:
                if self.lines[level][bdry_ln].below[1] != self.lines[level][ln].below[1]:
                    yield (bdry_ln, pt, ln)
        return

    def partial_root_group_elts(self, level, bdry_ln, pt, ln):
        fix_root_pt_map = {opt : opt for opt in self.lines[level][ln].adjacent}
        fix_root_ln_map = {oln : oln for oln in self.points[level][pt].adjacent}
        other_pts = set(other_pt for other_pt in self.lines[level][bdry_ln].adjacent
                        if not self.points[level][pt].below[1] == self.points[level][other_pt].below[1])
        base_pt = other_pts.pop()
        for other_pt in other_pts:
            pt_map, ln_map = dict(fix_root_pt_map), dict(fix_root_ln_map)
            pt_map[base_pt] = other_pt
            yield pt_map, ln_map
        return
    
    def check_root_group(self, level, bdry_ln, pt, ln):
        for pt_map, ln_map in self.partial_root_group_elts(level, bdry_ln, pt, ln):
            pt_map, ln_map = self.span(pt_map, ln_map, level, level)
            if len(pt_map) < len(self.points[level]):
                return False
        return True

    def moufang(self, level):
        for bdry_ln, pt, ln in self.enough_roots(level):
            if not self.check_root_group(level, bdry_ln, pt, ln):
                return False
        return True
                
    def stabilizers(self, level, basename, everything = False, skip_gap = False):
        gap_filename = 'gap_graph/%s.gap' % basename
        out_filename = 'gap_graph/%s.out' % basename
        ds = diff_sets(self.diff_mat)[0]
        dsf = DiffSetFactory(self.q)
        m, s = dsf.aut_diff_set(ds, frobenius_only = True)[1]
        def encode_permutation(pimages, limages):
            return '[%s, %s]' % (', '.join('[["p", %d], ["p", %d]]' % (x, ximg)
                                           for x, ximg in enumerate(pimages)),
                                 ', '.join('[["l", %d], ["l", %d]]' % (x, ximg)
                                           for x, ximg in enumerate(limages)))
        pimages = tuple((m*x + s) % self.delta for x in range(self.delta))
        limages = tuple(m*x % self.delta for x in range(self.delta))
        if not os.path.exists('gap_graph'):
            os.mkdir('gap_graph')
        if not skip_gap:
            with open(gap_filename, 'w') as gap_file:
                self.gap_graph(level, out_file = gap_file)
                gap_file.write('everything := %s;\n' % str(everything).lower());
                gap_file.write('gencoded := %s;\n' % (encode_permutation(pimages, limages)))
                with open(HJELMSLEV_GAP, 'r') as stab_file:
                    for line in stab_file:
                        gap_file.write(line)
                with open(STABILIZER_GAP, 'r') as stab_file:
                    for line in stab_file:
                        gap_file.write(line)
                gap_file.write('QUIT_GAP();\n')
            os.system('gap %s > %s' % (gap_filename, out_filename))
        stabilizers = {}
        with open(out_filename, 'r') as out_file:
            content = out_file.read()
            matches = re.finditer('\n([^\n ]+) has order (\d+)', content)
            for match in matches:
                stabilizers[match.group(1)] = int(match.group(2))
        return stabilizers

def dict_list(lst):
    result = collections.defaultdict(lambda : [None for l in lst])
    for i, d in enumerate(lst):
        for key, val in d.items():
            result[key][i] = val
    return dict(result)

def update_list_dict(d, upd):
    for key, val in upd.items():
        if key not in d:
            d[key] = val
        else:
            d[key] = [val or oval  for (oval, val) in zip(d[key], val)]

# Running as a script

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Produce the Hjelmslev plane around a vertex and perform checks.')
    parser.add_argument('q', type = int, help = 'the thickness parameter q')
    parser.add_argument('--instances', type = int, nargs = '*', help = 'the indices of the exponent matrices to process as output by difference_set_equivalence.py; if none is given, all indices are processed')
    parser.add_argument('--latex', action = 'store_true', help = 'whether to output latex code')
    parser.add_argument('--tentative', type = int, default = -1, help = 'for the Moufang check and the projection check: give up after that many steps; in this case the possible answers are True, False and Maybe')
    parser.add_argument('--indices', type = int, nargs = '+', default = [ 0, 1, 2 ], help = 'the indices to check')
    parser.add_argument('--level', type = int, default = 2, help = 'the Hjelmslev plane of which level to consider')
    parser.add_argument('--split', action = 'store_true', help = 'check in how many ways the projection to projective plane splits')
    parser.add_argument('--stabilizers', action = 'store_true', help = 'whether to investigate stabilizers')
    parser.add_argument('--all_stabilizers', action = 'store_true', help = 'whether to even compute stabilizers that are not necessarily needed')
    parser.add_argument('--moufang', action = 'store_true', help = 'check whether the Hjelmslev plane satisfies the Moufang condition')
    parser.add_argument('--dump', action = 'store_true', help = 'print the chambers of the Hjelmslev plane')
    parser.add_argument('--lower_level', type = int, default = 1, help = 'the Hjelmslev plane of which lower level to consider')
    parser.add_argument('--write_lattices', action = 'store_true', help = 'whether to (over)write data in the lattice file')
    parser.add_argument('--write_individual', action = 'store_true', help = 'whether to write data in an individual file')
    parser.add_argument('--gap_graph', action = 'store_true', help = 'whether to print the gap graph')
    parser.add_argument('--metapost', action = 'store_true', help = 'whether to print the metapost graph')
    parser.add_argument('--skip_gap', action = 'store_true', help = 'do not actually run gap')
    args = parser.parse_args()
    diffmats = tuple(((0,0,0),) + rels for rels in unpickl('diff_mats_%d.pickled' % args.q))
    if not args.instances:
        args.instances = range(len(diffmats))
    lattices = unpickl('lattices_%d.pickled' % args.q)
    for i in args.instances:
        diffmat = diffmats[i]
        moufang = [None, None, None]
        split = [None, None, None]
        stabs = [{}, {}, {}]
        moufang = [None, None, None]
        for vi in args.indices:
            hj = Hjelmslev()
            hj.init(args.q, diffmat, vi, maxlevel = args.level)
            if args.split:
                split[vi]=len(hj.splits(args.level, args.lower_level))
            if args.moufang:
                moufang[vi] = hj.moufang(args.level)
            if args.stabilizers or args.all_stabilizers:
                stabs[vi] = hj.stabilizers(args.level, 'stabilizers_%d_%d_%d' % (args.q, i, vi), everything = args.all_stabilizers, skip_gap = args.skip_gap)
            if args.gap_graph:
                hj.gap_graph(args.level, args.lower_level, variable_postfix = '_%d_%d_%d'  % (args.q, i, vi))
            if args.metapost:
                hj.metapost(args.level)
            if args.dump:
                print(hj.graph(args.level))
        old_split = lattices[diffmat[1:]].get('split_%d_%d' % (args.lower_level, args.level), [None, None, None])
        lattices[diffmat[1:]]['split_%d_%d' % (args.lower_level, args.level)] = [max(old, new) for old, new in zip(old_split, split)]
        old_moufang = lattices[diffmat[1:]].get('moufang_%d' % args.level, [None, None, None])
        lattices[diffmat[1:]]['moufang_%d' % args.level] = [max(old, new) for old, new in zip(old_moufang, moufang)]
        stab_dict = dict_list(stabs)
        update_list_dict(lattices[diffmat[1:]], stab_dict)
        if args.gap_graph or args.metapost or args.dump:
            continue
        conn = ' & ' if args.latex else ' \t '
        print(conn.join([str(i),
                          pretty.Matrix.str(diffmat, latex = args.latex),
                          pretty.List(pretty.Space(pretty.Maybe,1)).str(moufang, latex = args.latex),
                          pretty.List(pretty.Space(pretty.Str,1.5)).str(split, latex = args.latex),
                          conn.join(pretty.List(pretty.Order).str(e, latex = args.latex) for k, e in sorted(stab_dict.items()))]))
        if args.write_individual:
            pickl('lattices_%d_%d.pickled' % (args.q, i), {diffmat[1:] : lattices[diffmat[1:]]})
    if args.write_lattices:
        pickl('lattices_%d.pickled' % args.q, lattices)
