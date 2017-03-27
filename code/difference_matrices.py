'''
This is a collection of functions to study difference sets up to multiplicative
and additive equivalence. Further functions are related to so-called difference
matrices that parametrize panel regular lattices.

When run as a script, the input is a prime power q and the output consists of
representatives of all (Desarguesian) difference matrices of order q up to
equivalence. The number of difference matrices grows very fast so one probably
doesn't want to run it with q > 5.

Written by Stefan Witzel in 2015/16.
'''

from numpy import array, empty_like
from base import pickl, unpickl, factorize
import numpy
from fractions import gcd, Fraction
from itertools import permutations, combinations_with_replacement, product
from math import factorial
import argparse
import pickle
import sys
import num_th
from collections import defaultdict

def order_3_elts(k):
    if k < 3:
        return 0
    return order_3_elts(k-1) + (1 + order_3_elts(k-3))*(k-1)*(k-2)

def bound_diff_mat_classes(p,r):
    k = p**r + 1
    return Fraction(1,6*27*r**3)*(factorial(k)**2 + 3*factorial(k) + 2*(order_3_elts(k)+1))

def eq_classes(reps, rel, verbose = True):
    result = set()
    while reps:
        rep = next(iter(reps))
        cls = rel(rep)
        result.add(frozenset(cls))
        reps -= cls
        if verbose:
            print('Representatives left: %d.' % len(reps))
    return result

# Code related to difference sets

diff_sets = {
    2 : (0, 1, 3), 
    3 : (0, 1, 4, 6), 
    4 : (0, 1, 4, 14, 16), 
    5 : (0, 1, 6, 18, 22, 29), 
    7 : (0, 1, 5, 7, 17, 35, 38, 49), 
    8 : (0, 1, 17, 39, 41, 44, 48, 54, 62), 
    9 : (0, 1, 3, 9, 27, 49, 56, 61, 77, 81), 
    11 : (0, 1, 8, 14, 30, 45, 47, 56, 66, 106, 109, 129), 
    13 : (0, 1, 7, 15, 37, 41, 46, 79, 100, 103, 123, 155, 171, 173), 
    16 : (0, 1, 24, 35, 38, 40, 53, 86, 108, 114, 118, 135, 144, 185, 210, 254, 266), 
    17 : (0, 1, 52, 71, 78, 110, 113, 128, 200, 212, 214, 234, 243, 259, 267, 280, 297, 303), 
    19 : (0, 1, 4, 52, 62, 64, 84, 131, 154, 162, 180, 187, 273, 294, 300, 309, 328, 339, 344, 368), 
    23 : (0, 1, 20, 107, 119, 125, 162, 167, 189, 233, 246, 344, 369, 379, 390, 405, 419, 457, 464, 473, 481, 520, 522, 550), 
    25 : (0, 1, 34, 58, 77, 129, 193, 220, 281, 295, 327, 377, 387, 389, 395, 400, 417, 426, 442, 498, 527, 565, 568, 572, 616, 631), 
    27 : (0, 1, 30, 48, 76, 150, 175, 208, 219, 246, 272, 287, 289, 308, 311, 367, 398, 402, 412, 418, 452, 475, 484, 536, 603, 652, 745, 750), 
    31 : (0, 1, 33, 86, 90, 132, 148, 168, 191, 213, 241, 251, 260, 262, 265, 446, 490, 507, 586, 615, 650, 656, 663, 690, 774, 792, 800, 872, 887, 926, 938, 963), 
    32 : (0, 12, 92, 96, 151, 154, 174, 175, 226, 286, 300, 335, 343, 354, 453, 459, 499, 501, 530, 540, 566, 596, 603, 630, 718, 723, 736, 751, 768, 812, 821, 837, 859), 
    37 : (0, 1, 34, 72, 84, 150, 187, 213, 235, 244, 403, 420, 444, 543, 641, 649, 717, 736, 749, 818, 829, 869, 931, 966, 973, 976, 991, 996, 1012, 1040, 1087, 1101, 1130, 1157, 1276, 1349, 1353, 1355), 
    41 : (0, 1, 61, 135, 173, 190, 197, 219, 244, 276, 297, 310, 320, 340, 390, 471, 489, 508, 566, 581, 597, 648, 653, 736, 775, 801, 834, 869, 1107, 1135, 1149, 1176, 1239, 1251, 1287, 1291, 1336, 1490, 1493, 1499, 1501, 1618), 
    43 : (0, 1, 8, 45, 54, 115, 190, 221, 231, 247, 302, 387, 426, 521, 535, 556, 632, 664, 790, 913, 1040, 1087, 1114, 1117, 1139, 1150, 1173, 1179, 1239, 1252, 1267, 1357, 1435, 1437, 1455, 1519, 1538, 1586, 1610, 1679, 1730, 1747, 1851, 1889), 
    47 : (0, 1, 39, 60, 115, 141, 146, 198, 268, 290, 303, 330, 337, 348, 378, 381, 398, 502, 504, 568, 763, 828, 837, 847, 853, 940, 989, 1013, 1056, 1109, 1181, 1237, 1260, 1314, 1408, 1423, 1587, 1601, 1629, 1633, 1718, 1730, 1941, 2004, 2012, 2041, 2123, 2159), 
    49 : (0, 1, 22, 26, 36, 42, 110, 342, 435, 484, 518, 535, 562, 639, 685, 740, 752, 842, 923, 988, 1003, 1012, 1042, 1231, 1289, 1307, 1380, 1418, 1513, 1544, 1676, 1709, 1762, 1770, 1867, 1878, 1907, 1930, 2055, 2092, 2137, 2140, 2142, 2199, 2206, 2278, 2321, 2334, 2353, 2381), 53 : (0, 1, 18, 90, 101, 354, 429, 490, 514, 612, 620, 622, 671, 731, 753, 797, 809, 849, 911, 1054, 1074, 1083, 1087, 1171, 1178, 1199, 1236, 1306, 1387, 1458, 1622, 1637, 1669, 1672, 1714, 1837, 1843, 1868, 1873, 1916, 1942, 1983, 2010, 2029, 2063, 2086, 2149, 2213, 2347, 2361, 2516, 2555, 2571, 2609), 
    59 : (0, 1, 28, 268, 363, 371, 681, 702, 761, 799, 848, 873, 980, 1025, 1030, 1059, 1149, 1160, 1184, 1214, 1253, 1260, 1262, 1275, 1345, 1426, 1489, 1508, 1566, 1810, 1863, 1866, 1883, 1899, 1930, 2016, 2028, 2071, 2137, 2155, 2350, 2354, 2425, 2497, 2541, 2593, 2603, 2635, 2793, 2844, 2850, 2949, 2972, 2986, 3109, 3169, 3279, 3367, 3393, 3501), 
    61 : (0, 1, 73, 159, 205, 343, 427, 507, 549, 568, 734, 791, 845, 876, 879, 884, 981, 1010, 1058, 1108, 1164, 1170, 1177, 1179, 1197, 1207, 1260, 1307, 1469, 1572, 1589, 1647, 1663, 1707, 1742, 1820, 1824, 1996, 2064, 2257, 2401, 2493, 2515, 2602, 2616, 2640, 2661, 2710, 2861, 2873, 3081, 3107, 3148, 3214, 3362, 3385, 3417, 3592, 3603, 3628, 3668, 3732), 
    64 : (0, 1, 80, 114, 116, 199, 274, 342, 346, 512, 662, 733, 798, 806, 955, 1138, 1273, 1299, 1377, 1438, 1580, 1623, 1689, 1710, 1726, 1815, 1827, 1903, 1933, 1995, 2015, 2024, 2123, 2390, 2501, 2660, 2673, 2712, 2835, 2868, 2931, 2937, 2975, 2990, 3021, 3031, 3085, 3274, 3360, 3441, 3489, 3492, 3511, 3534, 3539, 3566, 3744, 3768, 3779, 3793, 3913, 3920, 3980, 4104, 4144), 
    67 : (0, 1, 11, 45, 66, 96, 143, 222, 327, 416, 509, 524, 581, 631, 705, 722, 855, 860, 888, 1078, 1107, 1402, 1461, 1470, 1501, 1513, 1539, 1588, 1629, 1700, 1702, 1739, 1764, 1831, 1873, 2019, 2077, 2231, 2294, 2376, 2501, 2571, 2854, 2881, 2967, 2975, 3021, 3028, 3314, 3349, 3397, 3429, 3639, 3657, 3663, 3848, 3954, 3968, 4056, 4112, 4116, 4132, 4213, 4216, 4235, 4352, 4365, 4388), 
    71 : (0, 1, 204, 225, 345, 351, 392, 493, 609, 677, 700, 722, 762, 872, 931, 983, 1143, 1158, 1314, 1566, 1570, 1666, 1715, 1768, 1908, 2025, 2155, 2385, 2439, 2475, 2482, 2508, 2568, 2600, 2603, 2667, 2754, 2773, 2782, 2800, 2811, 2848, 3070, 3087, 3101, 3264, 3320, 3336, 3455, 3463, 3526, 3587, 3629, 3631, 3767, 3865, 3904, 3977, 3987, 4011, 4061, 4066, 4091, 4142, 4200, 4455, 4537, 4607, 4695, 4715, 5036, 5101), 
    73 : (0, 1, 48, 83, 334, 385, 484, 558, 591, 678, 796, 1043, 1046, 1087, 1143, 1165, 1331, 1440, 1469, 1483, 1489, 1667, 1703, 1943, 1995, 2019, 2058, 2244, 2301, 2303, 2378, 2399, 2409, 2424, 2436, 2490, 2540, 2563, 2657, 2689, 2817, 2857, 2885, 2970, 2987, 3032, 3116, 3127, 3146, 3180, 3188, 3206, 3357, 3580, 3835, 3851, 3996, 4107, 4199, 4208, 4212, 4279, 4416, 4421, 4486, 4611, 4618, 4721, 4884, 4953, 5136, 5222, 5260, 5315), 
    79 : (0, 1, 105, 116, 161, 170, 587, 655, 727, 981, 1118, 1135, 1213, 1229, 1335, 1401, 1405, 1419, 1478, 1480, 1642, 1664, 1766, 1853, 1878, 1910, 1945, 1995, 2008, 2189, 2220, 2280, 2317, 2398, 2405, 2413, 2451, 2544, 2572, 2578, 2608, 2791, 2993, 3019, 3048, 3249, 3301, 3352, 3375, 3527, 3665, 3677, 3864, 3971, 4011, 4244, 4380, 4385, 4494, 4533, 4543, 4576, 4623, 4699, 4881, 4939, 5161, 5232, 5280, 5381, 5400, 5532, 5730, 6001, 6004, 6025, 6045, 6087, 6114, 6222), 
    81 : (0, 1, 21, 88, 116, 196, 319, 322, 357, 382, 500, 519, 583, 636, 654, 667, 669, 907, 1175, 1191, 1276, 1285, 1390, 1465, 1647, 1655, 1819, 1876, 2064, 2068, 2108, 2115, 2227, 2281, 2320, 2326, 2331, 2505, 2657, 2796, 2813, 2885, 2943, 3040, 3089, 3162, 3185, 3231, 3287, 3369, 3537, 3783, 3844, 3892, 3914, 4147, 4174, 4184, 4287, 4361, 4489, 4714, 4874, 5025, 5208, 5260, 5284, 5352, 5384, 5418, 5627, 5717, 5776, 5947, 6117, 6198, 6228, 6511, 6523, 6537, 6566, 6602), 
    83 : (0, 1, 12, 87, 118, 199, 271, 295, 329, 372, 405, 442, 457, 573, 627, 724, 730, 747, 787, 1027, 1110, 1138, 1275, 1554, 1584, 1625, 1900, 1991, 2023, 2257, 2459, 2508, 2592, 2601, 2654, 2700, 2727, 2792, 2848, 2907, 3271, 3371, 3415, 3437, 3591, 3598, 3633, 3848, 4011, 4029, 4093, 4232, 4257, 4495, 4533, 4543, 4593, 4622, 4683, 4832, 5014, 5146, 5226, 5321, 5335, 5552, 5688, 5704, 5708, 5782, 5833, 5901, 5937, 5940, 6042, 6201, 6209, 6214, 6626, 6681, 6748, 6750, 6769, 6795), 
    89 : (0, 1, 32, 125, 160, 246, 734, 764, 1143, 1269, 1328, 1443, 1451, 1532, 1549, 1604, 1694, 1758, 1995, 2095, 2230, 2339, 2360, 2387, 2406, 2447, 2642, 2644, 2681, 2756, 2914, 2929, 2967, 3117, 3161, 3183, 3187, 3197, 3260, 3265, 3316, 3587, 3620, 3632, 3645, 3694, 3797, 4064, 4244, 4251, 4262, 4291, 4345, 4437, 4440, 4557, 4659, 4701, 5013, 5033, 5184, 5272, 5351, 5408, 5473, 5591, 5600, 5696, 5772, 5806, 6013, 6192, 6397, 6466, 6482, 6825, 6967, 6995, 7038, 7151, 7291, 7315, 7481, 7533, 7649, 7672, 7733, 7783, 7874, 7880), 
    97 : (0, 1, 41, 253, 365, 538, 631, 646, 702, 787, 809, 891, 958, 1242, 1318, 1382, 1401, 1424, 1463, 1710, 1757, 1895, 2205, 2216, 2374, 2733, 3104, 3195, 3320, 3423, 3687, 3804, 3821, 3824, 3870, 4138, 4150, 4239, 4372, 4382, 4518, 4547, 4615, 4673, 4717, 4724, 4726, 4801, 4921, 5021, 5131, 5323, 5344, 5382, 5442, 5474, 5731, 6006, 6063, 6288, 6410, 6462, 6475, 6489, 6505, 6604, 6610, 6628, 6665, 6673, 6698, 6900, 6935, 7173, 7317, 7471, 7543, 7700, 7748, 7774, 7861, 8013, 8041, 8137, 8255, 8328, 8333, 8495, 8683, 8719, 8769, 8773, 8988, 9019, 9135, 9169, 9249, 9376), 
    101 : (1, 101, 239, 324, 680, 688, 832, 848, 942, 1146, 1216, 1248, 1344, 1422, 1473, 1528, 1555, 1559, 1608, 1805, 1815, 1945, 1978, 2412, 2413, 2415, 2437, 2510, 2584, 2629, 2761, 2888, 2914, 2957, 2995, 3156, 3204, 3224, 3345, 3409, 3533, 3600, 3642, 3708, 4021, 4211, 4258, 4270, 4299, 4304, 4310, 4375, 4531, 5056, 5361, 5483, 5705, 5809, 5830, 6231, 6238, 6501, 6531, 6593, 6643, 6744, 6862, 6946, 7082, 7154, 7237, 7399, 7476, 7489, 7512, 7635, 7670, 7724, 7784, 7863, 7954, 8149, 8164, 8713, 8847, 8992, 9001, 9112, 9149, 9168, 9483, 9501, 9540, 9666, 9683, 9727, 9741, 9907, 10023, 10086, 10173, 10201)
    }
    
class DiffSetFactory:
    def __init__(self, q):
        [(p,e)] = factorize(q)
        self.p = p
        self.e = e
        self.q = q
        self.delta = q**2 + q + 1
        self.add = tuple((1, x) for x in range(self.delta))
        self.mul = tuple((u, 0) for u in range(self.delta) if gcd(u, self.delta) == 1)
        self.aff = tuple((u, x) for u, _ in self.mul for _, x in self.add)

# Code related to difference sets

    @staticmethod
    def possibly_ordered(ordered):
        def ordered_diff_set(x):
            return tuple(x)
        def unordered_diff_set(x):
            return tuple(sorted(x))
        if ordered:
            return ordered_diff_set
        else:
            return unordered_diff_set

    def transform_diff_set(self, diff_set, trafo, ordered = False):
        ordered_diff_set = self.possibly_ordered(ordered)
        u, x = trafo
        return ordered_diff_set((u * d + x) % self.delta for d in diff_set)
    
    def orbit_diff_set(self, diff_set, trafos = None, ordered = False):
        if trafos == None:
            trafos = self.aff
        return set(self.transform_diff_set(diff_set, trafo, ordered) for trafo in trafos)

    def aut_diff_set(self, diff_set, frobenius_only = False):
        diff_set = tuple(sorted(diff_set))
        if not frobenius_only:
            return set(trafo for trafo in self.aff
                       if self.transform_diff_set(diff_set, trafo) == diff_set)
        else:
            for _, x in self.add:
                if self.transform_diff_set(diff_set, (self.p, x)) == diff_set:
                    break
            else:
                raise ValueError('No frobenius found: %d, %s' % (self.q, diff_set))
            return tuple((self.p**i, sum(x*self.p**j for j in range(i))) for i in range(3*self.e))

    def minimize_diff_set(self, diff_set, trafos, ordered = False):
        outcome_trafo = defaultdict(lambda : [])
        for trafo in trafos:
            outcome_trafo[self.transform_diff_set(diff_set, trafo, ordered)].append(trafo)
#        print(diff_set, ordered, min(outcome_trafo), outcome_trafo[min(outcome_trafo)])
        return outcome_trafo[min(outcome_trafo)]

# Code related to difference matrices

    def difference_matrices(self, diff_set):
        return set(frozenset(tuple(row) for row in zip(diff_set, col1, col2))
                   for col1 in permutations(diff_set)
                   for col2 in permutations(diff_set))

    def transform_diff_mat(self, diff_mat, col_perm, trafo_triple):
        return frozenset(tuple((u * row[i] + x) % self.delta
                         for (u,x), i in zip(trafo_triple, col_perm))
                         for row in diff_mat)

    def diff_mat_diff_sets(self, diff_mat):
        return tuple(tuple(row[i] for row in diff_mat) for i in range(3))
    
    def orbit_diff_mat(self, diff_mat, col_perms = None, trafos = None):
        if col_perms == None:
            col_perms = permutations(range(3))
        if trafos == None:
            trafos = tuple(self.aut_diff_set(ds) for ds in self.diff_mat_diff_sets(diff_mat))
        return set(self.transform_diff_mat(diff_mat, perm, trafo_triple)
                   for perm in col_perms for trafo_triple in product(*trafos))

    @staticmethod
    def diff_mat_key(diff_mat):
        return tuple(sorted(x for row in diff_mat for x in row)) + tuple(sorted(diff_mat))

    @staticmethod
    def actual_matrix(diff_mat):
        return tuple(sorted(diff_mat))
    
    def based_diff_mat(self, diff_mat):
        diff_mat = self.actual_matrix(diff_mat)
        return tuple(tuple((e - b) % self.delta for e, b in zip(row, diff_mat[0])) for row in diff_mat)

    def based_diff_mats(self, diff_mat):
        diff_mat = self.actual_matrix(diff_mat)
        return set(frozenset(tuple((e - b) % self.delta for e, b in zip(row, brow)) for row in diff_mat) for brow in diff_mat)

    def minimize_diff_mat(self, diff_mat):
        diff_sets = [tuple(row[i] for row in diff_mat) for i in range(3)]
        minimize_unordered = [self.minimize_diff_set(diff_set, self.mul, ordered = False) for diff_set in diff_sets]
#        minimize_ordered = [self.minimize_diff_set(diff_set, trafos, ordered = True) for diff_set, trafos in zip(diff_sets, minimize_unordered)]
        minimized_diff_mat = min(self.actual_matrix(mat) for mat in self.orbit_diff_mat(diff_mat, col_perms = ((0, 1, 2),), trafos = minimize_unordered))
        return minimized_diff_mat
        
# Running as a script


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Produce list of difference matrices up to equivalence.')
    parser.add_argument('q', type = int, help = 'for which q to produce the difference matrices')
    parser.add_argument('--filename', default = 'diff_mats_%d.pickled', help = 'file to which to write the difference matrices, %%d will be replaced by q')
    parser.add_argument('--auts', action = 'store_true', help = 'whether to compute the automorphism groups of the difference matrices')
    parser.add_argument('--write_lattices', action = 'store_true', help = 'whether to (over)write the lattice file')
    parser.add_argument('--lattice_filename', default = 'lattices_%d.pickled', help = 'filename for the lattice file')
    args = parser.parse_args()
    dsf = DiffSetFactory(args.q)
    diff_set = min(tuple(sorted(ds)) for ds in dsf.orbit_diff_set(diff_sets[args.q]))
    mats = dsf.difference_matrices(diff_set)
    classes = eq_classes(mats, dsf.orbit_diff_mat)
    reps = tuple(sorted(dsf.actual_matrix(min(cls, key = dsf.diff_mat_key)) for cls in classes))
    breps = tuple(dsf.minimize_diff_mat(dsf.based_diff_mat(rep)) for rep in reps)
    try:
        filename = args.filename % args.q
    except TypeError:
        filename = args.filename
    pickl(filename, tuple(brep[1:] for brep in breps))
    res = unpickl(filename)
    print('%d matrices' % len(res))
    try:
        lattices = unpickl(args.lattice_filename % args.q)
    except IOError:
        lattices = { brep[1:] : dict() for brep in breps}
    for rep, brep in zip(reps, breps):
        print(brep)
        lattices[brep[1:]]['diff_mat'] = rep
        lattices[brep[1:]]['based_diff_mat'] = brep
        if args.auts:
            auts = aut_diff_mat(brep)
            lattices[brep[1:]]['aut_dm'] = auts
            print(auts)
    if args.write_lattices:
        pickl(args.lattice_filename % args.q, lattices)
