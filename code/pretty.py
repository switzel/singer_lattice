import collections
import re

class Str:
    @staticmethod
    def str(content, latex = False):
        if content == None:
            return ''
        else:
            return str(content)

    default = None

class List:
    def __init__(self, types, sep = ' \t ', latex_sep = ' ', l = 3, raise_long = None):
        def list_str(lst, latex = False):
            if not latex:
                return sep.join(types.str(x, latex = False) for x in lst)
            else:
                if self.raise_long != None:
                    return latex_sep.join('\\raisebox{%.1fem}{%s}' % (((i + raise_long) % 2)-.5, types.str(x, latex = True)) for i, x in enumerate(lst))
                else:
                    return latex_sep.join(types.str(x, latex = True) for x in lst)
        self.raise_long = raise_long
        self.str = list_str
        self.default = [types.default for i in range(l)]

class Matrix(Str):
    @staticmethod
    def str(content, latex = False):
        if content and latex:
            return '\\left(\\begin{smallarray}{ccc}'+\
                   '\\\\'.join('&'.join(str(e) for e in row) for row in content) +\
                   '\\end{smallarray}\\right)'
        else:
            return Str.str(content, latex = latex)

class Homology:
    @staticmethod
    def str(content, latex = False):
        if not content:
            return '0'
        if latex:
            freq = collections.Counter(content)
            return r'\oplus'.join(r'(\Z/%d\Z)^{%d}' % (k,e) if e > 1 else r'\Z/%d\Z' % k for (k, e) in freq.items())
        else:
            return ' '.join(str(e) for e in content)

    default = tuple()

class Space(Str):
    @staticmethod
    def str(content, latex = False):
        if latex:
            return '\makebox[3em][c]{$%s$}' % (content, )
        else:
            return Str.str(content, latex = latex)

    def __init__(self, typ = Str, width = 3):
        def this_str(content, latex = False):
            if latex:
                return '\makebox[%.1fem][c]{$%s$}' % (width, typ.str(content, latex = latex))
            else:
                return typ.str(content, latex = latex)
        self.str = this_str
            
class Maybe(Str):
    @staticmethod
    def str(content, latex = False):
        return { True: '+', False: '-', None: ' ' }[content]

class AutDiffMat(Str):
    @staticmethod
    def str(content, latex = False):
        if content == None:
            return ''
        aff, symm = content
        gp_strs = []
        if len(aff) > 1:
            gp_strs.append('C_%d' % len(aff))
        if len(symm) == 1:
            pass
        elif len(symm) == 2:
            gp_strs.append(str(symm[1]))
        elif len(symm) == 3:
            gp_strs.append('C_3')
        elif len(symm) == 6:
            gp_strs.append('S_3')
        else:
            raise ValueError, 'Invalid permutation group %s' % (symm,)
        if not gp_strs:
            gp_strs.append('1')
        if latex:
            return ' \\rtimes '.join(gp_strs)
        else:
            return ' ><| '.join(gp_strs)

class Order(Str):
    @staticmethod
    def str(content, latex = False):
        if content == None or content >= 0:
            return Str.str(content, latex = latex)
        else:
            return '%s %s' % ('\ge' if latex else '>=',\
                              Str.str(-content, latex = latex))

##def matrix_str(mat, latex = False):
##    if mat == None:
##        return ''
##    if latex:
##        result = '\\left(\\begin{smallarray}{ccc}'
##        result += '\\\\'.join('&'.join(str(e) for e in row) for row in mat)
##        result += '\\end{smallarray}\\right)'
##        return result
##    else:
##        return str(mat)
##
##def homology_str(lst, latex = False):
##    if lst == None:
##        return ''
##    if latex:
##        if list(lst) != []:
##            freq = collections.Counter(lst)
##            return r'\oplus'.join(r'(\Z/%d\Z)^{%d}' % (k,e) if e > 1 else r'\Z/%d\Z' % k for (k, e) in freq.items())
##        else:
##            return '0'
##    else:
##        return ' '.join(str(e) for e in lst)
##
##def space_str(a, latex = False):
##    if a == None:
##        a = ''
##    if latex:
##        return '\makebox[3em][c]{$%s$}' % (a,)
##    else:
##        return str(a)
##
##def plain_str(a, latex = False):
##    if a == None:
##        return ''
##    return str(a)
##
##def maybe_str(a, latex = False):
##    if latex:
##        trans = { True : '\makebox[1em][c]{+}',
##                  False : '\makebox[1em][c]{-}',
##                  None : '\makebox[1em][c]{?}'}
##    else:
##        trans = { True : 'True', False : 'False', None : 'Maybe' }
##    try:
##        return trans[a]
##    except KeyError:
##        raise ValueError("Only, True, False, None (=Maybe) allowed, got %s" % a)
##
##def list_str(print_fcn):
##    def result(lst, latex = False):
##        if latex:
##            return ' '.join(print_fcn(x, latex = True) for x in lst)
##        else:
##            return '\t '.join(print_fcn(x, latex = False) for x in lst)
##    return result
##
##def aut_diff_mat_str(aut_dm, latex = False):
##    if aut_dm == None:
##        return ''
##    aff, symm = aut_dm
##    gp_strs = []
##    if len(aff) > 1:
##        gp_strs.append('C_%d' % len(aff))
##    if len(symm) == 1:
##        pass
##    elif len(symm) == 2:
##        gp_strs.append(str(symm[1]))
##    elif len(symm) == 3:
##        gp_strs.append('C_3')
##    elif len(symm) == 6:
##        gp_strs.append('S_3')
##    else:
##        raise ValueError, 'Invalid permutation group %s' % (symm,)
##    if not gp_strs:
##        gp_strs.append('1')
##    if latex:
##        return ' \\rtimes '.join(gp_strs)
##    else:
##        return ' ><| '.join(gp_strs)
##
##def stab_str(stab, latex = False):
##    if stab == None:
##        if latex:
##            return '\makebox[1.7em][c]{}'
##        else:
##            return ''
##    if stab < 0:
##        if latex:
##            return '\makebox[1.7em][c]{$\\ge %d$}' % (-stab)
##        else:
##            return '>= %d' % (-stab)
##    else:
##        if latex:
##            return '\makebox[1.7em][c]{$%d$}' % stab
##        else:
##            return '%d' % stab
##
####    match = re.match('K_(\d+)^(\d+)', key)
####    if not match:
####        raise ValueError, 'Unknown stabilizer: %s' % key
####    n = int(match.group(1))
####    k = int(match.group(2))
####    if n != 1:
####        raise NotImplemented, 'Do not know how to print stabilizers of larger balls: %s' %key
####    if stab == 1:
####        result = ''
####    else:
####        if latex:
####            result = 'C_{%d}' % stab
####        else:
####            result = 'C_%d' % stab
####    if k != 2:
####        result += 'P_{%d}^{%d}' % (n, k-1)
####    if result == '':
####        result = '1'
####    return result
