import sympy
import re
sympy.init_printing(use_latex=True)

def matrix_to_rational(matrix_a):
    return sympy.Matrix([[sympy.nsimplify(elem) for elem in row] for row in matrix_a])

# Regular expression for row operations. Matches the following three kinds of expression:
#
# R<int> <=> R<int>
# R<int> * (<factor>) => R<int>
# R<int> * (<factor>) + R<int> => R<int>
#
# where:
# <int> is a decimal integer
# <factor> is either a <decimal> or <decimal>/<decimal>, where <decimal> is a decimal number.
pattern = r'''^
 \s* R(?P<source_row> \d+)                          # Source row
 \s* (
      ( <=> \s* R(?P<swap_row> \d+))                # Swap row
      | 
      ( 
        \*                                          # Literal * 
        \s* \(                                      # Literal (
        \s* (?P<num>[+-]?((\d+(\.\d*)?)|(\.\d+)))   # Numerator
        (                                           # Optional denominator
          \s* /                                     # Literal /
          \s* (?P<den>[+-]?((\d+(\.\d*)?)|(\.\d+))) # Denominator
        )?                                          
        \s* \)                                      # Literal )
        (                                           # Optional add row
          \s* \+                                    # Literal +
          \s* R(?P<add_row> \d+)                    # Add row
        )?
        \s* => \s* R(?P<target_row> \d+)            # Target row
      )
     )
\s* $
'''
re_rop = re.compile(pattern, re.I | re.X)

t_swap = 0
t_scale = 1
t_scale_add = 2

def rop_compile(rop_str, symbolic=True):
    def error_msg(error_string):
        return (error_string + ': {}').format(rop_str)
    m = re_rop.match(rop_str)
    if m is None:
        raise ValueError(error_msg('invalid row operation'))
    md = m.groupdict()
    nsource = int(md['source_row']) - 1
    if md['swap_row'] is not None:
        nswap = int(md['swap_row']) - 1
        return (t_swap, nsource, nswap, 0)
    nnum = eval(md['num'])
    nden = 1 if md['den'] is None else eval(md['den'])
    if nden == 0:
        raise ValueError(error_msg('zero denominator in row operation'))
    c = sympy.Rational(nnum, nden) if symbolic else nnum / nden
    ntarget = int(md['target_row']) - 1
    if md['add_row'] is None:
        if nsource != ntarget:
            raise ValueError(error_msg('source and target rows are not the same'))
        if c == 0:
            raise ValueError(error_msg('zero scaling factor in row operation'))
        return (t_scale, nsource, ntarget, c)
    nadd = int(md['add_row']) - 1
    if nadd != ntarget:
        raise ValueError(error_msg('add and target rows are not the same'))
    if nadd == nsource:
        raise ValueError(error_msg('add and source rows are the same'))
    return (t_scale_add, nsource, nadd, c)

def rop_swap(matrix_a, i, j, inplace=False):
    matrix_b = matrix_a if inplace else matrix_a[:, :]
    matrix_b[i, :], matrix_b[j, :] = matrix_b[j, :], matrix_b[i, :]
    return matrix_b


def rop_scale(matrix_a, i, c, inplace=False):
    matrix_b = matrix_a if inplace else matrix_a[:, :]
    matrix_b[i, :] = c * matrix_b[i, :]
    return matrix_b


def rop_scale_add(matrix_a, i, c, j, inplace=False):
    matrix_b = matrix_a if inplace else matrix_a[:, :]
    matrix_b[j, :] += c * matrix_b[i, :]
    return matrix_b

def do_rop(matrix_a, t, r1, r2, c, inplace=False):
    case_dict = {
        t_swap: lambda a, r1, r2, c, inplace: rop_swap(a, r1, r2, inplace),
        t_scale: lambda a, r1, r2, c, inplace: rop_scale(a, r1, c, inplace),
        t_scale_add: lambda a, r1, r2, c, inplace: rop_scale_add(a, r1, c, r2, inplace)
    }
    return case_dict[t](matrix_a, r1, r2, c, inplace)

def rop(matrix_a, *rop_seq, symbolic=True, inplace=False):
    matrix_b = matrix_a if inplace else matrix_a[:, :]
    for rop_str in rop_seq:
        t, r1, r2, c = rop_compile(rop_str, symbolic)
        do_rop(matrix_b, t, r1, r2, c, inplace=True)
    return matrix_b