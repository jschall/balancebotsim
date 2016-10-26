from sympy import *

def quatderiv(quat, omega):
    P = omega[0]
    Q = omega[1]
    R = omega[2]

    return Matrix([[ 0, -P, -Q, -R],
                   [ P,  0,  R, -Q],
                   [ Q, -R,  0,  P],
                   [ R,  Q, -P,  0]]) * Matrix(quat)

def vector_in_frame(F, inlist):
    return F.x*inlist[0]+F.y*inlist[1]+F.z*inlist[2]

def count_subexpression(subexpr, expr):
    if hasattr(expr, "__getitem__"):
        return sum(map(lambda x: count_subexpression(subexpr, x), expr))
    else:
        return expr.count(subexpr)

def extractSubexpressions(inexprs, prefix='X', threshold=0, prev_subx=[]):
    subexprs, outexprs = cse(inexprs, symbols=numbered_symbols('__TMP__'), order='none')

    subexprs = prev_subx+subexprs

    for i in reversed(range(len(subexprs))):
        from sympy.logic.boolalg import Boolean
        ops_saved = (count_subexpression(subexprs[i][0], [[x[1] for x in subexprs], outexprs])-1)*subexprs[i][1].count_ops()
        if ops_saved < threshold or isinstance(subexprs[i][1], Boolean):
            sub = dict([subexprs.pop(i)])
            subexprs = [(x[0],x[1].xreplace(sub)) for x in subexprs]
            outexprs = [x.xreplace(sub) for x in outexprs]

    for i in range(len(subexprs)):
        newSym = Symbol('%s[%u]' % (prefix,i+len(prev_subx)))
        sub = {subexprs[i][0]:newSym}
        subexprs[i] = (newSym,subexprs[i][1])
        subexprs = [(x[0],x[1].xreplace(sub)) for x in subexprs]
        outexprs = [x.xreplace(sub) for x in outexprs]

    outexprs = [x.as_mutable() if type(x) is ImmutableDenseMatrix else x for x in outexprs]

    return tuple(outexprs+[subexprs])
