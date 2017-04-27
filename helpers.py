from sympy import *
from sympy.physics.mechanics import *
import math

def new_sym(prefix, sym_list, n=1, dynlevel=None):
    if isinstance(dynlevel,int):
        sym_list.extend(dynamicsymbols('%s%u:%u' % (prefix, 0, n), dynlevel))
    else:
        sym_list.extend(symbols('%s%u:%u' % (prefix, 0, n)))

    if n==1:
        return sym_list[-1]
    else:
        return sym_list[-n:]

def contact(dist, smoothing_dist):
    # "smoothed" step function
    return (1.+erf(dist/smoothing_dist))*0.5

def safe_normalize(v):
    ret = zeros(*v.shape)
    for i in range(len(v)):
        ret[i] = Piecewise((v[i]/v.norm(), v.norm() > 0.), (1. if i==0 else 0., True))

    return ret

def slope_intercept(p1,p2):
    assert len(p1)==2 and len(p2)==2
    return ((p1[1] - p2[1])/(p1[0] - p2[0]), (p1[0]*p2[1] - p1[1]*p2[0])/(p1[0] - p2[0]))

def coulomb_friction_model(v, normal_force, static_coeff, kinetic_coeff, smoothing_vel):
    v1, f1, v2, f2 = (smoothing_vel, static_coeff, 2*smoothing_vel, kinetic_coeff)
    m1, b1 = slope_intercept((0,0),(v1,f1))
    m2, b2 = slope_intercept((v1,f1),(v2,f2))

    return Piecewise((0, normal_force<0), (normal_force, True))*-Piecewise(
        (       -f2,              v < -v2),
        (-(m2*v+b2), And(v >= -v2, v < -v1)),
        (   m1*v+b1, And(v >= -v1, v <  v1)),
        (   m2*v+b2, And(v >=  v1, v <  v2)),
        (        f2, True))

def quat_multiply(q1, q2):
    return Matrix([q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3],
                   q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2],
                   q1[0]*q2[2] - q1[1]*q2[3] + q1[2]*q2[0] + q1[3]*q2[1],
                   q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1] + q1[3]*q2[0]])

def quatderiv(quat, omega):
    return 0.5*quat_multiply(quat,Matrix([0,omega[0],omega[1],omega[2]]))

def skew(v):
    return Matrix([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def quat_to_matrix(q):
    q = Matrix(q)
    return (q[0]**2-(q[1:,0].T*q[1:,0])[0])*eye(3) + 2.*(q[1:,0]*q[1:,0].T) + 2.*q[0]*skew(q[1:,0])

def quat_321_roll(quat):
    qr = quat[0]
    qi = quat[1]
    qj = quat[2]
    qk = quat[3]

    norm = sqrt(qr**2+qi**2+qj**2+qk**2)
    qr /= norm
    qi /= norm
    qj /= norm
    qk /= norm

    return math.atan2(2*(qr*qi+qj*qk), 1-2*(qi**2+qj**2))

def quat_321_pitch(quat):
    qr = quat[0]
    qi = quat[1]
    qj = quat[2]
    qk = quat[3]

    norm = sqrt(qr**2+qi**2+qj**2+qk**2)
    qr /= norm
    qi /= norm
    qj /= norm
    qk /= norm

    return math.asin(2*(qr*qj-qk*qi))

def vector_in_frame(F, inlist):
    return F.x*inlist[0]+F.y*inlist[1]+F.z*inlist[2]

def count_expression(count_expr, within_expr):
    if hasattr(within_expr, "__getitem__"):
        return sum(map(lambda x: count_expression(count_expr, x), within_expr))
    else:
        return within_expr.count(count_expr)

def extractSubexpressions(inexprs, prefix='X', threshold=1, prev_subx=[]):
    # Perform common subexpression extraction on inexprs and existing subexpressions
    subexprs, outexprs = cse(inexprs+[x[1] for x in prev_subx], symbols=numbered_symbols('__TMP__'), order='none')

    subexprs = list(zip([x[0] for x in prev_subx], outexprs[len(inexprs):]))+subexprs
    outexprs = outexprs[:len(inexprs)]

    done = False
    while not done:
        done = True

        for i in reversed(range(len(subexprs))):
            from sympy.logic.boolalg import Boolean
            ops_saved = (count_expression(subexprs[i][0], [[x[1] for x in subexprs], outexprs])-1)*subexprs[i][1].count_ops()
            if ops_saved < threshold or isinstance(subexprs[i][1], Boolean):
                sub = dict([subexprs.pop(i)])
                subexprs = [(x[0],x[1].xreplace(sub)) for x in subexprs]
                outexprs = [x.xreplace(sub) for x in outexprs]
                done = False
                break

    for i in range(len(subexprs)):
        newSym = Symbol('%s_%u' % (prefix,i))
        sub = {subexprs[i][0]:newSym}
        subexprs[i] = (newSym,subexprs[i][1])
        subexprs = [(x[0],x[1].xreplace(sub)) for x in subexprs]
        outexprs = [x.xreplace(sub) for x in outexprs]

    outexprs = [x.as_mutable() if type(x) is ImmutableDenseMatrix else x for x in outexprs]

    return tuple(outexprs+[subexprs])

def tube_inertia_xx_yy(m, h, r1, r2):
    return 1./12. * m*(3*(r1**2+r2**2)+h**2)

def tube_inertia_zz(m, h, r1, r2):
    return 1./2. * m*(r1**2+r2**2)
