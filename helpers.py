from sympy import *
import math

def contact(dist, smoothing_dist):
    # "smoothed" step function
    return (1.+erf(dist/smoothing_dist))*0.5

def traction_force(kinetic_friction_coeff, contact_vel, normal_force, smoothing_vel):
    # TODO: stiction model
    magnitude = kinetic_friction_coeff*abs(normal_force)*erf(contact_vel.norm()/smoothing_vel)
    direction = -contact_vel / Piecewise((contact_vel.norm(), contact_vel.norm() > 0.), (1.,True))
    return magnitude*direction

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

def quat_312_roll(quat):
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

def quat_312_pitch(quat):
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
        if threshold == 0:
            ops_saved = 1
        else:
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

def tube_inertia_xx_yy(m, h, r1, r2):
    return 1./12. * m*(3*(r1**2+r2**2)+h**2)

def tube_inertia_zz(m, h, r1, r2):
    return 1./2. * m*(r1**2+r2**2)
