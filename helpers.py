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
