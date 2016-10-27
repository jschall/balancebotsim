from helpers import *
from sympy import *
from sympy.physics.mechanics import *
from pydy.system import System
from pydy.codegen.ode_function_generators import LambdifyODEFunctionGenerator
import numpy.linalg
from datetime import datetime
import numpy as np

# Symbols
c = symbols('c0:20') # Constants
q = dynamicsymbols('q0:10')
qd = dynamicsymbols('q0:10',1)
u = dynamicsymbols('u0:9')
f = dynamicsymbols('f0:2')

# Define identifiers for constants
gravity = c[0]
pole_suspension_freq = c[1]
pole_suspension_zeta = c[2]
pole_length = c[3]
pole_mass = c[4]
cart_mass = c[5]
wheel_mass = c[6]
wheel_radius = c[7]
wheel_base = c[8]
ground_contact_freq = c[9]
ground_contact_zeta = c[10]

pole_inertia_xx = c[11]
pole_inertia_yy = c[12]
pole_inertia_zz = c[13]

wheel_inertia_xx = c[14]
wheel_inertia_yy = c[15]
wheel_inertia_zz = c[16]

cart_inertia_xx = c[17]
cart_inertia_yy = c[18]
cart_inertia_zz = c[19]

# Define identifiers for generalized coordinates
cart_quat = Matrix(q[0:4])
cart_pos = Matrix(q[4:7])
pole_theta = q[7]
lwheel_theta = q[8]
rwheel_theta = q[9]

# Define identifiers for generalized coordinate derivatives
cart_quat_dot = Matrix(qd[0:4])
cart_pos_dot = Matrix(qd[4:7])
pole_theta_dot = qd[7]
lwheel_theta_dot = qd[8]
rwheel_theta_dot = qd[9]

# Define identifiers for generalized speeds
cart_ang_vel = Matrix(u[0:3])
cart_vel = Matrix(u[3:6])
pole_omega = u[6]
lwheel_omega = u[7]
rwheel_omega = u[8]

# Define identifiers for force inputs
lwheel_motor_torque = f[0]
rwheel_motor_torque = f[1]
#lwheel_contact_force = Matrix(f[2:5])
#rwheel_contact_force = Matrix(f[5:8])

# Define newtonian reference frame and origin
N = ReferenceFrame('N')
O = Point('O')
O.set_vel(N, 0)

# Define lists of kinematic equations, forces and bodies that will be appended to
kde_list = []
force_list = []
body_list = []

# Define contact model constants
total_mass = pole_mass+cart_mass+wheel_mass*2
ground_contact_m = total_mass/2
ground_contact_k = ground_contact_m*(ground_contact_freq * 2.*pi)**2
ground_contact_c = 2.*ground_contact_zeta*sqrt(ground_contact_k*ground_contact_m)

# Define cart
cart_frame = ReferenceFrame('cart_frame')
cart_frame.orient(N, 'Quaternion', cart_quat)
cart_frame.set_ang_vel(N, vector_in_frame(N, cart_ang_vel))
cart_inertia = inertia(cart_frame, cart_inertia_xx, cart_inertia_yy, cart_inertia_zz)
cart_masscenter = O.locatenew('cart_masscenter', vector_in_frame(N, cart_pos))
cart_masscenter.set_vel(N, vector_in_frame(N, cart_vel))
cart_body = RigidBody('cart_body', cart_masscenter, cart_frame, cart_mass, (cart_inertia, cart_masscenter))

# Add cart to eqns
force_list.append((cart_masscenter, cart_mass*gravity*N.z))
kde_list.extend(list(cart_pos_dot-cart_vel)+list(cart_quat_dot-quatderiv(cart_quat, cart_ang_vel)))
body_list.append(cart_body)

# Define pole
pole_frame = ReferenceFrame('pole_frame')
pole_frame.orient(cart_frame, 'Axis', [pole_theta, cart_frame.x])
pole_inertia = inertia(pole_frame, pole_inertia_xx, pole_inertia_yy, pole_inertia_zz)
pole_masscenter = cart_masscenter.locatenew('pole_masscenter', -0.5*pole_length*pole_frame.z)
pole_masscenter.v2pt_theory(cart_masscenter, N, pole_frame)
pole_body = RigidBody('pole_body', pole_masscenter, pole_frame, pole_mass, (pole_inertia, pole_masscenter))
pole_top_point = pole_masscenter.locatenew('pole_top_point', -0.5*pole_length*pole_frame.z)
pole_top_point.v2pt_theory(pole_masscenter,N,pole_frame)
pole_top_pos = pole_top_point.pos_from(O).to_matrix(N)
pole_top_vel = pole_top_point.vel(N).to_matrix(N).subs(pole_theta_dot, pole_omega)
pole_top_ground_normal_force = 0.
pole_top_ground_normal_force = (-ground_contact_k*pole_top_pos[2] + -ground_contact_c*pole_top_vel[2])*((erf(10000*pole_top_pos[2])+1.)*0.5)

# Define suspension constants
pole_suspension_m = (pole_inertia+inertia_of_point_mass(pole_mass, pole_masscenter.pos_from(cart_masscenter), pole_frame)).dot(pole_frame.x).to_matrix(pole_frame)[0]
pole_suspension_k = pole_suspension_m*(pole_suspension_freq * 2.*pi)**2
pole_suspension_c = 2.*pole_suspension_zeta*sqrt(pole_suspension_k*pole_suspension_m)

# Add pole to eqns
force_list.append((pole_masscenter, pole_mass*gravity*N.z))
force_list.append((pole_frame, (-pole_suspension_k*pole_theta + -pole_suspension_c*pole_omega)*pole_frame.x))
force_list.append((cart_frame, -(-pole_suspension_k*pole_theta + -pole_suspension_c*pole_omega)*pole_frame.x))
force_list.append((pole_top_point, N.z*pole_top_ground_normal_force))
kde_list.append(pole_theta_dot-pole_omega)
body_list.append(pole_body)

# Define the ground direction vector
# This is done by subtracting the component of the N.z vector along cart_frame.y from N.z, resulting in only the component in the cart_frame x-z plane, and normalizing
ground_direction_vector = (N.z - N.z.dot(cart_frame.y)*cart_frame.y).normalize()

# Define lwheel
lwheel_frame = ReferenceFrame('lwheel_frame')
lwheel_frame.orient(cart_frame, 'Axis', [lwheel_theta, cart_frame.y])
lwheel_inertia = inertia(lwheel_frame, wheel_inertia_xx, wheel_inertia_yy, wheel_inertia_zz)
lwheel_masscenter = cart_masscenter.locatenew('lwheel_masscenter', -0.5*wheel_base*cart_frame.y)
lwheel_masscenter.v2pt_theory(cart_masscenter, N, lwheel_frame)
lwheel_body = RigidBody('lwheel_body', lwheel_masscenter, lwheel_frame, wheel_mass, (lwheel_inertia, lwheel_masscenter))
lwheel_contact_point = lwheel_masscenter.locatenew('lwheel_contact_point', ground_direction_vector*wheel_radius)
lwheel_contact_point.v2pt_theory(lwheel_masscenter,N,lwheel_frame)
lwheel_contact_pos = lwheel_contact_point.pos_from(O).to_matrix(N)
lwheel_contact_vel = lwheel_contact_point.vel(N).to_matrix(N)
lwheel_ground_normal_force = (-ground_contact_k*lwheel_contact_pos[2] + -ground_contact_c*lwheel_contact_vel[2])*((erf(10000*lwheel_contact_pos[2])+1.)*0.5)

# Add lwheel to eqns
force_list.append((lwheel_masscenter, wheel_mass*gravity*N.z))
force_list.append((lwheel_contact_point, N.z*lwheel_ground_normal_force))
force_list.append((lwheel_frame, lwheel_motor_torque*lwheel_frame.y))
force_list.append((cart_frame, -lwheel_motor_torque*lwheel_frame.y))
kde_list.append(lwheel_theta_dot-lwheel_omega)
body_list.append(lwheel_body)

# Define rwheel
rwheel_frame = ReferenceFrame('rwheel_frame')
rwheel_frame.orient(cart_frame, 'Axis', [rwheel_theta, cart_frame.y])
rwheel_inertia = inertia(rwheel_frame, wheel_inertia_xx, wheel_inertia_yy, wheel_inertia_zz)
rwheel_masscenter = cart_masscenter.locatenew('rwheel_masscenter', 0.5*wheel_base*cart_frame.y)
rwheel_masscenter.v2pt_theory(cart_masscenter, N, rwheel_frame)
rwheel_body = RigidBody('rwheel_body', rwheel_masscenter, rwheel_frame, wheel_mass, (rwheel_inertia, rwheel_masscenter))
rwheel_contact_point = rwheel_masscenter.locatenew('rwheel_contact_point', ground_direction_vector*wheel_radius)
rwheel_contact_point.v2pt_theory(rwheel_masscenter,N,rwheel_frame)
rwheel_contact_pos = rwheel_contact_point.pos_from(O).to_matrix(N)
rwheel_contact_vel = rwheel_contact_point.vel(N).to_matrix(N)
rwheel_ground_normal_force = (-ground_contact_k*rwheel_contact_pos[2] + -ground_contact_c*rwheel_contact_vel[2])*((erf(10000*rwheel_contact_pos[2])+1.)*0.5)

# Add rwheel to eqns
force_list.append((rwheel_masscenter, wheel_mass*gravity*N.z))
force_list.append((rwheel_contact_point, N.z*rwheel_ground_normal_force))
force_list.append((rwheel_frame, rwheel_motor_torque*rwheel_frame.y))
force_list.append((cart_frame, -rwheel_motor_torque*rwheel_frame.y))
kde_list.append(rwheel_theta_dot-rwheel_omega)
body_list.append(rwheel_body)

# Set constants
constants = {
        gravity: 9.80655,
        pole_suspension_freq: 2.,
        pole_suspension_zeta: 0.4,
        pole_length: 1.83,
        pole_mass: 0.39,
        cart_mass: 1.0,
        wheel_mass: 0.1,
        wheel_radius: .1,
        wheel_base: .2,
        ground_contact_freq: 10.,
        ground_contact_zeta: 0.25,
    }

constants.update({
        pole_inertia_xx: tube_inertia_xx_yy(pole_mass, pole_length, 0.03175, 0.03429).xreplace(constants),
        pole_inertia_yy: tube_inertia_xx_yy(pole_mass, pole_length, 0.03175, 0.03429).xreplace(constants),
        pole_inertia_zz: tube_inertia_zz(pole_mass, pole_length, 0.03175, 0.03429).xreplace(constants),

        wheel_inertia_xx: tube_inertia_xx_yy(wheel_mass, .02, wheel_radius*0.9, wheel_radius).xreplace(constants),
        wheel_inertia_yy: tube_inertia_zz(wheel_mass, .02, wheel_radius*0.9, wheel_radius).xreplace(constants),
        wheel_inertia_zz: tube_inertia_xx_yy(wheel_mass, .02, wheel_radius*0.9, wheel_radius).xreplace(constants),

        cart_inertia_xx: tube_inertia_xx_yy(cart_mass, wheel_base, wheel_radius*0.3, wheel_radius*0.9).xreplace(constants),
        cart_inertia_yy: tube_inertia_zz(cart_mass, wheel_base, wheel_radius*0.3, wheel_radius*0.9).xreplace(constants),
        cart_inertia_zz: tube_inertia_xx_yy(cart_mass, wheel_base, wheel_radius*0.3, wheel_radius*0.9).xreplace(constants)
        })

KM = KanesMethod(N, q_ind=q, u_ind=u, kd_eqs=kde_list)
KM.kanes_equations(force_list, body_list)

s = System(KM, times=np.linspace(0,20,20000))

s.constants=constants

r = 0.
p = 0.
y = 0.
s.initial_conditions = {
        cart_quat[1]: sin(r/2)*cos(p/2)*cos(y/2) - cos(r/2)*sin(p/2)*sin(y/2),
        cart_quat[2]: cos(r/2)*sin(p/2)*cos(y/2) + sin(r/2)*cos(p/2)*sin(y/2),
        cart_quat[3]: cos(r/2)*cos(p/2)*sin(y/2) - sin(r/2)*sin(p/2)*cos(y/2),
        cart_quat[0]: cos(r/2)*cos(p/2)*cos(y/2) + sin(r/2)*sin(p/2)*sin(y/2),

        pole_theta: 0.0,
        cart_pos[2]: -0.05-wheel_radius.xreplace(constants),
        lwheel_omega: 0.
    }

s.specifieds = {
        lwheel_motor_torque: 0.,
        rwheel_motor_torque: 0.,
        #lwheel_contact_force[0]: 0.,
        #lwheel_contact_force[1]: 0.,
        #lwheel_contact_force[2]: 0.,
        #rwheel_contact_force[0]: 0.,
        #rwheel_contact_force[1]: 0.,
        #rwheel_contact_force[2]: 0.
    }



s.generate_ode_function(generator='cython')

print(1)
result = s.integrate()
print(1)

import matplotlib.pyplot as plt
t = s.times
plt.subplot(211)
plt.plot(t,[quat_312_roll(x[0:4]) for x in result], t,[quat_312_pitch(x[0:4]) for x in result])
plt.subplot(212)
plt.plot(t,[x[6] for x in result])
plt.show()
