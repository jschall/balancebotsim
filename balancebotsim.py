from helpers import *
from sympy import *
from sympy.physics.mechanics import *
from pydy.system import System
from pydy.codegen.ode_function_generators import LambdifyODEFunctionGenerator
from scipy.integrate import odeint
import numpy.linalg
from datetime import datetime
import numpy as np
import sys

from balancebot_constants import *

# Symbols
q_sym = [] # Generalized coordinates
qd_sym = [] # Generalized coordinate derivatives
u_sym = [] # Generalized speeds
in_sym = [] # Inputs

def new_q(name, n=1):
    return (new_sym(name, q_sym, n=n, dynlevel=0), new_sym(name, qd_sym, n=n, dynlevel=1))

def new_u(name, n=1):
    return new_sym(name, u_sym, n, dynlevel=0)

def new_in(name, n=1):
    return new_sym(name, in_sym, n, dynlevel=0)

# Define identifiers for generalized coordinates and derivatives
cart_quat, cart_quat_dot = new_q('quat', 4)
cart_pos, cart_pos_dot = new_q('pos', 3)
pole_theta, pole_theta_dot = new_q('pole_theta')
lwheel_theta, lwheel_theta_dot = new_q('lwheel_theta')
rwheel_theta, rwheel_theta_dot = new_q('rwheel_theta')

# Define identifiers for generalized speeds
cart_ang_vel = new_u('ang_vel', 3)
cart_vel = new_u('vel',3)
pole_omega = new_u('pole_omega')
lwheel_omega = new_u('lwheel_omega')
rwheel_omega = new_u('rwheel_omega')

# Define identifiers for force inputs
right_motor_torque = new_in('right_motor_torque')
left_motor_torque = new_in('left_motor_torque')

# Define newtonian reference frame and origin, gravity
N = ReferenceFrame('N')
O = Point('O')
O.set_vel(N, 0)
gravity = g*N.z

# Define lists of kinematic equations, forces and bodies that will be appended to
kdes = []
forces = []
bodies = []

# Define cart
cart_frame = ReferenceFrame('cart_frame')
cart_frame.orient(N, 'Quaternion', cart_quat)
cart_frame.set_ang_vel(N, vector_in_frame(cart_frame, cart_ang_vel))
cart_inertia = inertia(cart_frame, cart_inertia_xx, cart_inertia_yy, cart_inertia_zz)
cart_masscenter = O.locatenew('cart_masscenter', vector_in_frame(N, cart_pos))
cart_masscenter.set_vel(N, vector_in_frame(N, cart_vel))
cart_body = RigidBody('cart', cart_masscenter, cart_frame, cart_mass, (cart_inertia, cart_masscenter))

# Add cart to eqns
forces.append((cart_masscenter, cart_mass*gravity))
kdes.append(cart_pos_dot[0]-cart_vel[0])
kdes.append(cart_pos_dot[1]-cart_vel[1])
kdes.append(cart_pos_dot[2]-cart_vel[2])
kdes.extend(kinematic_equations(cart_ang_vel, cart_quat, 'quaternion'))
bodies.append(cart_body)

# Define pole
pole_frame = ReferenceFrame('pole_frame')
pole_frame.orient(cart_frame, 'Axis', [pole_theta, cart_frame.x])
pole_inertia = inertia(pole_frame, pole_inertia_xx_yy, pole_inertia_xx_yy, pole_inertia_zz)
pole_masscenter = cart_masscenter.locatenew('pole_masscenter', -0.5*pole_length*pole_frame.z)
pole_masscenter.v2pt_theory(cart_masscenter, N, pole_frame)
pole_body = RigidBody('pole', pole_masscenter, pole_frame, pole_mass, (pole_inertia, pole_masscenter))
pole_top_point = pole_masscenter.locatenew('pole_top_point', -0.5*pole_length*pole_frame.z)
pole_top_point.v2pt_theory(pole_masscenter,N,pole_frame)
pole_top_pos = pole_top_point.pos_from(O).to_matrix(N)
pole_top_vel = pole_top_point.vel(N).to_matrix(N).subs(pole_theta_dot, pole_omega) # TODO: Why is this subs necessary?
pole_top_ground_normal_force = (-ground_contact_k*pole_top_pos[2] + -ground_contact_c*pole_top_vel[2])*contact(pole_top_pos[2], contact_smoothing_dist)

# Define battery
battery_frame = pole_frame
battery_inertia = inertia(battery_frame, battery_inertia_xx_yy, battery_inertia_xx_yy, battery_inertia_zz)
battery_masscenter = cart_masscenter.locatenew('battery_masscenter', -0.5*battery_length*battery_frame.z)
battery_masscenter.v2pt_theory(cart_masscenter,N,battery_frame)
battery_body = RigidBody('battery', battery_masscenter, battery_frame, battery_mass, (battery_inertia, battery_masscenter))

# Define payload
payload_frame = pole_frame
payload_inertia = inertia(payload_frame, payload_inertia_xx, payload_inertia_yy, payload_inertia_zz)
payload_masscenter = cart_masscenter.locatenew('payload_masscenter', -payload_position_h*payload_frame.z + payload_position_x*payload_frame.x)
payload_masscenter.v2pt_theory(cart_masscenter,N,battery_frame)
payload_body = RigidBody('payload', payload_masscenter, payload_frame, payload_mass, (payload_inertia, payload_masscenter))

# Define suspension constants
pole_suspension_m = (pole_inertia+battery_inertia+payload_inertia+inertia_of_point_mass(pole_mass, pole_masscenter.pos_from(cart_masscenter), pole_frame)+inertia_of_point_mass(battery_mass, battery_masscenter.pos_from(cart_masscenter), battery_frame)+inertia_of_point_mass(payload_mass, payload_masscenter.pos_from(cart_masscenter), payload_frame)).dot(pole_frame.x).to_matrix(pole_frame)[0]
pole_suspension_k = pole_suspension_m*(pole_suspension_freq * 2.*pi)**2
pole_suspension_c = 2.*pole_suspension_zeta*sqrt(pole_suspension_k*pole_suspension_m)

# Add pole to eqns
forces.append((pole_masscenter, pole_mass*gravity))
forces.append((pole_frame, (-pole_suspension_k*pole_theta + -pole_suspension_c*pole_omega)*pole_frame.x))
forces.append((cart_frame, -(-pole_suspension_k*pole_theta + -pole_suspension_c*pole_omega)*pole_frame.x))
forces.append((pole_top_point, N.z*pole_top_ground_normal_force))
kdes.append(pole_theta_dot-pole_omega)
bodies.append(pole_body)

# Add battery to eqns
forces.append((battery_masscenter, battery_mass*gravity))
bodies.append(battery_body)

# Add payload to eqns
forces.append((payload_masscenter, payload_mass*gravity))
bodies.append(payload_body)

# Define the ground direction vector and a wall direction vector
# This is done by subtracting the component of the N.z vector along cart_frame.y from N.z, resulting in only the component in the cart_frame x-z plane, and normalizing
ground_direction_vector = (N.z - N.z.dot(cart_frame.y)*cart_frame.y).normalize()
#wall_direction_vector = (N.x - N.x.dot(cart_frame.y)*cart_frame.y).normalize()

# Define lwheel
lwheel_frame = ReferenceFrame('lwheel_frame')
lwheel_frame.orient(cart_frame, 'Axis', [lwheel_theta, cart_frame.y])
lwheel_inertia = inertia(lwheel_frame, wheel_inertia_xx_zz, wheel_inertia_yy, wheel_inertia_xx_zz)
lwheel_masscenter = cart_masscenter.locatenew('lwheel_masscenter', -0.5*wheel_base*cart_frame.y)
lwheel_masscenter.v2pt_theory(cart_masscenter, N, lwheel_frame)
lwheel_body = RigidBody('lwheel', lwheel_masscenter, lwheel_frame, wheel_mass, (lwheel_inertia, lwheel_masscenter))

# Set up lwheel contact model
lwheel_ground_contact_point = lwheel_masscenter.locatenew('lwheel_ground_contact_point', ground_direction_vector*wheel_radius)
lwheel_ground_contact_point.v2pt_theory(lwheel_masscenter,N,lwheel_frame)
lwheel_ground_contact_pos = lwheel_ground_contact_point.pos_from(O).to_matrix(N)
lwheel_ground_contact_vel = lwheel_ground_contact_point.vel(N).to_matrix(N).subs(lwheel_theta_dot, lwheel_omega) # TODO: Why is this subs necessary?
lwheel_ground_contact_force = zeros(3,1)
lwheel_ground_contact_force[2] = (-ground_contact_k*lwheel_ground_contact_pos[2] + -ground_contact_c*lwheel_ground_contact_vel[2])*contact(lwheel_ground_contact_pos[2], contact_smoothing_dist)
lwheel_ground_contact_force[0:2,0] = coulomb_friction_model(lwheel_ground_contact_vel[0:2,0].norm(), -lwheel_ground_contact_force[2], wheel_ground_mu_s, wheel_ground_mu_k, friction_smoothing_vel)*safe_normalize(lwheel_ground_contact_vel[0:2,0])

#lwheel_wall_contact_point = lwheel_masscenter.locatenew('lwheel_wall_contact_point', wall_direction_vector*wheel_radius)
#lwheel_wall_contact_point.v2pt_theory(lwheel_masscenter,N,lwheel_frame)
#lwheel_wall_contact_pos = lwheel_wall_contact_point.pos_from(O).to_matrix(N)
#lwheel_wall_contact_vel = lwheel_wall_contact_point.vel(N).to_matrix(N).subs(lwheel_theta_dot, lwheel_omega) # TODO: Why is this subs necessary?
#lwheel_wall_contact_force = zeros(3,1)
#lwheel_wall_contact_force[0] = (-ground_contact_k*lwheel_wall_contact_pos[0] + -ground_contact_c*lwheel_wall_contact_vel[0])*contact(lwheel_wall_contact_pos[0], contact_smoothing_dist)
#lwheel_wall_contact_force[1:3,0] = coulomb_friction_model(lwheel_wall_contact_vel[1:3,0].norm(), -lwheel_wall_contact_force[0], wheel_ground_mu_s, wheel_ground_mu_k, friction_smoothing_vel)*safe_normalize(lwheel_wall_contact_vel[1:3,0])

# Add lwheel to eqns
forces.append((lwheel_masscenter, wheel_mass*gravity))
forces.append((lwheel_ground_contact_point, vector_in_frame(N, lwheel_ground_contact_force)))
#forces.append((lwheel_wall_contact_point, vector_in_frame(N, lwheel_wall_contact_force)))
forces.append((lwheel_frame, left_motor_torque*lwheel_frame.y))
forces.append((cart_frame, -left_motor_torque*lwheel_frame.y))
kdes.append(lwheel_theta_dot-lwheel_omega)
bodies.append(lwheel_body)

# Define rwheel
rwheel_frame = ReferenceFrame('rwheel_frame')
rwheel_frame.orient(cart_frame, 'Axis', [rwheel_theta, cart_frame.y])
rwheel_inertia = inertia(rwheel_frame, wheel_inertia_xx_zz, wheel_inertia_yy, wheel_inertia_xx_zz)
rwheel_masscenter = cart_masscenter.locatenew('rwheel_masscenter', 0.5*wheel_base*cart_frame.y)
rwheel_masscenter.v2pt_theory(cart_masscenter, N, rwheel_frame)
rwheel_body = RigidBody('rwheel', rwheel_masscenter, rwheel_frame, wheel_mass, (rwheel_inertia, rwheel_masscenter))

# Set up rwheel contact model
rwheel_ground_contact_point = rwheel_masscenter.locatenew('rwheel_ground_contact_point', ground_direction_vector*wheel_radius)
rwheel_ground_contact_point.v2pt_theory(rwheel_masscenter,N,rwheel_frame)
rwheel_ground_contact_pos = rwheel_ground_contact_point.pos_from(O).to_matrix(N)
rwheel_ground_contact_vel = rwheel_ground_contact_point.vel(N).to_matrix(N).subs(rwheel_theta_dot, rwheel_omega) # TODO: Why is this subs necessary?
rwheel_ground_contact_force = zeros(3,1)
rwheel_ground_contact_force[2] = (-ground_contact_k*rwheel_ground_contact_pos[2] + -ground_contact_c*rwheel_ground_contact_vel[2])*contact(rwheel_ground_contact_pos[2], contact_smoothing_dist)
rwheel_ground_contact_force[0:2,0] = coulomb_friction_model(rwheel_ground_contact_vel[0:2,0].norm(), -rwheel_ground_contact_force[2], wheel_ground_mu_s, wheel_ground_mu_k, friction_smoothing_vel)*safe_normalize(rwheel_ground_contact_vel[0:2,0])

#rwheel_wall_contact_point = rwheel_masscenter.locatenew('rwheel_wall_contact_point', wall_direction_vector*wheel_radius)
#rwheel_wall_contact_point.v2pt_theory(rwheel_masscenter,N,rwheel_frame)
#rwheel_wall_contact_pos = rwheel_wall_contact_point.pos_from(O).to_matrix(N)
#rwheel_wall_contact_vel = rwheel_wall_contact_point.vel(N).to_matrix(N).subs(rwheel_theta_dot, rwheel_omega) # TODO: Why is this subs necessary?
#rwheel_wall_contact_force = zeros(3,1)
#rwheel_wall_contact_force[0] = (-ground_contact_k*rwheel_wall_contact_pos[0] + -ground_contact_c*rwheel_wall_contact_vel[0])*contact(rwheel_wall_contact_pos[0], contact_smoothing_dist)
#rwheel_wall_contact_force[1:3,0] = coulomb_friction_model(rwheel_wall_contact_vel[1:3,0].norm(), -rwheel_wall_contact_force[0], wheel_ground_mu_s, wheel_ground_mu_k, friction_smoothing_vel)*safe_normalize(rwheel_wall_contact_vel[1:3,0])

# Add rwheel to eqns
forces.append((rwheel_masscenter, wheel_mass*gravity))
forces.append((rwheel_ground_contact_point, vector_in_frame(N, rwheel_ground_contact_force)))
#forces.append((rwheel_wall_contact_point, vector_in_frame(N, rwheel_wall_contact_force)))
forces.append((rwheel_frame, right_motor_torque*rwheel_frame.y))
forces.append((cart_frame, -right_motor_torque*rwheel_frame.y))
kdes.append(rwheel_theta_dot-rwheel_omega)
bodies.append(rwheel_body)

KM = KanesMethod(N, q_ind=q_sym, u_ind=u_sym, kd_eqs=kdes)
KM.kanes_equations(forces, bodies)

kdd = KM.kindiffdict()
mm = KM.mass_matrix_full
fo = KM.forcing_full

subs = {}
subs.update(kdd)
subs.update(dict(zip(q_sym+u_sym, symbols('x[0:%u]' % (len(q_sym+u_sym),)))))
subs.update(dict(zip(in_sym, symbols('in[0:%u]' % (len(in_sym),)))))

mm = mm.xreplace(subs)
fo = fo.xreplace(subs)

#mm, fo, subx1 = extractSubexpressions([mm, fo], 'subx1')

eom = mm.LUsolve(fo)
eom, subx = extractSubexpressions([eom], 'subx')

print(count_ops(eom)+count_ops(subx))

output = {'eom': {'subx': subx, 'param': Matrix(q_sym+u_sym+in_sym), 'ret': eom}}
for b in bodies:
    output[str(b)+'_masscenter_pos'] = {'param': Matrix(q_sym+u_sym), 'ret': simplify(b.masscenter.pos_from(O).to_matrix(N))}
    output[str(b)+'_masscenter_vel'] = {'param': Matrix(q_sym+u_sym), 'ret': simplify(b.masscenter.vel(N).to_matrix(N))}
    output[str(b)+'_rot_body_to_newtonian'] = {'param': Matrix(q_sym+u_sym), 'ret': simplify(N.dcm(b.frame))}

with open('dyn.h', 'w') as f:
    f.truncate()
    f.write('#pragma once\n')
    f.write('#define N_STATES %u\n' % (len(q_sym+u_sym),))
    f.write('void dyn(const float* x, const float* in, float* ret) {\n')

    for name,expr in output['eom']['subx']:
        f.write('    float %s = %s;\n' % (name,ccode(expr)))

    f.write('\n')

    for i in range(len(output['eom']['ret'])):
        f.write('ret[%u] = %s\n' % (i,ccode(output['eom']['ret'][i])))

    f.write('}\n')
