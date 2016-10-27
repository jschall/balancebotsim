from helpers import *
from sympy import *
from sympy.physics.mechanics import *
from pydy.system import System
from pydy.codegen.ode_function_generators import LambdifyODEFunctionGenerator
from scipy.integrate import odeint
import numpy.linalg
from datetime import datetime
import numpy as np

# Symbols
p_sym = symbols('c0:23') # Constants
q_sym = dynamicsymbols('q0:10')
qd_sym = dynamicsymbols('q0:10',1)
u_sym = dynamicsymbols('u0:9')
in_sym = dynamicsymbols('in0:2')

# Define identifiers for constants
g = p_sym[0]
pole_suspension_freq = p_sym[1]
pole_suspension_zeta = p_sym[2]
pole_length = p_sym[3]
pole_mass = p_sym[4]
cart_mass = p_sym[5]
wheel_mass = p_sym[6]
wheel_radius = p_sym[7]
wheel_base = p_sym[8]
ground_contact_freq = p_sym[9]
ground_contact_zeta = p_sym[10]
wheel_ground_friction_coeff = p_sym[11]
contact_smoothing_dist = p_sym[12]
friction_smoothing_vel = p_sym[13]

pole_inertia_xx = p_sym[14]
pole_inertia_yy = p_sym[15]
pole_inertia_zz = p_sym[16]

wheel_inertia_xx = p_sym[17]
wheel_inertia_yy = p_sym[18]
wheel_inertia_zz = p_sym[19]

cart_inertia_xx = p_sym[20]
cart_inertia_yy = p_sym[21]
cart_inertia_zz = p_sym[22]

# Define identifiers for generalized coordinates
cart_quat = q_sym[0:4]
cart_pos = q_sym[4:7]
pole_theta = q_sym[7]
lwheel_theta = q_sym[8]
rwheel_theta = q_sym[9]

# Define identifiers for generalized coordinate derivatives
cart_quat_dot = qd_sym[0:4]
cart_pos_dot = qd_sym[4:7]
pole_theta_dot = qd_sym[7]
lwheel_theta_dot = qd_sym[8]
rwheel_theta_dot = qd_sym[9]

# Define identifiers for generalized speeds
cart_ang_vel = u_sym[0:3]
cart_vel = u_sym[3:6]
pole_omega = u_sym[6]
lwheel_omega = u_sym[7]
rwheel_omega = u_sym[8]

# Define identifiers for force inputs
lwheel_motor_torque = in_sym[0]
rwheel_motor_torque = in_sym[1]

# Define newtonian reference frame and origin, gravity
N = ReferenceFrame('N')
O = Point('O')
O.set_vel(N, 0)
gravity = g*N.z

# Define lists of kinematic equations, forces and bodies that will be appended to
kdes = []
forces = []
bodies = []

# Define contact model constants
total_mass = pole_mass+cart_mass+wheel_mass*2
ground_contact_m = total_mass/2
ground_contact_k = ground_contact_m*(ground_contact_freq * 2.*pi)**2
ground_contact_c = 2.*ground_contact_zeta*sqrt(ground_contact_k*ground_contact_m)

# Define cart
cart_frame = ReferenceFrame('cart_frame')
cart_frame.orient(N, 'Quaternion', cart_quat)
cart_frame.set_ang_vel(N, vector_in_frame(cart_frame, cart_ang_vel))
cart_inertia = inertia(cart_frame, cart_inertia_xx, cart_inertia_yy, cart_inertia_zz)
cart_masscenter = O.locatenew('cart_masscenter', vector_in_frame(N, cart_pos))
cart_masscenter.set_vel(N, vector_in_frame(N, cart_vel))
cart_body = RigidBody('cart_body', cart_masscenter, cart_frame, cart_mass, (cart_inertia, cart_masscenter))

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
pole_inertia = inertia(pole_frame, pole_inertia_xx, pole_inertia_yy, pole_inertia_zz)
pole_masscenter = cart_masscenter.locatenew('pole_masscenter', -0.5*pole_length*pole_frame.z)
pole_masscenter.v2pt_theory(cart_masscenter, N, pole_frame)
pole_body = RigidBody('pole_body', pole_masscenter, pole_frame, pole_mass, (pole_inertia, pole_masscenter))
pole_top_point = pole_masscenter.locatenew('pole_top_point', -0.5*pole_length*pole_frame.z)
pole_top_point.v2pt_theory(pole_masscenter,N,pole_frame)
pole_top_pos = pole_top_point.pos_from(O).to_matrix(N)
pole_top_vel = pole_top_point.vel(N).to_matrix(N).subs(pole_theta_dot, pole_omega) # TODO: Why is this subs necessary?
pole_top_ground_normal_force = (-ground_contact_k*pole_top_pos[2] + -ground_contact_c*pole_top_vel[2])*contact(pole_top_pos[2], contact_smoothing_dist)

# Define suspension constants
pole_suspension_m = (pole_inertia+inertia_of_point_mass(pole_mass, pole_masscenter.pos_from(cart_masscenter), pole_frame)).dot(pole_frame.x).to_matrix(pole_frame)[0]
pole_suspension_k = pole_suspension_m*(pole_suspension_freq * 2.*pi)**2
pole_suspension_c = 2.*pole_suspension_zeta*sqrt(pole_suspension_k*pole_suspension_m)

# Add pole to eqns
forces.append((pole_masscenter, pole_mass*gravity))
forces.append((pole_frame, (-pole_suspension_k*pole_theta + -pole_suspension_c*pole_omega)*pole_frame.x))
forces.append((cart_frame, -(-pole_suspension_k*pole_theta + -pole_suspension_c*pole_omega)*pole_frame.x))
forces.append((pole_top_point, N.z*pole_top_ground_normal_force))
kdes.append(pole_theta_dot-pole_omega)
bodies.append(pole_body)

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
lwheel_contact_vel = lwheel_contact_point.vel(N).to_matrix(N).subs(lwheel_theta_dot, lwheel_omega) # TODO: Why is this subs necessary?
lwheel_contact_force = zeros(3)
lwheel_contact_force[2] = (-ground_contact_k*lwheel_contact_pos[2] + -ground_contact_c*lwheel_contact_vel[2])*contact(lwheel_contact_pos[2], contact_smoothing_dist)
lwheel_contact_force[0:2,0] = traction_force(wheel_ground_friction_coeff, lwheel_contact_vel[0:2,0], lwheel_contact_force[2], friction_smoothing_vel)

# Add lwheel to eqns
forces.append((lwheel_masscenter, wheel_mass*gravity))
forces.append((lwheel_contact_point, vector_in_frame(N, lwheel_contact_force)))
forces.append((lwheel_frame, lwheel_motor_torque*lwheel_frame.y))
forces.append((cart_frame, -lwheel_motor_torque*lwheel_frame.y))
kdes.append(lwheel_theta_dot-lwheel_omega)
bodies.append(lwheel_body)

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
rwheel_contact_vel = rwheel_contact_point.vel(N).to_matrix(N).subs(rwheel_theta_dot, rwheel_omega) # TODO: Why is this subs necessary?
rwheel_contact_force = zeros(3)
rwheel_contact_force[2] = (-ground_contact_k*rwheel_contact_pos[2] + -ground_contact_c*rwheel_contact_vel[2])*contact(rwheel_contact_pos[2], contact_smoothing_dist)
rwheel_contact_force[0:2,0] = traction_force(wheel_ground_friction_coeff, rwheel_contact_vel[0:2,0], rwheel_contact_force[2], friction_smoothing_vel)

# Add rwheel to eqns
forces.append((rwheel_masscenter, wheel_mass*gravity))
forces.append((rwheel_contact_point, vector_in_frame(N, rwheel_contact_force)))
forces.append((rwheel_frame, rwheel_motor_torque*rwheel_frame.y))
forces.append((cart_frame, -rwheel_motor_torque*rwheel_frame.y))
kdes.append(rwheel_theta_dot-rwheel_omega)
bodies.append(rwheel_body)

# Set constants
constants = {
        g: 9.80665,
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
        wheel_ground_friction_coeff: 0.5,
        contact_smoothing_dist: 1e-5,
        friction_smoothing_vel: 1e-5,
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

undefconsts = set(p_sym)-set(constants.keys())
assert not undefconsts, 'undefined constants %s' % (undefconsts,)

KM = KanesMethod(N, q_ind=q_sym, u_ind=u_sym, kd_eqs=kdes)
KM.kanes_equations(forces, bodies)

t_sim = 100
bb_sys = System(KM, times=np.linspace(0,t_sim,t_sim*1000))

bb_sys.constants=constants

bb_sys.initial_conditions = {
        cart_quat[0]: 1.,
        cart_quat[1]: 0.,
        cart_quat[2]: 0.,
        cart_quat[3]: 0.,
        cart_pos[2]: -wheel_radius.xreplace(constants),
        cart_vel[0]: 0.,
        cart_vel[1]: 0.,
        cart_vel[2]: 0.,
        pole_theta: 0.,
        pole_omega: 0.,
        lwheel_omega: 0.,
        rwheel_omega: 0.,
    }

bb_sys.specifieds = {
        lwheel_motor_torque: 0.,
        rwheel_motor_torque: 0.,
    }

bb_sys.generate_ode_function(generator='cython')

dyn = bb_sys.evaluate_ode_function
x0 = bb_sys._initial_conditions_padded_with_defaults()
x0 = [x0[k] for k in bb_sys.states]

# Example integration
#result = odeint(dyn,x0,bb_sys.times,(bb_sys._specifieds_padded_with_defaults(), bb_sys._constants_padded_with_defaults()))
