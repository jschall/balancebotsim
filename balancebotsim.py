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
p_sym = [] # Constants
q_sym = [] # Generalized coordinates
qd_sym = [] # Generalized coordinate derivatives
u_sym = [] # Generalized speeds
in_sym = [] # Inputs

def new_p(n=1):
    return new_sym('p', p_sym, n)
def new_q(n=1):
    return (new_sym('q', q_sym, n=n, dynlevel=0), new_sym('q', qd_sym, n=n, dynlevel=1))
def new_u(n=1):
    return new_sym('u', u_sym, n, dynlevel=0)
def new_in(n=1):
    return new_sym('in', in_sym, n, dynlevel=0)

# Define identifiers for constants
g = new_p()
pole_suspension_freq = new_p()
pole_suspension_zeta = new_p()
pole_length = new_p()
pole_mass = new_p()
cart_mass = new_p()
wheel_mass = new_p()
wheel_radius = new_p()
wheel_base = new_p()
#wheel_bearing_friction_coeff = new_p()
#wheel_bearing_radius = new_p()
ground_contact_freq = new_p()
ground_contact_zeta = new_p()
wheel_ground_friction_coeff = new_p()
contact_smoothing_dist = new_p()
friction_smoothing_vel = new_p()

pole_inertia_xx = new_p()
pole_inertia_yy = new_p()
pole_inertia_zz = new_p()

wheel_inertia_xx = new_p()
wheel_inertia_yy = new_p()
wheel_inertia_zz = new_p()

cart_inertia_xx = new_p()
cart_inertia_yy = new_p()
cart_inertia_zz = new_p()

# Define identifiers for generalized coordinates and derivatives
cart_quat, cart_quat_dot = new_q(4)
cart_pos, cart_pos_dot = new_q(3)
pole_theta, pole_theta_dot = new_q()
lwheel_theta, lwheel_theta_dot = new_q()
rwheel_theta, rwheel_theta_dot = new_q()

# Define identifiers for generalized speeds
cart_ang_vel = new_u(3)
cart_vel = new_u(3)
pole_omega = new_u()
lwheel_omega = new_u()
rwheel_omega = new_u()

# Define identifiers for force inputs
lwheel_motor_torque = new_in()
rwheel_motor_torque = new_in()

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

        cart_inertia_xx: tube_inertia_xx_yy(cart_mass, wheel_base-0.03, wheel_radius*0.3, wheel_radius*0.9).xreplace(constants),
        cart_inertia_yy: tube_inertia_zz(cart_mass, wheel_base-0.03, wheel_radius*0.3, wheel_radius*0.9).xreplace(constants),
        cart_inertia_zz: tube_inertia_xx_yy(cart_mass, wheel_base-0.03, wheel_radius*0.3, wheel_radius*0.9).xreplace(constants)
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

get_cart_pos = lambdify([q_sym+u_sym], (cart_masscenter.pos_from(O)-0.5*(wheel_base-0.03)*cart_frame.y).to_matrix(N).xreplace(constants))
get_cart_axis = lambdify([q_sym+u_sym], ((wheel_base-0.03)*cart_frame.y).to_matrix(N).xreplace(constants))

get_pole_pos = lambdify([q_sym+u_sym], (pole_masscenter.pos_from(O)-0.5*pole_length*pole_frame.z).to_matrix(N).xreplace(constants))
get_pole_axis = lambdify([q_sym+u_sym], (pole_length*pole_frame.z).to_matrix(N).xreplace(constants))

get_lwheel_pos = lambdify([q_sym+u_sym], (lwheel_masscenter.pos_from(O)-0.5*.02*cart_frame.y).to_matrix(N).xreplace(constants))
get_lwheel_axis = lambdify([q_sym+u_sym], (.02*cart_frame.y).to_matrix(N).xreplace(constants))

get_rwheel_pos = lambdify([q_sym+u_sym], (rwheel_masscenter.pos_from(O)-0.5*.02*cart_frame.y).to_matrix(N).xreplace(constants))
get_rwheel_axis = lambdify([q_sym+u_sym], (.02*cart_frame.y).to_matrix(N).xreplace(constants))

from multiprocessing import Process, Queue
from visual import *

def vis_proc(q):
    def vpy(v):
        return vector(v[1], -v[2], -v[0])

    scene = display(width=1024, height=1024)

    theta = radians(30.)
    scene.forward=vpy((cos(theta),0.,sin(theta)))

    floor = cylinder(pos=(0,0,0), axis=vpy((0.,0.,0.1)), material=materials.wood, radius=5.)

    x = q.get()

    cart = cylinder(pos=vpy(get_cart_pos(x)), axis=vpy(get_cart_axis(x)), radius=0.9*wheel_radius.xreplace(constants))
    lwheel = cylinder(pos=vpy(get_lwheel_pos(x)), axis=vpy(get_lwheel_axis(x)), radius=wheel_radius.xreplace(constants))
    rwheel = cylinder(pos=vpy(get_rwheel_pos(x)), axis=vpy(get_rwheel_axis(x)), radius=wheel_radius.xreplace(constants))
    pole = cylinder(pos=vpy(get_pole_pos(x)), axis=vpy(get_pole_axis(x)), radius=0.03429/2)

    while(True):
        rate(60)
        x = q.get()

        cart.pos = vpy(get_cart_pos(x))
        cart.axis = vpy(get_cart_axis(x))

        lwheel.pos = vpy(get_lwheel_pos(x))
        lwheel.axis = vpy(get_lwheel_axis(x))

        rwheel.pos = vpy(get_rwheel_pos(x))
        rwheel.axis = vpy(get_rwheel_axis(x))

        pole.pos = vpy(get_pole_pos(x))
        pole.axis = vpy(get_pole_axis(x))


if __name__ == '__main__':
    import signal
    import sys

    q = Queue(3)
    p = Process(target=vis_proc, args=(q,))
    t = 0.
    dt = 1./60.
    x = x0
    q.put(x)
    p.start()

    def signal_handler(signal, frame):
        p.terminate()
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)

    while(True):
        times = np.linspace(t,t+dt,2)
        x = odeint(dyn,x,times,(bb_sys._specifieds_padded_with_defaults(), bb_sys._constants_padded_with_defaults()), rtol=1e-5, atol=1e-11)[-1]
        t = times[-1]
        q.put(x)
