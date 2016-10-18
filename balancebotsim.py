from sympy import *
from sympy.physics.mechanics import *
from helpers import quatderiv, vector_in_frame

N = ReferenceFrame('N')
O = Point('O')
O.set_vel(N,0)
gravity = Symbol('g')

# Define cart symbols
cart_mass = Symbol('cart_mass')
cart_pos = Matrix(dynamicsymbols('cart_pos[0:3]', 0))
cart_pos_dot = Matrix(dynamicsymbols('cart_pos[0:3]', 1))
cart_vel = Matrix(dynamicsymbols('cart_vel[0:3]', 0))
cart_quat = Matrix(dynamicsymbols('cart_quat[0:4]', 0))
cart_quat_dot = Matrix(dynamicsymbols('cart_quat[0:4]', 1))
cart_ang_vel = Matrix(dynamicsymbols('cart_omega[0:3]'))

# Define cart
cart_frame = ReferenceFrame('cart_frame')
cart_frame.orient(N, 'Quaternion', cart_quat)
cart_frame.set_ang_vel(N, vector_in_frame(N, cart_ang_vel))
cart_inertia = inertia(cart_frame, 1., 1., 1.)
cart_masscenter = Point('cart_masscenter')
cart_masscenter.set_pos(O, vector_in_frame(N, cart_pos))
cart_masscenter.set_vel(N, vector_in_frame(N, cart_vel))
cart_body = RigidBody('cart_body', cart_masscenter, cart_frame, cart_mass, (cart_inertia, cart_masscenter))

kinematic_differential_equations = list(cart_pos_dot-cart_vel)+list(cart_quat_dot-quatderiv(cart_quat, cart_ang_vel))
force_list = [(cart_masscenter, cart_mass*gravity*N.z)]
body_list = [cart_body]

KM = KanesMethod(N, q_ind=[cart_pos, cart_quat], u_ind=[cart_vel, cart_ang_vel], kd_eqs=kinematic_differential_equations)
KM.kanes_equations(force_list, body_list)
kdd = KM.kindiffdict()
mm = KM.mass_matrix_full
fo = KM.forcing_full

eom = mm.inv() * fo
eom = eom.subs(kdd)
eom.simplify()
mechanics_printing()
mprint(eom)


## Define pole symbols
#pole_mass = Symbol('pole_mass')
#pole_theta = dynamicsymbols('pole_theta', 0)
#pole_theta_dot = dynamicsymbols('pole_theta', 1)
#pole_omega = dynamicsymbols('pole_omega')

## Define lwheel symbols
#lwheel_mass = Symbol('lwheel_mass')
#lwheel_theta = dynamicsymbols('lwheel_theta', 0)
#lwheel_theta_dot = dynamicsymbols('lwheel_theta', 1)
#lwheel_omega = dynamicsymbols('lwheel_omega')
#lwheel_contact_force = Symbol('lwheel_contact_force')
#lwheel_friction_force = Symbol('lwheel_friction_force')

## Define rwheel symbols
#rwheel_mass = Symbol('rwheel_mass')
#rwheel_theta = dynamicsymbols('rwheel_theta', 0)
#rwheel_theta_dot = dynamicsymbols('rwheel_theta', 1)
#rwheel_omega = dynamicsymbols('rwheel_omega')
#rwheel_contact_force = Symbol('rwheel_contact_force')
#rwheel_friction_force = Symbol('rwheel_friction_force')
