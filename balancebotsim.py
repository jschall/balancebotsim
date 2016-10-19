from helpers import *
from sympy import *
from sympy.physics.mechanics import *

kde_list = []
force_list = []
body_list = []
q_list = []
u_list = []

N = ReferenceFrame('N')
O = Point('O')
O.set_vel(N, 0)
gravity = Symbol('g')

# Define cart symbols
cart_mass = Symbol('cart_mass')
cart_pos = Matrix(dynamicsymbols('cart_pos[0:3]', 0))
cart_pos_dot = Matrix(dynamicsymbols('cart_pos[0:3]', 1))
cart_vel = Matrix(dynamicsymbols('cart_vel[0:3]', 0))
cart_quat = Matrix(dynamicsymbols('cart_quat[0:4]', 0))
cart_quat_dot = Matrix(dynamicsymbols('cart_quat[0:4]', 1))
cart_ang_vel = Matrix(dynamicsymbols('cart_omega[0:3]'))

# Define pole symbols
pole_mass = Symbol('pole_mass')
pole_theta = dynamicsymbols('pole_theta', 0)
pole_theta_dot = dynamicsymbols('pole_theta', 1)
pole_omega = dynamicsymbols('pole_omega')

# Define lwheel symbols
lwheel_mass = Symbol('lwheel_mass')
lwheel_theta = dynamicsymbols('lwheel_theta', 0)
lwheel_theta_dot = dynamicsymbols('lwheel_theta', 1)
lwheel_omega = dynamicsymbols('lwheel_omega')
lwheel_contact_force = Symbol('lwheel_contact_force')
lwheel_friction_force = Symbol('lwheel_friction_force')

# Define rwheel symbols
rwheel_mass = Symbol('rwheel_mass')
rwheel_theta = dynamicsymbols('rwheel_theta', 0)
rwheel_theta_dot = dynamicsymbols('rwheel_theta', 1)
rwheel_omega = dynamicsymbols('rwheel_omega')
rwheel_contact_force = Symbol('rwheel_contact_force')
rwheel_friction_force = Symbol('rwheel_friction_force')

# Define cart
cart_frame = ReferenceFrame('cart_frame')
cart_frame.orient(N, 'Quaternion', cart_quat)
cart_frame.set_ang_vel(N, vector_in_frame(N, cart_ang_vel))
cart_inertia = inertia(cart_frame, 1., 1., 1.)
cart_masscenter = O.locatenew('cart_masscenter', vector_in_frame(N, cart_pos))
cart_masscenter.set_vel(N, vector_in_frame(N, cart_vel))
cart_body = RigidBody('cart_body', cart_masscenter, cart_frame, 1., (cart_inertia, cart_masscenter))

# Add cart to eqns
q_list.extend(list(cart_pos))
u_list.extend(list(cart_vel))
q_list.extend(list(cart_quat))
u_list.extend(list(cart_ang_vel))
force_list.append((cart_masscenter, 1.*gravity*N.z))
kde_list.extend(list(cart_pos_dot-cart_vel)+list(cart_quat_dot-quatderiv(cart_quat, cart_ang_vel)))
body_list.append(cart_body)

# Define pole
pole_frame = ReferenceFrame('pole_frame')
pole_frame.orient(cart_frame, 'Axis', [pole_theta, cart_frame.x])
pole_frame.set_ang_vel(cart_frame, pole_omega*cart_frame.x)
pole_inertia = inertia(pole_frame, 1., 1., 1.)
pole_masscenter = cart_masscenter.locatenew('pole_masscenter', -1.*pole_frame.z)
pole_masscenter.v2pt_theory(cart_masscenter, N, pole_frame)
pole_body = RigidBody('pole_body', pole_masscenter, pole_frame, 1., (pole_inertia, pole_masscenter))

# Add pole to eqns
q_list.append(pole_theta)
u_list.append(pole_omega)
force_list.append((pole_masscenter, 1.*gravity*N.z))
print(cart_frame.ang_vel_in(N).dot(cart_frame.x))
kde_list.append(pole_theta_dot-pole_omega)
body_list.append(pole_body)
print(len(kde_list), len(q_list), len(u_list))


# Define lwheel
#lwheel_frame = ReferenceFrame('lwheel_frame')
#lwheel_frame.orient(cart_frame, 'Axis', [lwheel_theta, cart_frame.y])
#lwheel_frame.set_ang_vel(lwheel_frame, lwheel_omega*lwheel_frame.y)
#lwheel_inertia = inertia(lwheel_frame, 1., 1., 1.)
#lwheel_masscenter = cart_masscenter.locatenew('lwheel_masscenter', -0.12*lwheel_frame.y)
#lwheel_masscenter.v2pt_theory(cart_masscenter, N, lwheel_frame)
#lwheel_body = RigidBody('lwheel_body', lwheel_masscenter, lwheel_frame, lwheel_mass, (lwheel_inertia, lwheel_masscenter))

KM = KanesMethod(N, q_ind=q_list, u_ind=u_list, kd_eqs=kde_list)
KM.kanes_equations(force_list, body_list)
kdd = KM.kindiffdict()
mm = KM.mass_matrix_full
fo = KM.forcing_full
print(count_ops(mm))
print(count_ops(fo))
mm = Matrix(simplify(mm))
fo = Matrix(simplify(fo))
print(count_ops(mm))
print(count_ops(fo))
eom = mm.LUsolve(fo)
print(count_ops(eom))
subx,eom = cse(eom)
print(count_ops(eom)+count_ops(subx))
pprint(eom)
#mechanics_printing()
#mprint(eom)





