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
lwheel_contact_p = Matrix(symbols('lwheel_contact_p[0:3]'))
lwheel_contact_f = Matrix(symbols('lwheel_contact_f[0:3]'))

# Define rwheel symbols
rwheel_mass = Symbol('rwheel_mass')
rwheel_theta = dynamicsymbols('rwheel_theta', 0)
rwheel_theta_dot = dynamicsymbols('rwheel_theta', 1)
rwheel_omega = dynamicsymbols('rwheel_omega')

# Define cart
cart_frame = ReferenceFrame('cart_frame')
cart_frame.orient(N, 'Quaternion', cart_quat)
cart_frame.set_ang_vel(N, vector_in_frame(N, cart_ang_vel))
cart_inertia = inertia(cart_frame, 1., 1., 1.)
cart_masscenter = O.locatenew('cart_masscenter', vector_in_frame(N, cart_pos))
cart_masscenter.set_vel(N, vector_in_frame(N, cart_vel))
cart_body = RigidBody('cart_body', cart_masscenter, cart_frame, cart_mass, (cart_inertia, cart_masscenter))

# Add cart to eqns
q_list.extend(list(cart_pos))
u_list.extend(list(cart_vel))
q_list.extend(list(cart_quat))
u_list.extend(list(cart_ang_vel))
force_list.append((cart_masscenter, cart_mass*gravity*N.z))
kde_list.extend(list(cart_pos_dot-cart_vel)+list(cart_quat_dot-quatderiv(cart_quat, cart_ang_vel)))
body_list.append(cart_body)

# Define pole
pole_frame = ReferenceFrame('pole_frame')
pole_frame.orient(cart_frame, 'Axis', [pole_theta, cart_frame.x])
pole_inertia = inertia(pole_frame, 1., 1., 1.)
pole_masscenter = cart_masscenter.locatenew('pole_masscenter', -1.*pole_frame.z)
pole_masscenter.v2pt_theory(cart_masscenter, N, pole_frame)
pole_body = RigidBody('pole_body', pole_masscenter, pole_frame, pole_mass, (pole_inertia, pole_masscenter))

# Add pole to eqns
q_list.append(pole_theta)
u_list.append(pole_omega)
force_list.append((pole_masscenter, pole_mass*gravity*N.z))
kde_list.append(pole_theta_dot-pole_omega)
body_list.append(pole_body)

# Define lwheel
lwheel_frame = ReferenceFrame('lwheel_frame')
lwheel_frame.orient(cart_frame, 'Axis', [lwheel_theta, cart_frame.y])
lwheel_inertia = inertia(lwheel_frame, 1., 1., 1.)
lwheel_masscenter = cart_masscenter.locatenew('lwheel_masscenter', -0.12*cart_frame.y)
lwheel_masscenter.v2pt_theory(cart_masscenter, N, lwheel_frame)
lwheel_body = RigidBody('lwheel_body', lwheel_masscenter, lwheel_frame, lwheel_mass, (lwheel_inertia, lwheel_masscenter))
lwheel_contact_point = lwheel_masscenter.locatenew('lwheel_contact_point', vector_in_frame(lwheel_frame, lwheel_contact_p))
lwheel_contact_point.v2pt_theory(lwheel_masscenter, N, lwheel_frame)
lwheel_contact_force = vector_in_frame(lwheel_frame, lwheel_contact_f)

# Add lwheel to eqns
q_list.append(lwheel_theta)
u_list.append(lwheel_omega)
force_list.append((lwheel_masscenter, lwheel_mass*gravity*N.z))
force_list.append((lwheel_contact_point, lwheel_contact_force))
kde_list.append(lwheel_theta_dot-lwheel_omega)
body_list.append(lwheel_body)


KM = KanesMethod(N, q_ind=q_list, u_ind=u_list, kd_eqs=kde_list)
KM.kanes_equations(force_list, body_list)
kdd = KM.kindiffdict()
mm = KM.mass_matrix_full
fo = KM.forcing_full

subx, outx = cse([mm,fo])
mm = outx[0]
fo = outx[1]
mm = Matrix(simplify(mm))
fo = Matrix(simplify(fo))
print(count_ops(mm) + count_ops(fo) + count_ops(subx))
eom = mm.LUsolve(fo)
#print(simplify(eom))
eom_subx, eom = cse(eom)
subx = subx+eom_subx
print(count_ops(eom)+count_ops(subx))

#subx,eom = cse(eom)
#print(count_ops(eom)+count_ops(subx))
#pprint(eom)
#mechanics_printing()
#mprint(eom)





