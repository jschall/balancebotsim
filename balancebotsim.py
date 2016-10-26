from helpers import *
from sympy import *
from sympy.physics.mechanics import *

# Define constants
pole_suspension_omega = 1.0
pole_suspension_zeta = 0.4

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
pole_length = Symbol('pole_length')
pole_mass = Symbol('pole_mass')
pole_theta = dynamicsymbols('pole_theta', 0)
pole_theta_dot = dynamicsymbols('pole_theta', 1)
pole_omega = dynamicsymbols('pole_omega')


# Define wheel symbols
wheel_mass = Symbol('wheel_mass')
wheel_radius = Symbol('wheel_radius')
wheel_base = Symbol('wheel_base')

# Define lwheel symbols
lwheel_theta = dynamicsymbols('lwheel_theta', 0)
lwheel_theta_dot = dynamicsymbols('lwheel_theta', 1)
lwheel_omega = dynamicsymbols('lwheel_omega')
lwheel_contact_force = Matrix(dynamicsymbols('lwheel_contact_force[0:3]'))
lwheel_motor_torque = dynamicsymbols('lwheel_motor_torque')

# Define rwheel symbols
rwheel_theta = dynamicsymbols('rwheel_theta', 0)
rwheel_theta_dot = dynamicsymbols('rwheel_theta', 1)
rwheel_omega = dynamicsymbols('rwheel_omega')
rwheel_contact_force = Matrix(dynamicsymbols('rwheel_contact_force[0:3]'))
rwheel_motor_torque = dynamicsymbols('rwheel_motor_torque')

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
pole_masscenter = cart_masscenter.locatenew('pole_masscenter', -0.5*pole_length*pole_frame.z)
pole_masscenter.v2pt_theory(cart_masscenter, N, pole_frame)
pole_body = RigidBody('pole_body', pole_masscenter, pole_frame, pole_mass, (pole_inertia, pole_masscenter))

# Define suspension constants
pole_suspension_m = (pole_inertia+inertia_of_point_mass(pole_mass, pole_masscenter.pos_from(cart_masscenter), pole_frame)).dot(pole_frame.x).to_matrix(pole_frame)[0]
pole_suspension_k = pole_suspension_m*pole_suspension_omega**2
pole_suspension_c = 2.*pole_suspension_zeta*sqrt(pole_suspension_k*pole_suspension_m)

# Add pole to eqns
q_list.append(pole_theta)
u_list.append(pole_omega)
force_list.append((pole_masscenter, pole_mass*gravity*N.z))
force_list.append((pole_frame, (-pole_suspension_k*pole_theta + -pole_suspension_c*pole_omega)*pole_frame.x))
force_list.append((cart_frame, -(-pole_suspension_k*pole_theta + -pole_suspension_c*pole_omega)*pole_frame.x))
kde_list.append(pole_theta_dot-pole_omega)
body_list.append(pole_body)

# Define the ground direction vector
# This is done by subtracting the component of the N.z vector along cart_frame.y from N.z, resulting in only the component in the cart_frame x-z plane, and normalizing
ground_direction_vector = (N.z - N.z.dot(cart_frame.y)*cart_frame.y) / N.z.dot(cart_frame.y)

# Define lwheel
lwheel_frame = ReferenceFrame('lwheel_frame')
lwheel_frame.orient(cart_frame, 'Axis', [lwheel_theta, cart_frame.y])
lwheel_inertia = inertia(lwheel_frame, 1., 1., 1.)
lwheel_masscenter = cart_masscenter.locatenew('lwheel_masscenter', -0.5*wheel_base*cart_frame.y)
lwheel_masscenter.v2pt_theory(cart_masscenter, N, lwheel_frame)
lwheel_body = RigidBody('lwheel_body', lwheel_masscenter, lwheel_frame, wheel_mass, (lwheel_inertia, lwheel_masscenter))
lwheel_contact_point = lwheel_masscenter.locatenew('lwheel_contact_point', ground_direction_vector*wheel_radius)
lwheel_contact_point.v2pt_theory(lwheel_masscenter,N,lwheel_frame)

# Add lwheel to eqns
q_list.append(lwheel_theta)
u_list.append(lwheel_omega)
force_list.append((lwheel_masscenter, wheel_mass*gravity*N.z))
force_list.append((lwheel_contact_point, vector_in_frame(N,lwheel_contact_force)))
force_list.append((lwheel_frame, lwheel_motor_torque*lwheel_frame.y))
force_list.append((cart_frame, -lwheel_motor_torque*lwheel_frame.y))
kde_list.append(lwheel_theta_dot-lwheel_omega)
body_list.append(lwheel_body)

# Define rwheel
rwheel_frame = ReferenceFrame('rwheel_frame')
rwheel_frame.orient(cart_frame, 'Axis', [rwheel_theta, cart_frame.y])
rwheel_inertia = inertia(rwheel_frame, 1., 1., 1.)
rwheel_masscenter = cart_masscenter.locatenew('rwheel_masscenter', 0.5*wheel_base*cart_frame.y)
rwheel_masscenter.v2pt_theory(cart_masscenter, N, rwheel_frame)
rwheel_body = RigidBody('rwheel_body', rwheel_masscenter, rwheel_frame, wheel_mass, (rwheel_inertia, rwheel_masscenter))
rwheel_contact_point = rwheel_masscenter.locatenew('rwheel_contact_point', ground_direction_vector*wheel_radius)
rwheel_contact_point.v2pt_theory(rwheel_masscenter,N,rwheel_frame)

# Add rwheel to eqns
q_list.append(rwheel_theta)
u_list.append(rwheel_omega)
force_list.append((rwheel_masscenter, wheel_mass*gravity*N.z))
force_list.append((rwheel_contact_point, vector_in_frame(N,rwheel_contact_force)))
force_list.append((rwheel_frame, rwheel_motor_torque*rwheel_frame.y))
force_list.append((cart_frame, -rwheel_motor_torque*rwheel_frame.y))
kde_list.append(rwheel_theta_dot-rwheel_omega)
body_list.append(rwheel_body)

KM = KanesMethod(N, q_ind=q_list, u_ind=u_list, kd_eqs=kde_list)
KM.kanes_equations(force_list, body_list)
kdd = KM.kindiffdict()
mm = KM.mass_matrix_full
fo = KM.forcing_full

rwheel_contact_pos = Matrix(simplify(rwheel_contact_point.pos_from(O).to_matrix(N)))
rwheel_contact_vel = Matrix(simplify(rwheel_contact_point.vel(N).to_matrix(N)))
lwheel_contact_pos = Matrix(simplify(lwheel_contact_point.pos_from(O).to_matrix(N)))
lwheel_contact_vel = Matrix(simplify(lwheel_contact_point.vel(N).to_matrix(N)))

print(".")
subx, outx = cse([mm, fo, rwheel_contact_pos, rwheel_contact_vel, lwheel_contact_pos, lwheel_contact_vel])
mm, fo, rwheel_contact_pos, rwheel_contact_vel, lwheel_contact_pos, lwheel_contact_vel = tuple(outx)

print(".")
mm = Matrix(simplify(mm))
fo = Matrix(simplify(fo))

print(".")
eom = mm.LUsolve(fo)

print(".")
eom_subx, outx = cse([eom])
eom = outx[0]
subx = subx+eom_subx
print(count_ops(eom)+count_ops(subx))

with open('out.srepr', 'wb') as f:
    f.truncate()
    f.write(srepr({'q_list':q_list, 'u_list':u_list, 'subx':subx, 'eom':eom, 'rwheel_contact_pos':rwheel_contact_pos, 'rwheel_contact_vel':rwheel_contact_vel, 'lwheel_contact_pos':lwheel_contact_pos, 'lwheel_contact_vel':lwheel_contact_vel}).encode('utf8'))




