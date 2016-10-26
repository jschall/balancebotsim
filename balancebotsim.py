from helpers import *
from sympy import *
from sympy.physics.mechanics import *
from pydy.system import System

# Define constants
gravity = 9.80655
pole_suspension_omega = 2. * 2.*pi
pole_suspension_zeta = 0.4
pole_length = 1.83
pole_mass = 0.39
cart_mass = 1.0
wheel_mass = 0.1
wheel_radius = .1
wheel_base = .2

kde_list = []
force_list = []
body_list = []
q = dynamicsymbols('q0:10')
qd = dynamicsymbols('q0:10',1)
u = dynamicsymbols('u0:9')
f = dynamicsymbols('f0:8')

N = ReferenceFrame('N')
O = Point('O')
O.set_vel(N, 0)

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
lwheel_contact_force = Matrix(f[0:3])
lwheel_motor_torque = f[3]
rwheel_contact_force = Matrix(f[4:7])
rwheel_motor_torque = f[7]


# Define cart
cart_frame = ReferenceFrame('cart_frame')
cart_frame.orient(N, 'Quaternion', cart_quat)
cart_frame.set_ang_vel(N, vector_in_frame(N, cart_ang_vel))
cart_inertia = inertia(cart_frame, 1., 1., 1.)
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
pole_inertia = inertia(pole_frame, 1., 1., 1.)
pole_masscenter = cart_masscenter.locatenew('pole_masscenter', -0.5*pole_length*pole_frame.z)
pole_masscenter.v2pt_theory(cart_masscenter, N, pole_frame)
pole_body = RigidBody('pole_body', pole_masscenter, pole_frame, pole_mass, (pole_inertia, pole_masscenter))

# Define suspension constants
pole_suspension_m = (pole_inertia+inertia_of_point_mass(pole_mass, pole_masscenter.pos_from(cart_masscenter), pole_frame)).dot(pole_frame.x).to_matrix(pole_frame)[0]
pole_suspension_k = pole_suspension_m*pole_suspension_omega**2
pole_suspension_c = 2.*pole_suspension_zeta*sqrt(pole_suspension_k*pole_suspension_m)

# Add pole to eqns
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
force_list.append((rwheel_masscenter, wheel_mass*gravity*N.z))
force_list.append((rwheel_contact_point, vector_in_frame(N,rwheel_contact_force)))
force_list.append((rwheel_frame, rwheel_motor_torque*rwheel_frame.y))
force_list.append((cart_frame, -rwheel_motor_torque*rwheel_frame.y))
kde_list.append(rwheel_theta_dot-rwheel_omega)
body_list.append(rwheel_body)

rwheel_contact_pos = rwheel_contact_point.pos_from(O).to_matrix(N)
rwheel_contact_vel = rwheel_contact_point.vel(N).to_matrix(N)
lwheel_contact_pos = lwheel_contact_point.pos_from(O).to_matrix(N)
lwheel_contact_vel = lwheel_contact_point.vel(N).to_matrix(N)

KM = KanesMethod(N, q_ind=q, u_ind=u, kd_eqs=kde_list)
KM.kanes_equations(force_list, body_list)

f_val = [0.,0.,0.,0.,0.,0.,0.,0.]

s = System(KM, times=[0.,0.01])
s.specifieds = {'symbols':f, 'values': f_val}
s.generate_ode_function(generator='cython')
print(s.integrate())

#kdd = KM.kindiffdict()
#mm = KM.mass_matrix_full
#fo = KM.forcing_full

#mm, fo, subx = extractSubexpressions([mm,fo], 'subx')

#mm = simplify(mm).as_mutable()
#fo = simplify(fo).as_mutable()

#eom = mm.LUsolve(fo)

#eom, rwheel_contact_pos, rwheel_contact_vel, lwheel_contact_pos, lwheel_contact_vel, subx = extractSubexpressions([eom, rwheel_contact_pos, rwheel_contact_vel, lwheel_contact_pos, lwheel_contact_vel], 'subx', prev_subx=subx)

#print(count_ops([eom, rwheel_contact_pos, rwheel_contact_vel, lwheel_contact_pos, lwheel_contact_vel, subx]))

#with open('out.srepr', 'wb') as f:
    #f.truncate()
    #f.write(srepr({'subx':subx, 'eom':eom, 'rwheel_contact_pos':rwheel_contact_pos, 'rwheel_contact_vel':rwheel_contact_vel, 'lwheel_contact_pos':lwheel_contact_pos, 'lwheel_contact_vel':lwheel_contact_vel}).encode('utf8'))




