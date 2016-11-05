from helpers import *

# Gravity force
g = 9.80665

# Pole consants
pole_radius = 0.017145
pole_wall_thickness = 0.00254
pole_suspension_freq = 1.5
pole_suspension_zeta = 0.5
pole_length = 1.83
pole_mass = 0.39
pole_inertia_xx_yy = tube_inertia_xx_yy(pole_mass, pole_length, pole_radius-pole_wall_thickness, pole_radius)
pole_inertia_zz = tube_inertia_zz(pole_mass, pole_length, pole_radius-pole_wall_thickness, pole_radius)

# Cart constants
cart_radius = .08
cart_width = .202
cart_mass = 0.6
cart_inertia_xx = tube_inertia_xx_yy(cart_mass, cart_width, 0, cart_radius)
cart_inertia_yy = tube_inertia_zz(cart_mass, cart_width, 0, cart_radius)
cart_inertia_zz = tube_inertia_xx_yy(cart_mass, cart_width, 0, cart_radius)

# Wheel constants
wheel_mass = 0.2
wheel_radius = .092
wheel_base = cart_width+.0255
wheel_inertia_xx_zz = 0.000804/2
wheel_inertia_yy = 0.000804

# Friction model constants
wheel_ground_mu_s = 1. # Static coefficient of friction
wheel_ground_mu_k = .8 # kinetic coefficient of friction
friction_smoothing_vel = 5e-5 # Velocity over which friction fades in. Larger values are more numerically stable.

# Contact model constants
total_mass = pole_mass+cart_mass+wheel_mass*2
wheel_ground_contact_freq = 15.
wheel_ground_contact_zeta = 0.25
wheel_ground_contact_m = total_mass/2 # Assume each wheel bears half of the mass of the vehicle
ground_contact_k = wheel_ground_contact_m*(wheel_ground_contact_freq * 2.*pi)**2
ground_contact_c = 2.*wheel_ground_contact_zeta*sqrt(ground_contact_k*wheel_ground_contact_m)
contact_smoothing_dist = 5e-5 # Distance over which ground contact behavior fades in. Larger values are more numerically stable.
