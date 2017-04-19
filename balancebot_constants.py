from helpers import *

# Gravity acceleration
g = 9.80655

# Pole consants
pole_radius = 0.017145        # from datasheet
pole_wall_thickness = 0.00254 # from datasheet
pole_suspension_freq = 1.5    # guessed
pole_suspension_zeta = 0.5    # guessed
pole_length = 1.83            # from datasheet
pole_mass = 0.39              # from datasheet
pole_inertia_xx_yy = tube_inertia_xx_yy(pole_mass, pole_length, pole_radius-pole_wall_thickness, pole_radius)
pole_inertia_zz = tube_inertia_zz(pole_mass, pole_length, pole_radius-pole_wall_thickness, pole_radius)

# Battery constants
battery_mass = .4
battery_length = .53
battery_radius = .01
battery_inertia_xx_yy = tube_inertia_xx_yy(battery_mass, battery_length, 0, battery_radius)
battery_inertia_zz = tube_inertia_zz(battery_mass, battery_length, 0, battery_radius)

# Payload constants
payload_mass = 0.8
payload_height = .25
payload_width = .19
payload_thickness = .013
payload_position_h = pole_length-payload_height/2.
payload_position_x = pole_radius+payload_thickness/2
payload_inertia_xx = payload_mass/12. * (payload_width**2+payload_height**2)
payload_inertia_yy = payload_mass/12. * (payload_thickness**2+payload_height**2)
payload_inertia_zz = payload_mass/12. * (payload_width**2+payload_thickness**2)

# Cart constants
cart_radius = .08
cart_width = .202
cart_mass = 0.4
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
wheel_ground_mu_s = 1.2 # Static coefficient of friction
wheel_ground_mu_k = 1.0 # kinetic coefficient of friction
friction_smoothing_vel = 1e-4 # Velocity over which friction fades in. Larger values are more numerically stable.

# Contact model constants
total_mass = pole_mass+cart_mass+payload_mass+battery_mass+wheel_mass*2
wheel_ground_contact_freq = 30.
wheel_ground_contact_zeta = 0.3
wheel_ground_contact_m = total_mass/2 # Assume each wheel bears half of the mass of the vehicle
ground_contact_k = wheel_ground_contact_m*(wheel_ground_contact_freq * 2.*pi)**2
ground_contact_c = 2.*wheel_ground_contact_zeta*sqrt(ground_contact_k*wheel_ground_contact_m)
contact_smoothing_dist = 1e-4 # Distance over which ground contact behavior fades in. Larger values are more numerically stable.
