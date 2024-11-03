import numpy as np
import math as math

class SpringModel:
    def __init__(self, height, width, length, inner_radius, outer_radius, density, epsilon_e, epsilon_d, material_modulus, friction_coefficient,brake_dia,theta,F_Brake):
        # Input parameters
        self.height = height
        self.width = width
        self.length = length
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.density = density
        self.epsilon_e = epsilon_e
        self.epsilon_d = epsilon_d
        self.material_modulus = material_modulus
        self.friction_coefficient = friction_coefficient
        self.brake_dia=brake_dia
        self.theta=theta
        self.fbrake= F_Brake

        # Derived parameters
        self.volume = self.calculate_volume()
        self.mass = self.calculate_mass()
        self.inertia = self.calculate_inertia()
        self.torsional_spring_constant = self.calculate_torsional_spring_constant()
        self.time_constant = self.calculate_time_constant()
        
    def calculate_volume(self):
        """ Calculate the volume of the spring. """
        return np.pi * self.height * (self.outer_radius**2 - self.inner_radius**2)

    def calculate_mass(self):
        """ Calculate the mass of the spring. """
        return self.volume * self.density

    def calculate_inertia(self):
        """ Calculate the moment of inertia of the spring. """
        return (1/2) * self.mass * ((self.outer_radius**2 + self.inner_radius**2) / 2)

    def calculate_torsional_spring_constant(self):
        """ Calculate the torsional spring constant. """
        return (self.material_modulus * (self.outer_radius - self.inner_radius) ** 3) / (2 * self.height)

    def calculate_time_constant(self):
        """ Calculate the time constant for the spring. """
        return self.inertia / (self.torsional_spring_constant)
    
    def capstan_ratio(self):
        Capstan_ratio= 2.7183^((-self.friction_coefficient*self.theta*math.pi)/180)
        return Capstan_ratio
   
    def calculate_btf(self):
        CR=self.capstan_ratio()
        return math.sqrt(1+CR^2-(2*CR*math.cos(self.theta*math.pi)/100))
   
    def calculate_forces(self):
        """ Calculate the forces acting on the spring. """
        F_radial = self.friction_coefficient * self.mass * 9.81  # Assuming gravity
        F_theta = 0.5 * F_radial  # Arbitrary for illustration; depends on your application
        F_axial = F_radial
        F_Pull =  self.capstan_ratio()*F_radial
        return F_radial, F_theta, F_axial,F_Pull
    
    def calculate_braketorque(self):
        F_radial, F_theta, F_axial,F_Pull = self.calculate_forces()
        btf= self.calculate_btf
        brake_torque= (self.friction_coefficient*(self.brake_dia/2)+max(btf*self.fbrake,F_Pull))
        return brake_torque
    
    
    def calculate_outputs(self):
        """ Calculate output parameters. """
        F_radial, F_theta, F_axial,F_Pull = self.calculate_forces()
        CapstanR= self.capstan_ratio()
        # Output parameters
        time_to_stop = self.time_constant  # Time to stop
        torque = F_radial * self.outer_radius  # Torque calculation
        freq_switching = 1 / (time_to_stop + self.time_constant)  # Frequency of switching
        reflected_inertia = self.inertia / (1 + (self.friction_coefficient / 2))  # Reflected inertia
        time_to_restore = self.time_constant  # Assume same as time constant
        period = time_to_stop + time_to_restore  # Total period
        bt_f=  self.calculate_btf()
        Brake_torque=  self.calculate_braketorque()
        Net_force= self.theta*CapstanR

        
        return {
            "Time to Stop": time_to_stop,
            "Torque": torque,
            "Force Radial": F_radial,
            "Force Theta": F_theta,
            "Force Axial": F_axial,
            "Capstan Ratio": CapstanR,
            "Force Pull" :  F_Pull,
            "Frequency of Switching": freq_switching,
            "Reflected Inertia": reflected_inertia,
            "Time to Restore": time_to_restore,
            "Total Period": period,
            "Brake Torque":Brake_torque,
            "BTF":bt_f,
            "Net Force Out":Net_force
        }


spring = SpringModel(
    height=0.1,              # Height in meters
    width=0.02,              # Width in meters
    length=0.05,             # Length in meters
    inner_radius=0.005,      # Inner radius in meters
    outer_radius=0.01,       # Outer radius in meters
    density=2700,            # Density in kg/m^3 (e.g., aluminum)
    epsilon_e=0.4,           # Electric permittivity
    epsilon_d=0.5,           # Dielectric constant
    material_modulus=70e9,
    brake_dia=0.1,   
    friction_coefficient=0.2,
    theta= 2638.54,# Angle of wrap
    F_brake= 1  
)

# Calculate output parameters
outputs = spring.calculate_outputs()

# Display results
for key, value in outputs.items():
    print(f"{key}: {value:}")
# Plots based on input parameters 
