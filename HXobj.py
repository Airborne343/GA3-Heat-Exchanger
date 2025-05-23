import numpy as np

class HeatExchanger:
    def __init__ (self, length: int, pitch, tube_count: int, baffle_count: int, type: str, passes: int, N_shell: int, baffle_height = 0.7, bundle_height = 0.7, rows = 4):
        # Required attributes
        self.tube_count= tube_count
        self.baffle_count= baffle_count
        self.type = type
        self.passes = passes
        self.N_shell = N_shell
        self.length = length
        self.pitch = pitch
        self.baffle_height = baffle_height
        self.bundle_height = bundle_height
        self.rows = rows

    # Optional attributes with default values
        self.tube_ID = 6/1000 # in meters
        self.tube_OD = 8/1000  # in meters
        self.D_shell = 0.064 # in meters
        self.baffle_width = self.length/(self.baffle_count + 1)

        self.area_nozzle = np.pi * (0.01)**2 # in meters
        self.area_shell = (self.D_shell/self.pitch)*(self.baffle_width)*(self.pitch - self.tube_OD) / self.N_shell
        self.area_tube = 0.25 * self.tube_ID**2 * np.pi
        self.area_tubes = self.tube_count * self.area_tube
        self.area_pipe = (0.25* self.D_shell**2 * np.pi)
        self.inner_surface_area = np.pi * self.tube_ID * self.length
        self.outer_surface_area = np.pi * self.tube_OD * self.length
        
        self.charc_D_shell = self.D_shell * (self.area_shell/self.area_pipe)
        self.sigma = self.area_tubes / self.area_pipe

        if self.baffle_height < 0.5 - self.bundle_height/2:
            self.crossflow_prop = 0
        elif self.baffle_height > 0.5 + self.bundle_height/2:
            self.crossflow_prop = 1
        else:
            self.crossflow_prop = (self.baffle_height - (0.5-self.bundle_height/2)) / (2* self.bundle_height/2)
 
        if self.type.lower() == "60":
            self.a = 0.20
            self.c = 0.20
            #outer conv coeff
            self.A = 0.627615
            self.B = -0.690640
            self.C = -0.0507472
            self.D = 0.0141049
            self.E = -0.000937714

        elif self.type.lower() == "0" or self.type.lower() == "45":
            self.a = 0.34
            self.c = 0.15
            #outer conv coeff
            self.A = 0.374177
            self.B = -0.671577
            self.C = -0.0784051
            self.D = 0.02507191
            self.E = -0.00145983
        else:
            raise ValueError("Invalid heat exchanger type. Use 'triangle' or 'square'.")
        
    #Fluid Constants
        self.heat_cap = 4179 #J/kg K
        self.density = 990.1 #kg/m^3
        self.water_heat_conductivity = 0.632 #W/mK
        self.tube_heat_conductivity = 386 #W/mK
        self.dynamic_viscosity = 6.51 * (10**(-4)) #kg/ms
        self.Prandtl_no = 4.31 #non-dim constant

    #Temperatures
        self.temp_cold = 20 #degree Celsius
        self.temp_hot = 60 #degree Celsius

    def summary(self):
        print(f"Heat Exchanger Summary:")
        print(f"Type: {self.type}")
        print(f"Tube count: {self.tube_count}")
        print(f"Baffle count: {self.baffle_count}")
        print(f"Tube outer diameter: {self.tube_OD} m")
        print(f"Tube inner diameter: {self.tube_ID} m")
        print(f"Tube length: {self.length} m")
        print(f"Shell diameter: {self.D_shell} m")