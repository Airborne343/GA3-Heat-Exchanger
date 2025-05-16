import numpy as np

class HeatExchanger:
    def __init__ (self, tube_count: int, baffle_count: int, type: str):
        # Required attributes
        self.tube_count= tube_count
        self.baffle_count= baffle_count
        self.type= type

    # Optional attributes with default values
        self.tube_OD = 8/1000  # in meters
        self.tube_ID = 6/1000 # in meters
        self.length = 0.35 # in meters
        self.D_shell = 0.064 # in meters
        self.pipearea = 0.25 * self.tube_ID**2 * np.pi
        self.sigma = self.pipearea / (0.25* self.D_shell**2 * np.pi)

    def summary(self):
        print(f"Heat Exchanger Summary:")
        print(f"  Type: {self.type}")
        print(f"  Tube count: {self.tube_count}")
        print(f"  Baffle count: {self.baffle_count}")
        print(f"  Tube outer diameter: {self.tube_OD} m")
        print(f"  Tube inner diameter: {self.ID} m")
        print(f"  Tube length: {self.length} m")
        print(f"  Shell diameter: {self.D_shell} m")
