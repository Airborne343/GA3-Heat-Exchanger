from dataclasses import dataclass, field

class HeatExchanger:
    # Required attributes
    tube_count: int
    baffle_count: int
    type: str  # e.g., 'shell-and-tube', 'plate', etc.

    # Optional attributes with default values
    tube_OD: float = field(default=8/1000)  # in meters
    tube_ID: float = field(default=6/1000)  # in meters
    length: float = field(default=0.35) # in meters
    D_shell: float = field(default=0.064) # in meters

    def summary(self):
        print(f"Heat Exchanger Summary:")
        print(f"  Type: {self.type}")
        print(f"  Tube count: {self.tube_count}")
        print(f"  Baffle count: {self.baffle_count}")
        print(f"  Tube outer diameter: {self.tube_OD} m")
        print(f"  Tube inner diameter: {self.ID} m")
        print(f"  Tube length: {self.length} m")
        print(f"  Shell diameter: {self.D_shell} m")

     def __post_init__(self):
        # Compute derived attribute
        