import numpy as np
import re

class HydroCoeff:
    """
    The HydroCoeff class represents the hydrodynamic coefficients of a system.

    It calculates the damping and added mass based on the phase lag, hydrodynamic force, motion amplitude, and angular frequency.

    Attributes:
        _phaselag (float): The phase lag of the system.
        _hydrodynamicforce (float): The hydrodynamic force acting on the system.
        _motionamplitude (float): The amplitude of motion of the system.
        _w (float): The angular frequency of the system.
        _damping (float): The calculated damping of the system.
        _addedmass (float): The calculated added mass of the system.
        _restoring_coefficient (float): The calculated restoring coefficient of the system.

    Methods:
        __init__(self, phaselag: float, hydrodynamicforce: float, motionamplitude: float, w: float, rho: float, g: float, Awp: float):
            Initializes the HydroCoeff object with the given phase lag, hydrodynamic force, motion amplitude, angular frequency, rho, g, and Awp.
            It also calculates the damping, added mass, and restoring coefficient.

        calculate_restoring_coefficient(self):
            Calculates the restoring coefficient based on the given parameters.

        calculate_damping(self):
            Calculates the damping based on the hydrodynamic force, phase lag, motion amplitude, angular frequency, and restoring coefficient.

        calculate_added_mass(self):
            Calculates the added mass based on the hydrodynamic force, phase lag, motion amplitude, angular frequency, and restoring coefficient.

        __str__(self):
            Returns a string representation of the HydroCoeff object.

        __repr__(self):
            Returns a string representation of the HydroCoeff object that can be used to recreate the object.

        phaselag(self) -> float:
            Returns the phase lag.

        phaselag(self, value: float):
            Sets the phase lag to the given value and recalculates the damping, added mass, and restoring coefficient.

        hydrodynamicforce(self) -> float:
            Returns the hydrodynamic force.

        hydrodynamicforce(self, value: float):
            Sets the hydrodynamic force to the given value and recalculates the damping, added mass, and restoring coefficient.

        motionamplitude(self) -> float:
            Returns the motion amplitude.

        motionamplitude(self, value: float):
            Sets the motion amplitude to the given value and recalculates the damping, added mass, and restoring coefficient.

        w(self) -> float:
            Returns the angular frequency.

        w(self, value: float):
            Sets the angular frequency to the given value and recalculates the damping, added mass, and restoring coefficient.

        damping(self) -> float:
            Returns the damping.

        addedmass(self) -> float:
            Returns the added mass.

        restoring_coefficient(self) -> float:
            Returns the restoring coefficient.
    """

    def __init__(self, phaselag: float, hydrodynamicforce: float, motionamplitude: float, w: float, rho: float, g: float, Awp: float):
        """
        Initializes the HydroCoeff object with the given phase lag, hydrodynamic force, motion amplitude, angular frequency, rho, g, and Awp.
        It also calculates the damping, added mass, and restoring coefficient.
        """
        self._phaselag = phaselag
        self._hydrodynamicforce = hydrodynamicforce
        self._motionamplitude = motionamplitude
        self._w = w
        self._rho = rho
        self._g = g
        self._Awp = Awp
        self._damping = None
        self._addedmass = None
        self._restoring_coefficient = None

        self.calculate_restoring_coefficient()
        self.calculate_damping()
        self.calculate_added_mass()

    def calculate_restoring_coefficient(self):
        """
        Calculates the restoring coefficient based on the given parameters.
        """
        self._restoring_coefficient = self._rho * self._g * self._Awp

    def calculate_damping(self):
        """
        Calculates the damping based on the hydrodynamic force, phase lag, motion amplitude, angular frequency, and restoring coefficient.
        """
        self._damping = - (self._hydrodynamicforce - self._restoring_coefficient * self._motionamplitude) * np.sin(self._phaselag) / (self._motionamplitude * self._w)

    def calculate_added_mass(self):
        """
        Calculates the added mass based on the hydrodynamic force, phase lag, motion amplitude, angular frequency, and restoring coefficient.
        """
        self._addedmass = (self._hydrodynamicforce - self._restoring_coefficient * self._motionamplitude) * np.cos(self._phaselag) / (self._motionamplitude * self._w ** 2)

    def __str__(self):
        """
        Returns a string representation of the HydroCoeff object.
        """
        return f"HydroCoeff: phaselag={self._phaselag}, hydrodynamicforce={self._hydrodynamicforce}, motionamplitude={self._motionamplitude}, w={self._w}, damping={self._damping}, addedmass={self._addedmass}, restoring_coefficient={self._restoring_coefficient}"

    def __repr__(self):
        """
        Returns a string representation of the HydroCoeff object that can be used to recreate the object.
        """
        return f"HydroCoeff(phaselag={self._phaselag}, hydrodynamicforce={self._hydrodynamicforce}, motionamplitude={self._motionamplitude}, w={self._w}, restoring_coefficient={self._restoring_coefficient})"

    @property
    def phaselag(self) -> float:
        """
        Returns the phase lag.
        """
        return self._phaselag

    @phaselag.setter
    def phaselag(self, value: float):
        """
        Sets the phase lag to the given value and recalculates the damping, added mass, and restoring coefficient.
        """
        self._phaselag = value
        self.calculate_restoring_coefficient()
        self.calculate_damping()
        self.calculate_added_mass()

    @property
    def hydrodynamicforce(self) -> float:
        """
        Returns the hydrodynamic force.
        """
        return self._hydrodynamicforce

    @hydrodynamicforce.setter
    def hydrodynamicforce(self, value: float):
        """
        Sets the hydrodynamic force to the given value and recalculates the damping, added mass, and restoring coefficient.
        """
        self._hydrodynamicforce = value
        self.calculate_restoring_coefficient()
        self.calculate_damping()
        self.calculate_added_mass()

    @property
    def motionamplitude(self) -> float:
        """
        Returns the motion amplitude.
        """
        return self._motionamplitude

    @motionamplitude.setter
    def motionamplitude(self, value: float):
        """
        Sets the motion amplitude to the given value and recalculates the damping, added mass, and restoring coefficient.
        """
        self._motionamplitude = value
        self.calculate_restoring_coefficient()
        self.calculate_damping()
        self.calculate_added_mass()

    @property
    def w(self) -> float:
        """
        Returns the angular frequency.
        """
        return self._w

    @w.setter
    def w(self, value: float):
        """
        Sets the angular frequency to the given value and recalculates the damping, added mass, and restoring coefficient.
        """
        self._w = value
        self.calculate_restoring_coefficient()
        self.calculate_damping()
        self.calculate_added_mass()

    @property
    def damping(self) -> float:
        """
        Returns the damping.
        """
        return self._damping

    @property
    def addedmass(self) -> float:
        """
        Returns the added mass.
        """
        return self._addedmass

    @property
    def restoring_coefficient(self) -> float:
        """
        Returns the restoring coefficient.
        """
        return self._restoring_coefficient
    

def process_line(line):
    # Extract time value
    time = float(line.split()[0])

    # Extract the numeric values within parentheses
    float_values = []
    tokens = re.findall(r'\(([^)]+)\)', line)
    for token in tokens:
        # Removing any lingering parentheses and splitting by spaces
        cleaned_values = token.replace('(', '').replace(')', '').split()
        float_values.extend([float(x) for x in cleaned_values])

    # Return time followed by the extracted float values
    return [time] + float_values


def createForceFile(forces_file: str) -> tuple:
    data = []

    with open(forces_file, "r") as datafile:
        for line in datafile:
            if not line.startswith("#"):
                data.append(tuple(process_line(line)))  # Convert to tuple

    # Define the dtype for the structured array
    dtype = [('time', float), 
             ('pressure_x', float), ('pressure_y', float), ('pressure_z', float),
             ('viscous_x', float), ('viscous_y', float), ('viscous_z', float),
             ('pressure_moment_x', float), ('pressure_moment_y', float), ('pressure_moment_z', float),
             ('viscous_moment_x', float), ('viscous_moment_y', float), ('viscous_moment_z', float)
            ]

    # Create a structured NumPy array
    data_array = np.array(data, dtype=dtype)

    # Extract the required data
    time = data_array['time']
    forceX = data_array['pressure_x'] #+ data_array['viscous_x']
    forceY = data_array['pressure_y'] #+ data_array['viscous_y']
    forceZ = data_array['pressure_z'] + data_array['viscous_z']

    return time, forceX, forceY, forceZ


