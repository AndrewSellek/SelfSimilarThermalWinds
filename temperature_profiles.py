"""""""""
TEMPERATURE PROFILES

The temperature profiles are described here.
Consists of a base class TemperatureStructure and examples of implementations
- Since a self-similar solution is being applied, the radial variation is always a power law of index t.
- Each implementation must provide the angular variation C(phi) as a function of cylindrical coordinates R and z.
- The solution is expressed relative to the base of the streamline, hence C(phi_b) must equal 1.
- Each implementation must also provide the derivative of ln(C) with respect to phi.
- Finally each implementation must define its "key".
- All keys should be listed in the dictionary temperatureKeys.
"""""""""
import numpy as np

### Base class describing the temperature/sound speed structure of the wind.
class TemperatureStructure(object):

    def __init__(self, t, phi_b):
        self._t = t
        self._phi_b = phi_b/180*np.pi

    @property
    def t(self):
        return self._t

    @property
    def key(self):
        # Must be specified in instances
        return NotImplementedError

    def C(self, R, z):
        # Must be specified in instances
        return NotImplementedError

    def dlnCdphi(self, R, z):
        # Must be specified in instances
        return NotImplementedError

    def cs_sq(self, R, z):
        r = np.sqrt(R**2 + z**2)
        return r**(-self.t) * self.C(R, z)

### Example implementations
### See Section 3.2 of Sellek et al. (in prep.) for details
class SphericalTemperature(TemperatureStructure):
    @property
    def key(self):
        return 's'

    def C(self, x, y):
        # Must be specified in instances
        return np.ones_like(x)

    def dlnCdphi(self, x, y):
        # Must be specified in instances
        return np.zeros_like(x)
class CylindricalTemperature(TemperatureStructure):
    @property
    def key(self):
        return 'c'

    def C(self, x, y):
        r = np.sqrt(x**2 + y**2)
        return (x/(r*np.cos(self._phi_b)))**(-self.t)

    def dlnCdphi(self, x, y):
        return self.t * y/x

### Key dictionary
temperatureKeys = {"s": SphericalTemperature, "c": CylindricalTemperature}

print("Imported temeprature profiles")
