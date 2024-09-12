"""""""""
Import Packages and Setup Python Environment
"""""""""
import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz

b_def   = 1.5
t_def   = 0.0
phi_def = 0.0
chi_def = 90

"""""""""
Define Streamline object class
"""""""""
class Streamline(object):

    # Initialise
    def __init__(self, streamlineData, header_dict, dataFile):
        self._x = streamlineData['x']
        self._y = streamlineData['y']
        self._u = streamlineData['u']
        self._dxdy = streamlineData['dx']

        try:
            self._Mach_b = header_dict['M_b']
        except:
            self._Mach_b = None

        try:
            self._b = header_dict['b']
        except:
            self._b = b_def

        try:
            self._t = header_dict['t']
        except:
            self._t = t_def

        try:
            self._phi_b = header_dict['phi_b']  # degrees
        except:
            self._phi = phi_def

        try:
            self._chi_b = header_dict['chi_b']  # pi
        except:
            self._chi_b = chi_def/180

        self._dataFile = dataFile

        self._path = cumtrapz(np.sin(self.theta),self.y,initial=0)
        self._flow_time = cumtrapz(np.sin(self.theta)/self.u,self.y,initial=0)

    # Spatial coordinates
    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def r(self):
        return np.sqrt(self._x**2+self._y**2)

    # Angles
    @property
    def phi(self):
        return np.arctan2(self._y,self._x)

    @property
    def theta(self):
        return np.arctan2(1,self._dxdy)

    @property
    def chi(self):
        return self.theta - self.phi

    # Velocity and components
    @property
    def u(self):
        return self._u

    @property
    def ux(self):
        return self._u * np.cos(self.theta)

    @property
    def uy(self):
        return self._u * np.sin(self.theta)

    @property
    def ur(self):
        return self._u * np.cos(self.chi)

    @property
    def utheta(self):
        return -self._u * np.sin(self.chi)

    # Mach number
    @property
    def Mach_b(self):
        return self._Mach_b

    # Model properties
    @property
    def b(self):
        return self._b

    @property
    def t(self):
        return self._t

    @property
    def phi_b(self):
        return self._phi_b    # degrees

    @property
    def chi_b(self):
        return self._chi_b    # pi

    @property
    def flow_time(self):
        return self._flow_time

    # Component interpolators
    def lookup_ux(self, x, y):
        lookup = interp1d(self.phi, self.ux, bounds_error=False)
        return lookup(np.arctan2(y,x))

    def lookup_uy(self, x, y):
        lookup = interp1d(self.phi, self.uy, bounds_error=False)
        return lookup(np.arctan2(y,x))

    def lookup_ur(self, phi):
        lookup = interp1d(self.phi, self.ur, bounds_error=False)
        return lookup(phi)

    def lookup_utheta(self, phi):
        lookup = interp1d(self.phi, self.utheta, bounds_error=False)
        return lookup(phi)

    def lookup_u(self, phi):
        lookup = interp1d(self.phi, self.u, bounds_error=False)
        return lookup(phi)

    def lookup_rtilde(self, phi):
        lookup = interp1d(self.phi, self.r, bounds_error=False)
        return lookup(phi)

    def lookup_rb(self, r, phi):
        return r/self.lookup_rtilde(phi)

    # Datafile
    @property
    def dataFile(self):
        return self._dataFile

def read_streamline(dataFile):
    # Read file
    try:
        streamlineData = np.genfromtxt(dataFile, names=True, skip_header=1, comments='#')
        with open(dataFile) as f:
            header_dict = {}
            header = f.readline().strip()
            print("Analytic streamline found: {}\n".format(dataFile), header)
            header_split = header.split(",")
            for i in range(len(header_split)):
                if i ==0:
                    header_split[i] = header_split[i].split(" ")[-1]
                header_split[i] = header_split[i].strip()
                header_split[i] = header_split[i].split("=")
                try:
                    header_dict[header_split[i][0]]=float(header_split[i][1])
                except ValueError:
                    if header_split[i][0]=='chi_b':
                        header_dict[header_split[i][0]]=float(header_split[i][1].split("p")[0])
                    elif  header_split[i][0]=='t':
                        header_dict[header_split[i][0]]=float(header_split[i][1].split("(")[0])
                    else:
                        continue
    except ValueError:
        # Catch old versions that don't have a header
        streamlineData = np.genfromtxt(dataFile, names=True, skip_header=0, comments='#')
        print("Analytic streamline found:")
        header_dict = None
    except OSError:
        print("File ", dataFile, " could not be found")
        return None

    # Store as streamline object
    stream = Streamline(streamlineData, header_dict, dataFile)

    return stream

def main():
    # Test script to make a basic plot of the streamline on the R-z plane
    """""""""
    Argument Setup
    """""""""
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", "-f", nargs='+', type=str, default=None, help='Provide *.dat file(s) for analytic streamline.')
    args = parser.parse_args()

    # Plot the streamlines
    plt.figure()
    for dataFile in args.file:
        stream = read_streamline(dataFile)
        if stream:
            plt.plot(stream.x, stream.y, label='$\mathcal{{M}} = {:03.3f}$'.format(stream.Mach_b))

    # Decorate and display plot
    plt.xlim([0,10])
    plt.ylim([0,10])
    plt.xlabel('$\\tilde{R}$')
    plt.ylabel('$\\tilde{z}$')
    plt.legend()
    plt.show()

if __name__=="__main__":
    main()
