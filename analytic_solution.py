"""""""""
Import Packages and Setup Python Environment
"""""""""
import numpy as np
import argparse
import os
from scipy.interpolate import interp1d
from scipy.special import erfc, erfcinv
import multiprocessing
import time
import itertools

"""""""""
Custom Tools
"""""""""
from table_read import *

"""""""""
Scenario Classes
"""""""""
from temperature_profiles import *

class WindBase(object):
    # Describes the base of the wind and its geometry
    def __init__(self, b, phi_b, chi_b):
        self._b = b
        self._phi_b = phi_b             # degrees
        self._chi_b = chi_b             # degrees

    @property
    def b(self):
        return self._b

    @property
    def phi_b(self):
        return self._phi_b/180*np.pi    # radians

    @property
    def phi_deg(self):
        return self._phi_b              # degrees

    @property
    def chi_b(self):
        return self._chi_b/180*np.pi    # radians

    @property
    def chi_deg(self):
        return self._chi_b              # degrees

    @property
    def theta_b(self):
        return self.chi_b + self.phi_b  # radians

"""""""""
Main Routines
"""""""""

def main():
    """""""""
    Argument Setup
    """""""""
    parser = argparse.ArgumentParser()

    # Density profile, temperature and base angles
    parser.add_argument("--b", "-b", type=float, nargs='+', default=[1.5], help='Density Profile Index')
    parser.add_argument("--t", "-t", type=float, nargs='+', default=[0.0], help='Temperature Profile Index')
    parser.add_argument("--k", "-k", type=str, nargs='+', default=['s'], choices=temperatureKeys.keys(), help='Temperature Key: s=spherical, c=cylindrical, {custom}')
    parser.add_argument("--phib", "-p", type=float, nargs='+', default=[0.0],  help='Inclination of the wind base (degrees) [Default: 0]')
    parser.add_argument("--chib", "-c", type=float, nargs='+', default=[90.0], help='Angle of the wind at the base (degrees) [Default: 90]')

    # Initial value and outer radius of solution
    parser.add_argument("--yend", "-y", type=float, default=100, help='Largest y to integrate out [100]')
    parser.add_argument("--resolution", "-n", type=float, default=1e-5, help='dy to use for integration [1e-5]')
    parser.add_argument("--tryMach", "-M", type=float, nargs='+', default=[0.1], help='Mach number(s) to try')

    # Search and refine modes
    parser.add_argument("--search", "-S", action='store_true', help='First find the maximum.')
    parser.add_argument("--refine", "-R", action='store_true', help='Refine search from a previous successful value.')
    parser.add_argument("--file", "-f", type=str, default='launch_Mach.dat', help='File to store results.')

    # Run parameters
    parser.add_argument("--noProcesses", "-N", type=int, default=4, help='Maximum number of multiprocessing processes to run at once.')
    parser.add_argument("--verbose", "-v", action='store_true', help='Print more updates such as Mach number being tried.')

    # Interpret
    global verbose
    args     = parser.parse_args()
    search   = args.search
    refine   = args.refine
    saveFile = args.file
    verbose  = args.verbose

    # Check for non-default modes
    if search or refine:
        """Get ready to save"""
        # Check for savefile and open, initializing with header if need be
        results_exist = os.path.exists(saveFile)
        filehandle = open(saveFile, 'a')
        if results_exist==False:
            print("Creating output file")
            filehandle.write("b\tT_key\tt\tphi_b\tchi_b/pi\tM_b\tdy\ty_end\n")
        filehandle.close()
    if refine:
        if len(args.k)>1:
            raise Exception("Table read cannot deal with temperature key range")
        launch_Mach = read_Mach_table(saveFile, choose_k=args.k[0])
    else:
        launch_Mach = None

    # Run pool of solvers for each combination of input parameters
    pool = multiprocessing.Pool(processes=min(args.noProcesses, len(args.b)*len(args.t)*len(args.phib)*len(args.chib)*len(args.k)))
    results = [pool.apply_async(solve_streamline, args=(b, p, c, t, k, args, search, refine, launch_Mach,)) for b, t, p, c, k in itertools.product(args.b, args.t, args.phib, args.chib, args.k)]

    # Deal with results
    if search or refine:
        for p in results:
            result = p.get()
            if result==0:
                print("Temperature structure error")
            else:
                print(time.asctime())
                print("Saving: ", result) 
                filehandle = open(saveFile, 'a')
                filehandle.write(result)
                filehandle.close()

    # Finish
    pool.close()
    pool.join()

def search_up(base, cs, dy, yend, estMach, dM, adaptResolution=True):
    # Estimate Mach number to precision dM
    negativeWarn = False
    if estMach==0:
        if verbose:
            print("M=0 not valid: setting M=dM={}".format(dM))
        estMach = dM

    # Increase Mach number until f<0 encountered
    while not negativeWarn:
        if verbose:
            print("Test {:3.3}".format(estMach))
        if dy > 1e-2 * estMach**2 and adaptResolution:
            dyAdapt = 10**(-2 + np.floor(2*np.log10(estMach)))
            print("Higher resolution needed, increasing to dy={} for M={}".format(dyAdapt, estMach))
            negativeWarn = Mach_iteration(base, cs, dyAdapt, yend, estMach)
        else:
            negativeWarn = Mach_iteration(base, cs, dy, yend, estMach)
        estMach += dM

    # Return Mach number to highest valid value
    estMach -= 2*dM

    return estMach

def search_down(base, cs, dy, yend, estMach, dM, adaptResolution=True):
    # Refine Mach number to precision dM
    negativeWarn = True

    # Decrease Mach number until f<0 no longer encountered
    while negativeWarn:
        if verbose:
            print("Test {:3.3}".format(estMach))
        if dy > 1e-2 * estMach**2 and adaptResolution:
            dyAdapt = 10**(-2 + np.floor(2*np.log10(estMach)))
            print("Higher resolution needed, increasing to dy={} for M={}".format(dyAdapt, estMach))
            negativeWarn = Mach_iteration(base, cs, dyAdapt, yend, estMach)
        else:
            negativeWarn = Mach_iteration(base, cs, dy, yend, estMach)
        estMach += dM

    # Return Mach number to highest valid value
    estMach -= dM

    return estMach

def solve_streamline(b, p, c, t, k, args, search, refine, launch_Mach=None):
    base = WindBase(b, p, c)
    cs = temperatureKeys[k](t, p)

    if search:
        # First search for M_b,max 
        print("Search for", base.b, cs.t, base.phi_b, base.chi_b)
        initMach = args.tryMach[0]
        maxMach = initMach
        for dM in [0.1, 0.01, 0.001]:
            maxMach = search_up(base, cs, args.resolution, args.yend, maxMach, dM)
        calcMach = [maxMach]

    elif refine:
        # First refine estimate for M_b,max 
        print("Refine for", base.b, cs.t, base.phi_b, base.chi_b)
        initMach = launch_Mach[base.b][cs.t][base.phi_deg][float('{:.3f}'.format(base.chi_b/np.pi))]
        maxMach = search_down(base, cs, args.resolution, args.yend, initMach, -0.001)
        calcMach = [maxMach]

    else:
        # Calculate only specified cases
        print("Calculate for b={}, t={}, phi_b={}, chi_b={}".format(base.b, cs.t, base.phi_deg, base.chi_b))
        calcMach = args.tryMach

    for M in calcMach:
	# Solve and save desired streamline(s)
        if M>0:
            print("Mb={}...".format(M))
            Mach_iteration(base, cs, args.resolution, args.yend, M, save=True)
        else:
            print("Mb<=0 is not valid. Skipping...")

    if search or refine:
        return "{:3.2f}\t{}\t{:3.2f}\t{:2.0f}\t{:3.3f}\t{:3.3f}\t{}\t{}\n".format(base.b, cs.key, cs.t, base.phi_deg, base.chi_b/np.pi, maxMach, args.resolution, args.yend)
    else:
        return 1

"""""""""
Functions for updating u
Equations 24-28 of Sellek et al. (in prep.)
"""""""""
def f(x, y, dx, u, cs, Mach):
    f1 = -Mach**2 * u * (Mach**2*u**2/cs.cs_sq(x, y)-1) * (x - y*dx) / ((1+dx**2)**(0.5) * (x*dx + y))
    f2 =  Mach**2 * u * (x*dx + y) / ((1+dx**2)**(0.5) * (x - y*dx))
    f = f1 + f2
    #try:
    if (f <= 0):
        print("M = {} breaks when (x,y)=({},{})".format(Mach, x, y))
        return f, True
    #except ValueError:
    #    pass
    return f, False
def g(x, y, dx, u, base, cs, Mach):
    g1 = (base.b + cs.t) * (1 + dx**2)**(0.5) / (x - y*dx) * cs.cs_sq(x, y)
    g21 = - Mach**2 * u**2 * dx * (x - y*dx) / ((1+dx**2)**(0.5) * x * (x*dx+y))
    g22 = - Mach**2 * u**2 * (x - y*dx) / (1+dx**2)**(0.5) / (x**2+y**2) * ( cs.t - cs.dlnCdphi(x,y) * (x - y*dx) / (x*dx + y) )
    return g1 + g21 + g22
def dudy(x, y, dx, u, base, cs, Mach):
    G = g(x, y, dx, u, base, cs, Mach)
    F, negativeWarn = f(x, y, dx, u, cs, Mach)
    return G/F, negativeWarn

"""""""""
Functions for updating x
"""""""""
def A(x, y, dx, base):
    # Equation 15
    return x/np.cos(base.phi_b) * (x - y*dx) / (1+dx**2)**(0.5)/np.sin(base.chi_b)
def dAdy(u, du, x, y, dx, base, cs, Mach):
    # Equation 23
    Area = A(x, y, dx, base)
    invrho = Area*u
    return invrho * ( (Mach**2/cs.cs_sq(x, y) - 1/u**2) * du + 1/u * ( -cs.t + cs.dlnCdphi(x,y) * (x - y*dx) / (x*dx + y) ) * (x*dx+y) / (x**2+y**2) )
def d2xdy2(x, y, dx, dA, base):
    # Equation 18
    t1 = (1+dx**2) * (x - y*dx) * dx
    t2 = (1+dx**2)**(1.5) * dA  * np.cos(base.phi_b) * np.sin(base.chi_b)
    div = x * (x*dx + y)
    return (t1 - t2) / div

"""""""""
Solution for a given Mach number
"""""""""
def Mach_iteration(base, cs, dy, yend, Mach, save=False, save_dy=1e-2):
    negativeWarn=False

    """Use y (cylindrical z) as independent variable"""
    y0 = np.sin(base.phi_b)
    y  = np.arange(y0, yend+dy, dy)

    """Use x (cylindrical radius), dx and u as dependent variables and set to initial values"""
    x   = np.cos(base.phi_b)
    dx  = np.cos(base.theta_b)/np.sin(base.theta_b)
    u   = 1
    if save:
        yarr  = [y0]
        xarr  = [x]
        dxarr = [dx]
        uarr  = [u]
    
    """Loop through y"""
    for N in range(1,len(y)):
        yi = y[N]

        if N==1 and base.chi_b == np.pi/2:
            """If launched normal, evaluate using known limits at the base to avoid large numbers"""
            #print("Using known limits at the base")
            dA  = dx / np.cos(base.phi_b)
            du  = dA / (Mach**2-1)
            d2x = (base.b+cs.t) / (Mach**2 * np.cos(base.phi_b)**3)                

        """Evaluate gradents in velocity, area and dx"""
        du, negativeWarn = dudy(x, yi, dx, u, base, cs, Mach)   # derivative of u with respect to y
        dA  = dAdy(u, du, x, yi, dx, base, cs, Mach)            # derivative of A with respect to y
        d2x = d2xdy2(x, yi, dx, dA, base)                       # derivative of dx with respect to y

        """Advance u"""
        u  += du*dy

        """Advance x and x'"""
        x  += dx*dy+0.5*d2x*dy**2
        dx += d2x*dy

        """Add to arrays"""
        if save and (N % int(save_dy/dy) == 0):
            xarr.append(x)
            dxarr.append(dx)
            uarr.append(u)
            yarr.append(yi)
        elif save and negativeWarn:
            xarr.append(x)
            dxarr.append(dx)
            uarr.append(u)
            yarr.append(yi)

        if negativeWarn:
            break

    if save:
        """Save streamline"""
        output_data = np.array([xarr,yarr,dxarr,uarr])
        datafile = 'streamline_data_b{:03.0f}_t{}{:03.0f}_p{:02.0f}_c{:03.0f}_M{:04.0f}.dat'.format(100*base.b,cs.key,100*cs.t,base.phi_deg,100*base.chi_b/np.pi,Mach*1000)
        print("Saving", datafile)
        np.savetxt(datafile, output_data.T, delimiter='\t', header='Generated for b={}, t={}({}), phi_b={}, chi_b={}pi, M_b={}\nx\ty\tdx\tu'.format(base.b,cs.t,cs.key,base.phi_deg,base.chi_b/np.pi,Mach))

    return negativeWarn

if __name__ == "__main__":
    main()
