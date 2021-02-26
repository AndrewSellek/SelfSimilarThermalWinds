# Self Similar Thermal Winds
Author: Andrew Sellek
Date: February 2020

This package is designed to calculate the structure of self-similar thermal winds launched from discs as explored in Sellek, Clarke and Booth (in prep.).
Such models have been shown to provide robust predictions for the structure of winds, particularly in context describing winds from protoplanetary discs.

The package contains tools for three purposes:
* To calculate the maximum Mach number with which the wind can launch.
* To calculate the Cartesian equation of a streamline, along with its velocity structure.
* To read the outputs of the above.

The scenarios are described, primarily, by four parameters:
* The gradient b of the density profile at the base, assumed to be a power law \\( \rho \propto r^{-b} \\)
* The gradient t of the imposed temperature profile, assumed to be a power law \\( T \propto r^{-t} \\)
* The inclination of the wind base to the midplane \\( \phi\_{\rm b} \\)
* The angle at which the wind emerges from the base \\( \chi\_{\rm b} \\)
In addition, any wind solution also depends on the Mach number at the base \\( \mathcal{M}\_{\rm b} \\)

## Analytic Solution (analytic\_solution.py)
Uses the equations in the form given in Sellek et al. (in prep.) to integrate the streamline solution using a first order Euler method.
The solver uses x to represnet the cylindrical radius R and y to represent the cylindrical coordinate z.

There are three modes of operation:
* **\[Default\]** Calculate the streamline solution for a given set of wind parameters and given Mach number(s).
* **\[Search\]** Find the maximum Mach number allowed in a given scenario to 3 decimal places.
* **\[Refine\]** Refine a previous estimate of the maximum Mach number, for example using a higher resolution or a larger integration limit than previously.

### Command Line Arguments
The wind parameters are specified using the following arguments; the code solves for all combinations:
* **\[-b\]** Density Profile Index \[Default: b=1.5\]
* **\[-t\]** Temperature Profile Index \[Default: t=0.0 (isothermal)\]
* **\[-p\]** Inclination of the wind base (degrees) \[Default: 0 degrees\]
* **\[-c\]** Angle of the wind at the base (degrees) \[Default: 90 degrees\]
* **\[-k\]** Temperature Key (see below)

The mode may be selected using:
* **\[-M\]** To specify the Mach number(s) to solve. In search mode the first value only will be used as a starting point for the search. \[Default: 0.1\]
* **\[-S\]** To use Search
* **\[-R\]** To use Refine
* **\[-f\]** To specify the file from which to read/write maximum Mach numbers \[Default: 'launch_Mach.dat'\]

The integration of the solution depends on:
* **\[-y\]** The y value at which to stop the integration \[Default: 100\]
* **\[-n\]** The step in y to use for the integration \[Default: 10^-5\] (N.B. it is recommended to use a resolution smaller than dy = 10^-2 \* Mb^2, the code will by default enforce this limit)

The final arguments are:
* **\[-N\]** The number of separate parallel processes to run \[Default: 4\]
* **\[-v\]** Verbose 

### Temperature Structures (temperature\_profiles.py)
Although the code defaults to a spherically symmetric temperature profile (key='s'), the user may choose instead a cylindrical profile (key='c'):
    python analytic_solution.py -k c

The temperature stucture is treated using classes, so a user may define their own temperature class with a custom angular dependence in the file temperature_profiles.py

Each implementation must provide the angular variation C(phi) as a function of cylindrical coordinates R and z. The solution is expressed relative to the base of the streamline, hence C(phi_b) must equal 1. Each implementation must also provide the derivative of ln(C) with respect to phi.

Finally each implementation must define its "key". All keys must be listed in the dictionary temperatureKeys in order to be used.

### Example
Winds with b=1.0, an inclination of 36 degrees and t=0.5 using both spherical and cylindrical geometry for the temperature.
Solve for the maximum Mach number by integrating to y=1000 and store in my\_results.dat
    python analytic\_solution.py -b 1.0 -p 36 -t 0.5 -k s c -S -y 1000 -f my\_results.dat

## Table Reader (table\_read.py) 
Reads the tables of maximum Mach numbers produced by analytic_solution.py

    read\_Mach\_table(tableFile)
Returns a dictionary of values read from a file 'tableFile' containing the maximum Mach number of a range of scenarios.
    search\_Mach\_table(Mach\_table, b, t, phi\_b, chi\_b)
Searches the table for the maximum Mach number corresponding to the parameters b, t, phi\_b, chi\_b

A sample table of values, containing all those used to make make Figure 1 of Sellek et al. (in prep.) is given in examples/launch\_Mach.dat

## Streamline Reader (streamline\_read.py)
Reads the streamline_data\*.dat files produced by analytic\_solution.py

    read_streamline(dataFile)
Returns a streamline object. Its attributes include the wind parameters for which it was run, the velocity and its components, an methods to return the radius/velocity at requested values of the elevation phi.

Sample streamlines are given in the folder examples/streamlines
