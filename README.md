# vigilant-potato

This script is a graphic solution to determine the critical condition in the saltwater upconing in a leaky confined aquifer based on analytical solution derived from Bower, et al (1999).  The critical rise and critical pumping rate of saltwater upconing are solved when equal pressure and equal pressure gradient conditions both exist.  Above the critical condition, the exerted pumping will induce an unstable freshwater-saltwater cone.

Two functions were solved using SciPy package-Modified Bessel function of the second kind of integer order zero.  Script was set up to solve for parameters sets posed by Figure 4 (Bower, et al., 1999) for comparison and validation.  Script was also set up to solve for user input parameters defined in the “Parameters.txt”:

#Distance from the top of the aquifer to the bottom of the well screen (m)
#Thickness of the aquifer from the bottom of the confining unit to the initial condition of the saltwater-freshwater interface (m)
#Distance from the top of the aquifer to the top of the screen (m)
#Radial distance from the pumping well (rw) (m)
#Vertical hydraulic conductivity of the aquifer (m/day)
#Horizontal hydraulic conductivity of the aquifer (m/day)
#Vertical hydraulic conductivity of the overlying confining unit (m/day)
#Thickness of the overlying confining unit (m)
#Transmissivity of the aquifer (m2/day)
#Summation factor (Define the accuracy of the solution, with a greater summation factor, accuracy improves and runtime increases) (dimensionless)

A graphical representation of the parameters can be found in Figure 3 (Bower, et al., 1999).  I would love to hear any feedback and comments.

Reference:
#Bower, J. W., Motz, L. H., and Durden, D. W. (1999). Analytical solution for determining the critical condition of saltwater upconing in a leaky artesian aquifer. Journal of Hydrology, 221(1-2), 43-54.
