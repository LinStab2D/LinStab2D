# LINSTAB2D - Stability and resolvent analysis of compressible viscous flows in MATLAB 


This package contains a MATLAB tool for analyzing 2D compressible and incompressible flows. The code was designed to be user-friendly, with high-level functions for mesh creation, imposing boundary conditions, etc, reducing the necessary time for setting up a new analysis and the opportunity for coding errors.

The code is capable of performing

* Temporal/global stability analysis,  
* Spatial/temporal-spatial stability analysis,
* Resolvent analysis.


The code uses finite difference schemes to construct sparse matrices representing the linearized compressible Navier-Stokes operator. It uses iterative methods to obtain eigenvalues and vectors, when performing stability analysis, and the leading gains and forcing/response modes, when performing resolvent analysis.

The code relies only on the basic MATLAB installation and should run in most MATLAB versions, although it has only been tested in MATLAB 2020. Octave is also partially supported (see notes at the end of the document for details).

The code can be downloaded or cloned from GitHub
> git clone https://github.com/eduardomartini/LinStab2D.git

or be downloaded from MATLAB central
> LINK TO BE INCLUDED  

All the examples should run out of the box. The following examples are available:

* Example 0 : A validation of the curvilinear formulation.
* Example 1 : Resolvent analysis of a Blausius boundary layer.
* Example 2 : Resolvent analysis of a round Jet.
* Example 3 : Temporal and spatial locally parallel analysis of a streaky supersonic boundary layer.
* Example 4 : Stability analysis of the flow around a cavity.
* Example 5 : Resolvent analysis around a paraboloid body.

These examples provide a good point of departure for the analysis of other flows. 


## Bugs, requests, and suggestions

If you identify any bug in the code or have any request/suggestion, please create an issue on the project's GitHub page.

## OCTAVE support
Most of the code should work without issues in Octave. There are two exceptions: *plotting* and *eigs* functions. Writing custom plotting functions is straightforward. The incompatibility with *eigs* is due to the Octave implementation not accepting complex-valued matrices.

To overcome this limitation, the system needs to be expanded to represent the complex vector as an extended vector, where the real and imaginary parts are represented in different degrees of freedom (DOFs). This functionality is already implemented in the *resolvent.m* file, which can serve as a template for future modifications. A customized version of _eigs.m_ can be included in the latter to improve compatibility.


# Citing the code
When using the code, please cite: 

> TO BE ADDED
