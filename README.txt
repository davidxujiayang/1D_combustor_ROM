Simulation & Reduced-Order Modeling code for a quasi-1D model rocket combustor
Geometry & experiment setup see: Harvazinski, et al. "Coupling between 
hydrodynamics, acoustics, and heat release in a self-excited unstable combustor"

=================================================================================
How to use
All subroutines are called in start.m. You can choose to run selected parts 
by modifying control variables, RUN_STEADY, etc.

=================================================================================
Important functions
read_shape:			Creates geometry info for the system
Solve_Steady_State:	Obtains steady state result
FOM:				Does full-order unsteady simulation & generates 
					snapshot for ROM
ROM:				Does reduced-order simulation
                    Note: in ROM, CFL>1 can be used -> acceleration
FluxFunction:		Computes ROE flux
jacoRQ:				Computes the Jacobian of the residual using Van Leer method

=================================================================================
Xu, Jiayang		davidxu@umich.edu
5 March 2017
The University of Michigan

