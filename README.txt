This directory contains data and code related to the publication:

Caldag, H. O., & Yesilyurt, S. (2019). Trajectories of magnetically-actuated 
helical swimmers in cylindricalchannels at low Reynolds numbers. Journal of 
Fluids and Structures, 90, 164-176. https://doi.org/10.1016/j.jfluidstructs.2019.06.005

-----------------------------------------------------------------------------------

DESCRIPTION OF THE CODES

The codes require COMSOL Multiphysics linked with a MATLAB installation. This is called
Livelink for Matlab and more information is available below:

https://www.comsol.com/livelink-for-matlab

This configuration allows building and running a COMSOL model directly from MATLAB.

Brief descriptions for the scripts are given below:

- helical_swim_sim.m is the main script to be run. It sets up the model, runs FEM simulations
and kinematically integrates the velocities to obtain trajectories.

- initial_build.m prepares the base COMSOL-based finite-element method (FEM) model to be used in
simulations.

- sim_run.m is called in-between the time steps. It updates the base model with parameters for the
current time step and solves the system, returning velocity values to helical_swim_sim.m to be 
kinematically integrated for setting up the next time step.

-----------------------------------------------------------------------------------

EXPLANATION ON DATA

Folders contain data relevant to the plots in the article. The files are named as combinations of letters
and numbers which contain information about the physical configuration.

Not all filenames have all the identifiers listed below. In that case, the value used should be the default one,
the value in helical_swim_sim.m.

D: Channel diameter
L: Number of wavelengths of the helical tail
Q: Flowrate
f: Magnetic field rotation frequency (Negative f indicates puller-mode, positive f indicates pusher-mode)
z0: Initial position in the z-direction
B: Magnetic field strength (indicator for Mason number, default is 3600)
dt: Time step size
thax: Misalignment angle (in degrees)
NG: Indicates gravity is not included in the computations (No Gravity)

-----------------------------------------------------------------------------------

DATA STRUCTURE

Folder D1.6-L4 contains refcode_fordatainfo.m file which includes the detailed information on what each column
on .dat files represent. Note that even though .dat file size can go up to 49 columns, not all of .dat files contain
that many columns.

DESCRIPTION OF THE VARIABLES IN .dat FILES

The descriptions are also available at the end of helical_swim_sim.m
Each row in data corresponds to a time instant simulated.

Contact hakanosmanc@gmail.com for any questions.

% Number	|	Name		|	Description
% ==========================================================================================
% 1		| 	allt		|	Vector containing all simulated time steps
% 2		|	U		|	Linear velocity in the x-direction
% 3		|	V		|	Linear velocity in the y-direction
% 4		|	W		|	Linear velocity in the z-direction
% 5		|	avx		|	Angular velocity in the x-direction
% 6		|	avy		|	Angular velocity in the y-direction
% 7		|	avz		|	Angular velocity in the z-direction
% 8		|	rb(1)		|	x-coordinate of the swimmer
% 9		|	rb(2)		|	y-coordinate of the swimmer
% 10		|	rb(3)		|	z-coordinate of the swimmer
% 11		|	Fvx		|	Viscous force, x-component
% 12		|	Fvy		|	Viscous force, y-component
% 13		|	Fvz		|	Viscous force, z-component
% 14		|	Ftz		|	Tail weight (buoyancy force subtracted)
% 15 		|	Fhz		|	Head weight (buoyancy force subtracted)
% 16		|	Tvzx		|	Viscous torque in the x-direction due to the forces in the z-direction
% 17		|	Tvyx		|	Viscous torque in the x-direction due to the forces in the y-direction
% 18		|	Tmx		|	Magnetic torque in the x-direction
% 19		|	Tgx		|	Gravity-induced torque in the x-direction
% 20		|	Tvxy		|	Viscous torque in the y-direction due to the forces in the x-direction
% 21		|	Tvzy		|	Viscous torque in the y-direction due to the forces in the z-direction
% 22		|	Tmy		|	Magnetic torque in the y-direction	
% 23		|	Tgy		|	Gravity-induced torque in the y-direction
% 24		|	Tvyz		|	Viscous torque in the z-direction due to the forces in the y-direction
% 25		|	Tvxz		|	Viscous torque in the z-direction due to the forces in the x-direction
% 26		|	Tmz		|	Magnetic torque in the z-direction
% 27		|	Tgz		|	Gravity-induced torque in the z-direction
% 28		|	Fvhx		|	Viscous force in the x-direction (for head only)
% 29		|	Fvhy		|	Viscous force in the y-direction (for head only)
% 30		|	Fvhz		|	Viscous force in the z-direction (for head only)
% 31		|	Fvtx		|	Viscous force in the x-direction (for tail only)
% 32		|	Fvty		|	Viscous force in the y-direction (for tail only)
% 33		|	Fvtz		|	Viscous force in the z-direction (for tail only)
% 34		|	Tvhzx		|	Viscous torque in the x-direction due to the forces in the z-direction (head only)
% 35		|	Tvhyx		|	Viscous torque in the x-direction due to the forces in the y-direction (head only)
% 36		|	Tvhxy		|	Viscous torque in the y-direction due to the forces in the x-direction (head only)
% 37		|	Tvhzy		|	Viscous torque in the y-direction due to the forces in the z-direction (head only)
% 38		|	Tvhyz		|	Viscous torque in the z-direction due to the forces in the y-direction (head only)
% 39		|	Tvhxz		|	Viscous torque in the z-direction due to the forces in the x-direction (head only)
% 40		|	Tvtzx		|	Viscous torque in the x-direction due to the forces in the z-direction (tail only)
% 41		|	Tvtyx		|	Viscous torque in the x-direction due to the forces in the y-direction (tail only)
% 42		|	Tvtxy		|	Viscous torque in the y-direction due to the forces in the x-direction (tail only)
% 43		|	Tvtzy		|	Viscous torque in the y-direction due to the forces in the z-direction (tail only)
% 44		|	Tvtyz		|	Viscous torque in the z-direction due to the forces in the y-direction (tail only)
% 45		|	Tvtxz		|	Viscous torque in the z-direction due to the forces in the x-direction (tail only)
% 46		|	wcall		|	Checks if there is wall contact (whole geometry)
% 47		|	wchead		|	Checks if there is wall contact (head only)
% 48		|	wctail		|	Checks if there is wall contact (tail only)
% 49		|	Fwr		|	Repulsive force (in case of contact) in the radial direction
