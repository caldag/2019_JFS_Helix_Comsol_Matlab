# 2019_JFS_Helix_Comsol_Matlab
Contains code and data related to the publication: Caldag, H. O., &amp; Yesilyurt, S. (2019). Trajectories of magnetically-actuated  % helical swimmers in cylindricalchannels at low Reynolds numbers. Journal of  % Fluids and Structures, 90, 164-176. https://doi.org/10.1016/j.jfluidstructs.2019.06.005

-----------------------------------------------------------------------------------

DESCRIPTION OF THE CODES

The codes require COMSOL Multiphysics linked with a MATLAB installation. This is called
Livelink for Matlab and more information is available below:

https://www.comsol.com/livelink-for-matlab

This configuration allows building and running a COMSOL model directly from MATLAB.

Note that the final versions of the codes uploaded here could NOT be tested on Livelink, you may receive some errors.

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

Note that even though .dat file size can go up to 49 columns, not all of .dat files contain
that many columns. README.txt and helical_swim_sim.m contain the order of parameters 
and descriptions for all columns.

Contact hakanosmanc@gmail.com for any questions.
