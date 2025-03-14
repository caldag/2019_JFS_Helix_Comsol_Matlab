% Part of the code package related to the publication:
%
% Caldag, H. O., & Yesilyurt, S. (2019). Trajectories of magnetically-actuated 
% helical swimmers in cylindricalchannels at low Reynolds numbers. Journal of 
% Fluids and Structures, 90, 164-176. https://doi.org/10.1016/j.jfluidstructs.2019.06.005
%
% This script is to be used with COMSOL with MATLAB software that comes with
% COMSOL Multiphysics. The program opens up a MATLAB interface that can 
% understand COMSOL-specific commands.
%
% This script makes initializations and the kinematic time-stepping between
% steady-state finite-element solutions of swimming of a helical robot in
% low Reynolds number environment.
%
% The code goes through the time steps and calls the function helical_sim. The
% inputs sent are position and orientation of the robot at time t, alongside the
% physical and geometric parameters of the system.
%
% helical_sim builds and runs the COMSOL-based finite-element model and returns
% the instantaneous linear and angular velocities. Then the instantaneous velocities are
% used to advance the swimmer position and orientation to time t+dt.

%% Geometric and physical setup
clearvars; clc; close all;
lam = 1; % Helix wavelength
nlam = 4; % Number of rotations of the helix
Lh = 1.5; % Robot head length (cylindrical head)
Db = 0.8; % Robot head diameter (cylindrical head)
rt  = 0.1; % Minor radius of the tail
B = Db/2-2*rt; % Major radius of the tail
voltail = pi*rt*rt*sqrt(4*pi*pi*B*B+lam*lam)*nlam; % Tail volume
volhead= pi*(Db/2)^2*Lh; % Head volume
rhohead = 7; % Head density (denser than tail because of the magnet)
rhotail= 1; % Tail density
x0c = (rhotail*voltail*(nlam/2*lam)-rhohead*volhead*Lh/2)/(rhotail*voltail+rhohead*volhead);
% x0c is the center-of-mass of the whole swimmer

Rch =1.6; % Channel radius
Lch = 10; % Channel length 

grav = 9.81; % Gravitational acceleration (m^2/s)
freq = 15; % Magnetic field rotation frequency (Hz)
omega = 2*pi*sign(freq); % Angular frequency (non-dimensional)
Q = 0; % Channel flow rate
Vin = Q*0.0207/abs(freq); % Corresponding inlet velocity

dt = 1/40; % Time step
t0 = 0; % Initial time
tF = 80; % Final time
allt = [t0:dt:tF]; % Time array

SAVE = 1; % Boolean toggle to save output or not
COMSOL_SAVE = 0; % Boolean toggle to save individual COMSOL simulations for each time step (Costly and takes up space)
COMSOL_PLOT = 0; % Boolean toggle to do plotting
GRAV = 1;   % Boolean toggle to include or exclude gravity
WALL = 1;   % Boolean toggle to include or exclude wall reaction force

run = strcat('TD',num2str(2*Rch),'_L',num2str(nlam*lam),'_Q',num2str(Q),'_f',num2str(freq)); % Filename that retains main system parameters
if ~isdir(run); mkdir(run); end % All relevant data saved in a directory with this name

K = max(size(allt));
U = zeros(size(allt)); % Stores the velocity values in the x-direction
V = zeros(size(allt)); % Stores the velocity values in the y-direction
W = zeros(size(allt)); % Stores the velocity values in the z-direction
avx = zeros(size(allt)); % Stores the angular velocity values in the x-direction
avy = zeros(size(allt)); % Stores the angular velocity values in the x-direction
avz = zeros(size(allt)); % Stores the angular velocity values in the x-direction
angx = zeros(size(allt)); % Misalignment angle (with respect to the x- axis)
theta = zeros(size(allt)); % Tangential component of position in the cylindrical frame
rball = zeros(3,K); % % Stores the x, y, and z coordinates of the center-of-mass of the swimmer

Fvx=U;Fvy=U;Fvz=U;Ftz=U;Fhz=U;Tvzx=U;Tvyx=U;Tmx=U;Tgx=U;
Tvxy=U;Tvzy=U;Tmy=U;Tgy=U;Tvyz=U;Tvxz=U;Tmz=U;Tgz=U;Fvhx=U;
Fvhy=U;Fvhz=U;Fvtx=U;Fvty=U;
Fvtz=U;Tvhzx=U;Tvhyx=U;Tvhxy=U;Tvhzy=U;Tvhyz=U;Tvhxz=U;
Tvtzx=U;Tvtyx=U;Tvtxy=U;Tvtzy=U;Tvtyz=U;Tvtxz=U;
wcall=U;wchead=U;wctail=U;Fwr=U;

u = zeros(6,1); % Instantaneous velocity vector containing the linear and angular velocities
U0 = 0; V0 = 0; W0 = 0; % Initial velocities in the x-, y- and z-directions.
vel0 = [U0;V0;W0]; % Initial velocity vector
vel1= vel0; vel2 = vel0; % These will be velocities at earlier time steps, we only initialize them
R=zeros(3); % Rotation matrix from the global frame to the swimmer frame
R0 = eye(3); % Initial rotation matrix set as identity matrix

Rst=zeros(3,3,K);
wxR1 = zeros(3); wxR2 = zeros(3); % Cross products of angular velocity vector with the rotation matrix at earlier time steps
% will be used to advance the orientation of the swimmer
rb0 = [0;0;0]; % Initial position
ang0 = 0; % Rotation angle of the robot around the main propulsion axis (the x- axis)
ma0 = [0;1;0]; % Initial magnetization vector

AB = 0; % Toggle to enable Adams-Bashforth integration
CN=1; % Toggle to enable Crank-Nicholson integration
% Set only one of them to 1

if (AB && CN) || (~AB && ~CN)
	error('Error in selecting a time-stepper! Set either AB or CN equal to 1 and the other to 0.')
end

model=initialbuild(nlam,lam,Db,rt,Lh,Rch,Lch,rb0,R,ma0,omega,0,grav,...
         Vin,rhohead,rhotail,abs(freq),GRAV,WALL);
		 % This is the first time the finite-element model is built in COMSOL
		 % The first simulation is outside of the script to fit with the loop structure below.
k=0;
for k = k1:K % In-between the stated time steps
    pst = allt(k); % The time instance simulated
    wxR0 = [cross(u(4:6),R0(:,1))';cross(u(4:6),R0(:,2))';cross(u(4:6),R0(:,3))']';
    vel0 = u(1:3);
	
	% Below is the Adams-Bashforth integration to advance position and orientation
	% Notice that simpler formulations are used at the initial time steps
     if k==1 || ~AB
         RR = R0  + dt*wxR0;          % Rotation matrix
         rb = rb0 + vel0*dt;          % position
         ang = ang0 + u(4)*dt*180/pi; % angular rotation of the helix
    end
    if AB && k==2
        RR = R0 + dt*(3/2*wxR0-1/2*wxR1);               % Rotation matrix
        rb = rb0 + dt*(3/2*vel0-1/2*vel1);              % position
        ang = ang0 + dt*(3/2*u(4)-1/2*avx(k-1))*180/pi; % angular rotation of the helix
    end
    if AB && k>=3
        RR = R0 + dt*((23/12)*wxR0-(4/3)*wxR1+(5/12)*wxR2);             % Rotation matrix
        rb = rb0 + dt*(23/12*vel0-4/3*vel1+5/12*vel2);                  % position
        ang = ang0 + dt*(23/12*u(4)-4/3*avx(k-1)+5/12*avx(k-2))*180/pi; % helix rotation
    end
    
	% Below is the formulation based on Crank-Nicholson
    if CN
       omegav=u(4:6);
       OHM=[0 -omegav(3) omegav(2); omegav(3) 0 -omegav(1); -omegav(2) omegav(1) 0];
       RN1=inv(eye(3)-dt/2*OHM);
       RR(:,1)=RN1*(eye(3)+dt/2*OHM)*R0(:,1);
       RR(:,2)=RN1*(eye(3)+dt/2*OHM)*R0(:,2);
       RR(:,3)=RN1*(eye(3)+dt/2*OHM)*R0(:,3);

            if k==1
                rb = rb0 + vel0*dt;          % position
                ang = ang0 + u(4)*dt*180/pi;
            elseif k==2
                rb = rb0 + dt*(3/2*vel0-1/2*vel1);              % position
                ang = ang0 + dt*(3/2*u(4)-1/2*avx(k-1))*180/pi; % angular rotation of the helix
            else
                rb = rb0 + dt*(23/12*vel0-4/3*vel1+5/12*vel2);    % position
                ang = ang0 + dt*(23/12*u(4)-4/3*avx(k-1)+5/12*avx(k-2))*180/pi; % helix rotation
            end
    end
    
    [RL,RS,RV] = svd(RR); R = RL*RV'; % Normalized rotation matrix
    Rst(:,:,k)=R;
    ha = R*[1;0;0];
    ma = R*ma0; % Rotation of the magnetization vector
    ja = R*[0;1;0];
    rot = rem(ang,360); % helix rotation in degrees
	
	% beta, alpha and gamma are orientation angles
    beta = acos(R(1,1))*180/pi;
    alpha =  atan2(-R(3,1),R(2,1))*180/pi;
    gamma = atan2(R(1,3),R(1,2))*180/pi;

     [u,model] = sim_run(model,nlam,lam,Db,rt,Lh,Rch,Lch,rb,R,ma,omega,pst,grav,...
         Vin,rhohead,rhotail,abs(freq),GRAV,WALL,COMSOL_SAVE,COMSOL_PLOT,run);
		 % Running the next time step with new position and orientation information
		 
    wxR2 = wxR1; wxR1 = wxR0; % Move old data back one step to variables keeping older data
    vel2 = vel1; vel1 = vel0; vel0 = u(1:3);
    R0 = R; % Set the new "initial" rotation matrix, position and helix rotation
    rb0 = rb;
    ang0 = ang;

    U(k) = u(1);
    V(k) = u(2);
    W(k) = u(3);
    avx(k) = u(4);
    avy(k) = u(5);
    avz(k) = u(6); % Assign velocities solved into the vectors
    rball(:,k) = rb;
    theta(k) = atan2(rb(3),rb(2));
    angx(k) = ang; % Assign position and orientation into the vectors
    
    Fvx(k)=u(7);Fvy(k)=u(8);Fvz(k)=u(9);Ftz(k)=u(10);Fhz(k)=u(11);
	Tvzx(k)=u(12);Tvyx(k)=u(13);Tmx(k)=u(14);Tgx(k)=u(15);Tvxy(k)=u(16);
	Tvzy(k)=u(17);Tmy(k)=u(18);Tgy(k)=u(19);Tvyz(k)=u(20);Tvxz(k)=u(21);Tmz(k)=u(22);
	Tgz(k)=u(23);Fvhx(k)=u(24);Fvhy(k)=u(25);Fvhz(k)=u(26);Fvtx(k)=u(27);Fvty(k)=u(28);
	Fvtz(k)=u(29);Tvtzx(k)=u(30);Tvhyx(k)=u(31);Tvhxy(k)=u(32);Tvhzy(k)=u(33);Tvhyz(k)=u(34);Tvhxz(k)=u(35);
    Tvtzx(k)=u(36);Tvtyx(k)=u(37);Tvtxy(k)=u(38);Tvtzy(k)=u(39);Tvtyz(k)=u(40);Tvtxz(k)=u(41);wcall(k)=u(42);
    wchead(k)=u(43);wctail(k)=u(44);Fwr(k)=u(45); % All other force and torque components
	% Find what each of them correspond to at the end of this script.
	
end
dat = [allt' U' V' W' avx' avy' avz' rball' Fvx' Fvy' Fvz' Ftz' Fhz' Tvzx' Tvyx' Tmx' Tgx'...
        Tvxy' Tvzy' Tmy' Tgy' Tvyz' Tvxz' Tmz' Tgz' Fvhx' Fvhy' Fvhz' Fvtx' Fvty'...
        Fvtz' Tvhzx' Tvhyx' Tvhxy' Tvhzy' Tvhyz' Tvhxz' Tvtzx' Tvtyx' Tvtxy' Tvtzy' Tvtyz'...
        Tvtxz' wcall' wchead' wctail' Fwr']; % This is the order of the columns in recorded data
		
		
save(strcat('work_data_',run,'.mat'), '-regexp','^(?!(model)$).'); % This saves the workspace in Matlab
if SAVE; save(strcat('./traj_data_',run,'.dat'),'dat','-ascii'); end % This only saves the matrix dat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION OF THE VARIABLES IN dat

% dat matrix is the output of the model and contains 49 different columns for 49 different parameters/variables.
% Each row in data corresponds to a time instant simulated.
% The .dat files provided with this package may contain lower (or higher!) number of columns
% These 49 columns are what can be confirmed to be true. Data on columns larger than 49 should be ignored.

% Number	|	Name	|	Description
% ==========================================================================================
% 1			| 	allt	|	Vector containing all simulated time steps
% 2			|	U		|	Linear velocity in the x-direction
% 3			|	V		|	Linear velocity in the y-direction
% 4			|	W		|	Linear velocity in the z-direction
% 5			|	avx		|	Angular velocity in the x-direction
% 6			|	avy		|	Angular velocity in the y-direction
% 7			|	avz		|	Angular velocity in the z-direction
% 8			|	rb(1)	|	x-coordinate of the swimmer
% 9			|	rb(2)	|	y-coordinate of the swimmer
% 10		|	rb(3)	|	z-coordinate of the swimmer
% 11		|	Fvx		|	Viscous force, x-component
% 12		|	Fvy		|	Viscous force, y-component
% 13		|	Fvz		|	Viscous force, z-component
% 14		|	Ftz		|	Tail weight (buoyancy force subtracted)
% 15 		|	Fhz		|	Head weight (buoyancy force subtracted)
% 16		|	Tvzx	|	Viscous torque in the x-direction due to the forces in the z-direction
% 17		|	Tvyx	|	Viscous torque in the x-direction due to the forces in the y-direction
% 18		|	Tmx		|	Magnetic torque in the x-direction
% 19		|	Tgx		|	Gravity-induced torque in the x-direction
% 20		|	Tvxy	|	Viscous torque in the y-direction due to the forces in the x-direction
% 21		|	Tvzy	|	Viscous torque in the y-direction due to the forces in the z-direction
% 22		|	Tmy		|	Magnetic torque in the y-direction	
% 23		|	Tgy		|	Gravity-induced torque in the y-direction
% 24		|	Tvyz	|	Viscous torque in the z-direction due to the forces in the y-direction
% 25		|	Tvxz	|	Viscous torque in the z-direction due to the forces in the x-direction
% 26		|	Tmz		|	Magnetic torque in the z-direction
% 27		|	Tgz		|	Gravity-induced torque in the z-direction
% 28		|	Fvhx	|	Viscous force in the x-direction (for head only)
% 29		|	Fvhy	|	Viscous force in the y-direction (for head only)
% 30		|	Fvhz	|	Viscous force in the z-direction (for head only)
% 31		|	Fvtx	|	Viscous force in the x-direction (for tail only)
% 32		|	Fvty	|	Viscous force in the y-direction (for tail only)
% 33		|	Fvtz	|	Viscous force in the z-direction (for tail only)
% 34		|	Tvhzx	|	Viscous torque in the x-direction due to the forces in the z-direction (head only)
% 35		|	Tvhyx	|	Viscous torque in the x-direction due to the forces in the y-direction (head only)
% 36		|	Tvhxy	|	Viscous torque in the y-direction due to the forces in the x-direction (head only)
% 37		|	Tvhzy	|	Viscous torque in the y-direction due to the forces in the z-direction (head only)
% 38		|	Tvhyz	|	Viscous torque in the z-direction due to the forces in the y-direction (head only)
% 39		|	Tvhxz	|	Viscous torque in the z-direction due to the forces in the x-direction (head only)
% 40		|	Tvtzx	|	Viscous torque in the x-direction due to the forces in the z-direction (tail only)
% 41		|	Tvtyx	|	Viscous torque in the x-direction due to the forces in the y-direction (tail only)
% 42		|	Tvtxy	|	Viscous torque in the y-direction due to the forces in the x-direction (tail only)
% 43		|	Tvtzy	|	Viscous torque in the y-direction due to the forces in the z-direction (tail only)
% 44		|	Tvtyz	|	Viscous torque in the z-direction due to the forces in the y-direction (tail only)
% 45		|	Tvtxz	|	Viscous torque in the z-direction due to the forces in the x-direction (tail only)
% 46		|	wcall	|	Checks if there is wall contact (whole geometry)
% 47		|	wchead	|	Checks if there is wall contact (head only)
% 48		|	wctail	|	Checks if there is wall contact (tail only)
% 49		|	Fwr		|	Repulsive force (in case of contact) in the radial direction


