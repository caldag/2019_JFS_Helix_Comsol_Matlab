% Part of the code package related to the publication:
%
% Caldag, H. O., & Yesilyurt, S. (2019). Trajectories of magnetically-actuated 
% helical swimmers in cylindricalchannels at low Reynolds numbers. Journal of 
% Fluids and Structures, 90, 164-176. https://doi.org/10.1016/j.jfluidstructs.2019.06.005
%
% This script, called from traj_helical.m, builds the initial FEM model to be used
% in the simulation
%
% Majority of the code is generated automatically from a complete COMSOL Multiphysics model file.

function model = initial_build(nlam,lam,Db,rt,Lh,rch,Lch,rb,R,ma,omega,pstime,...
    grav,Vin,rhohead,rhotail,freq,GRAV,WALL)

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.hist.disable

model.modelPath('C:\'); % Insert the address of storage here

model.label('swimmer_new.mph'); % FEM file name


% Below we assign all the parameters sent from traj_helical.m
% See traj_helical.m for descriptions

model.param.set('freq', num2str(freq));
model.param.set('Db', num2str(Db));
model.param.set('nlam',num2str(nlam));
model.param.set('rt',num2str(rt));
model.param.set('lamt',num2str(lam));
model.param.set('omega',num2str(omega));
model.param.set('xb', num2str(0));
model.param.set('yb', num2str(rb(2)));
model.param.set('zb', num2str(rb(3)));
model.param.set('hax', num2str(R(1,1)));
model.param.set('hay', num2str(R(2,1)));
model.param.set('haz', num2str(R(3,1)));
model.param.set('jax', num2str(R(1,2)));
model.param.set('jay', num2str(R(2,2)));
model.param.set('jaz', num2str(R(3,2)));
model.param.set('kax', num2str(R(1,3)));
model.param.set('kay', num2str(R(2,3)));
model.param.set('kaz', num2str(R(3,3)));
model.param.set('mx', num2str(ma(1)));
model.param.set('my', num2str(ma(2)));
model.param.set('mz', num2str(ma(3)));
model.param.set('avx', num2str(2*pi));
model.param.set('pstime',num2str(pstime));
model.param.set('rhotail', num2str(rhotail));
model.param.set('rhohead', num2str(rhohead));
model.param.set('rch', num2str(rch));
model.param.set('Lch', num2str(Lch));
model.param.set('Lh', num2str(Lh));
model.param.set('Vin', num2str(abs(Vin)));
model.param.set('grav',num2str(grav));
model.param.set('IFGRAV',num2str(GRAV));
model.param.set('IFWALL',num2str(WALL));
model.param.set('B', 'Db/2-rt*1.2');
model.param.set('rho', '1');
model.param.set('x0', 'xb');
model.param.set('y0', 'yb');
model.param.set('z0', 'zb');
model.param.set('Lh', '1.5');
model.param.set('Btail', 'Db/2-rt*1.1');

model.param.set('lscale', '1e-3'); % Length scale of the system
model.param.set('rholiq', '1264'); % Liquid density
model.param.set('muliq', '1.412'); % Liquid viscosity
model.param.set('Rescale', 'rholiq*lscale*lscale*freq/muliq'); % Reynolds number
model.param.set('Tmx', 'my*Bz-mz*By');
model.param.set('Tmy', 'mz*Bx-mx*Bz');
model.param.set('Tmz', 'mx*By-my*Bx'); % Cartesian components of the magnetic torque
model.param.set('Bx', '0');
model.param.set('By', 'B0*cos(omega*pstime)');
model.param.set('Bz', 'B0*sin(omega*pstime)'); % Cartesian components of the magnetic field
model.param.set('voltail', 'pi*rt*rt*sqrt(4*pi*pi*B*B+lamt*lamt)*nlam'); % Tail volume
model.param.set('volhead', 'pi*Db^2/4*Lh'); % Head volume
model.param.set('x0b', 'Lh/2'); % Head center-of-mass
model.param.set('Lt', 'nlam*lamt-rt'); % tail length
model.param.set('x0t', 'Lt/2+Lh'); % Center-of-mass of the tail in the x- direction (assuming fully aligned with x-axis)
model.param.set('x0com', '(rhotail*voltail*x0t+rhohead*volhead*x0b)/(rhotail*voltail+rhohead*volhead)'); % Total center-of-mass

model.param.set('visc', '1/Rescale'); % Viscosity as a non-dimensional measure
model.param.set('delw', 'rt'); % Wall contact threshold. If the swimmer is delw close to the wall, it will experience a repulsive force
model.param.set('Ksp', '0');
model.param.set('Cwall', '1e3');
model.param.set('tscale', '1/freq'); % Time scale
model.param.set('Ftz', '-voltail*grav/lscale*tscale*tscale*(rhotail-rho)');
model.param.set('Fhz', '-volhead*grav/lscale*tscale*tscale*(rhohead-rho)'); % Gravitational forces on the tail and the head
model.param.set('Fhy', '0');
model.param.set('Fhx', '0');
model.param.set('Fty', '0');
model.param.set('Ftx', '0'); % Other external forces are zero
model.param.set('xbt', 'xb+hax*x0t');
model.param.set('ybt', 'yb+hay*x0t');
model.param.set('zbt', 'zb+haz*x0t');
model.param.set('xbh', 'xb+hax*x0b');
model.param.set('ybh', 'yb+hay*x0b');
model.param.set('zbh', 'zb+haz*x0b'); % Actual centers of mass of head and tail (orientation in consideration)
model.param.set('Tgx', '(ybh-y0)*Fhz - (zbh-z0)*Fhy + (ybt-y0)*Ftz - (zbt-z0)*Fty');
model.param.set('Tgy', '(zbh-z0)*Fhx - (xbh-x0)*Fhz + (zbt-z0)*Ftx - (xbt-x0)*Ftz');
model.param.set('Tgz', '(xbh-x0)*Fhy - (ybh-y0)*Fhx + (xbt-x0)*Fty - (ybt-y0)*Ftx'); % Gravitational torque components

model.param.set('psi', 'atan2(jaz,kaz)');
model.param.set('thet', '-asin(haz)');
model.param.set('phi', 'atan2(hay/cos(thet),hax/cos(thet))'); % Orientation angles

% The section below establishes the geometric setup
model.component.create('mod1', false);
model.component('mod1').geom.create('geom1', 3);
model.result.table.create('tbl1', 'Table');
model.component('mod1').mesh.create('mesh1');
model.component('mod1').geom('geom1').selection.create('csel1', 'CumulativeSelection');
model.component('mod1').geom('geom1').selection('csel1').label('head');
model.component('mod1').geom('geom1').selection.create('csel2', 'CumulativeSelection');
model.component('mod1').geom('geom1').selection('csel2').label('tail');
model.component('mod1').geom('geom1').selection.create('csel3', 'CumulativeSelection');
model.component('mod1').geom('geom1').selection('csel3').label('swimmer');
model.component('mod1').geom('geom1').create('wp1', 'WorkPlane');
model.component('mod1').geom('geom1').feature('wp1').set('quickplane', 'yx');
model.component('mod1').geom('geom1').feature('wp1').set('unite', true);
model.component('mod1').geom('geom1').feature('wp1').geom.create('r1', 'Rectangle');
model.component('mod1').geom('geom1').feature('wp1').geom.feature('r1').set('pos', [0 0]);
model.component('mod1').geom('geom1').feature('wp1').geom.feature('r1').set('size', {'Db/2' 'Lh'});
model.component('mod1').geom('geom1').feature('wp1').geom.create('fil1', 'Fillet');
model.component('mod1').geom('geom1').feature('wp1').geom.feature('fil1').set('radius', 'rt');
model.component('mod1').geom('geom1').feature('wp1').geom.feature('fil1').selection('point').set('r1(1)', [2 3]);
model.component('mod1').geom('geom1').create('rev1', 'Revolve');
model.component('mod1').geom('geom1').feature('rev1').set('contributeto', 'csel1');
model.component('mod1').geom('geom1').feature('rev1').set('selresult', true);
model.component('mod1').geom('geom1').feature('rev1').set('selresultshow', 'all');
model.component('mod1').geom('geom1').feature('rev1').set('angtype', 'full');
model.component('mod1').geom('geom1').feature('rev1').selection('input').set({'wp1'});
model.component('mod1').geom('geom1').create('hel2', 'Helix');
model.component('mod1').geom('geom1').feature('hel2').set('contributeto', 'csel2');
model.component('mod1').geom('geom1').feature('hel2').set('selresult', true);
model.component('mod1').geom('geom1').feature('hel2').set('selresultshow', 'all');
model.component('mod1').geom('geom1').feature('hel2').set('pos', {'Lh-rt' '0' '0'});
model.component('mod1').geom('geom1').feature('hel2').set('axis', [1 0 0]);
model.component('mod1').geom('geom1').feature('hel2').set('turns', 'nlam');
model.component('mod1').geom('geom1').feature('hel2').set('rmaj', 'Btail');
model.component('mod1').geom('geom1').feature('hel2').set('rmin', 'rt');
model.component('mod1').geom('geom1').feature('hel2').set('axialpitch', 'lamt');
model.component('mod1').geom('geom1').feature('hel2').set('chirality', 'left');
model.component('mod1').geom('geom1').feature('hel2').set('endcaps', 'perpspine');
model.component('mod1').geom('geom1').create('sph1', 'Sphere');
model.component('mod1').geom('geom1').feature('sph1').set('contributeto', 'csel2');
model.component('mod1').geom('geom1').feature('sph1').set('selresult', true);
model.component('mod1').geom('geom1').feature('sph1').set('selresultshow', 'all');
model.component('mod1').geom('geom1').feature('sph1').set('pos', {'Lh-rt+nlam*lamt' '(Btail)*sin((nlam)*2*pi)' '-(Btail)*cos((nlam)*2*pi)'});
model.component('mod1').geom('geom1').feature('sph1').set('r', 'rt');
model.component('mod1').geom('geom1').create('uni2', 'Union');
model.component('mod1').geom('geom1').feature('uni2').set('contributeto', 'csel3');
model.component('mod1').geom('geom1').feature('uni2').set('selresult', true);
model.component('mod1').geom('geom1').feature('uni2').set('selresultshow', 'all');
model.component('mod1').geom('geom1').feature('uni2').set('intbnd', false);
model.component('mod1').geom('geom1').feature('uni2').selection('input').set({'hel2' 'rev1' 'sph1'});
model.component('mod1').geom('geom1').create('mov1', 'Move');
model.component('mod1').geom('geom1').feature('mov1').set('displx', '-x0com');
model.component('mod1').geom('geom1').feature('mov1').selection('input').set({'uni2'});
model.component('mod1').geom('geom1').create('rot1', 'Rotate');
model.component('mod1').geom('geom1').feature('rot1').set('rot', 'psi');
model.component('mod1').geom('geom1').feature('rot1').set('pos', [0 0 0]);
model.component('mod1').geom('geom1').feature('rot1').set('axis', [1 0 0]);
model.component('mod1').geom('geom1').feature('rot1').selection('input').set({'mov1'});
model.component('mod1').geom('geom1').create('rot2', 'Rotate');
model.component('mod1').geom('geom1').feature('rot2').set('rot', 'thet');
model.component('mod1').geom('geom1').feature('rot2').set('pos', [0 0 0]);
model.component('mod1').geom('geom1').feature('rot2').set('axis', [0 1 0]);
model.component('mod1').geom('geom1').feature('rot2').selection('input').set({'rot1'});
model.component('mod1').geom('geom1').create('rot3', 'Rotate');
model.component('mod1').geom('geom1').feature('rot3').set('rot', 'phi');
model.component('mod1').geom('geom1').feature('rot3').set('pos', [0 0 0]);
model.component('mod1').geom('geom1').feature('rot3').set('axis', [0 0 1]);
model.component('mod1').geom('geom1').feature('rot3').selection('input').set({'rot2'});
model.component('mod1').geom('geom1').create('mov2', 'Move');
model.component('mod1').geom('geom1').feature('mov2').set('displx', 'xb');
model.component('mod1').geom('geom1').feature('mov2').set('disply', 'yb');
model.component('mod1').geom('geom1').feature('mov2').set('displz', 'zb');
model.component('mod1').geom('geom1').feature('mov2').selection('input').set({'rot3'});
model.component('mod1').geom('geom1').create('cyl2', 'Cylinder');
model.component('mod1').geom('geom1').feature('cyl2').set('pos', {'-Lch/2+(Lt+Lh)/2-x0com' '0' '0'});
model.component('mod1').geom('geom1').feature('cyl2').set('axis', [1 0 0]);
model.component('mod1').geom('geom1').feature('cyl2').set('r', 'rch');
model.component('mod1').geom('geom1').feature('cyl2').set('h', 'Lch');
model.component('mod1').geom('geom1').create('dif1', 'Difference');
model.component('mod1').geom('geom1').feature('dif1').selection('input').set({'cyl2'});
model.component('mod1').geom('geom1').feature('dif1').selection('input2').set({'mov2'});
model.component('mod1').geom('geom1').run;

% Below are the variables of the model, mostly related to wall contact force
model.component('mod1').variable.create('var1');
model.component('mod1').variable('var1').set('dwall', 'rch-rsw'); % Distance to wall
model.component('mod1').variable('var1').set('rsw', 'sqrt(y*y+z*z)'); % Radial position
model.component('mod1').variable('var1').set('theta', 'atan2(z,y)'); % Angular position on the radial plane
model.component('mod1').variable('var1').set('velr', 'v*cos(theta)+w*sin(theta)'); % Radial velocity
model.component('mod1').variable('var1').set('wallcheck', '(dwall<delw)'); % Boolean showing if repulsive force should be exerted
model.component('mod1').variable('var1').set('fwally', 'fwallr*cos(theta)'); % Cartesian components of the repulsive forces
model.component('mod1').variable('var1').set('fwallz', 'fwallr*sin(theta)');
model.component('mod1').variable('var1').set('fwallr', '-max(0,Kspv)*wallcheck-Ksp*(delw-dwall)'); % Radial repulsive force
model.component('mod1').variable('var1').set('Tstressr', 'v_lm*cos(theta)+w_lm*sin(theta)'); % Radial stress exerted by the swimmer to the fluid
model.component('mod1').variable('var1').set('Kspv', '-(intop1(Tstressr)/(intop1(wallcheck)+0.1))'); 
% Integrate the stresses on regions where the swimmer is close to the wall, this will be the repulsive force
% A 0.1 is added to the denominator to prevent division by zero
model.component('mod1').variable('var1').selection.named('geom1_csel3_bnd');

model.component('mod1').cpl.create('intop1', 'Integration');
model.component('mod1').cpl.create('intop2', 'Integration');
model.component('mod1').cpl.create('intop3', 'Integration');
model.component('mod1').cpl('intop1').selection.named('geom1_csel3_bnd');
model.component('mod1').cpl('intop2').selection.named('geom1_csel1_bnd');
model.component('mod1').cpl('intop3').selection.named('geom1_csel2_bnd');

% Introducing the physical system (flow, pressure constraint etc.)
model.component('mod1').physics.create('ge', 'GlobalEquations', 'geom1');
model.component('mod1').physics.create('spf2', 'CreepingFlow', 'geom1');
model.component('mod1').physics('spf2').identifier('spf2');
model.component('mod1').physics('spf2').field('turbulentkineticenergy').field('k2');
model.component('mod1').physics('spf2').field('turbulentdissipationrate').field('ep2');
model.component('mod1').physics('spf2').field('specificdissipationrate').field('om2');
model.component('mod1').physics('spf2').field('reciprocallength').field('G2');
model.component('mod1').physics('spf2').field('correctedvelocity').field('uc2');
model.component('mod1').physics('spf2').field('correctedvelocity').component({'uc2x' 'uc2y' 'uc2z'});
model.component('mod1').physics('spf2').field('correctedpressure').field('pc2');
model.component('mod1').physics('spf2').field('turbulentkinematicviscosity').field('nutilde2');
model.component('mod1').physics('spf2').field('dimensionless1').field('yPlus2');
model.component('mod1').physics('spf2').field('dimensionless2').field('uPlus2');
model.component('mod1').physics('spf2').field('dimensionless3').field('zeta2');
model.component('mod1').physics('spf2').field('dimensionless4').field('alpha2');
% Define inlet only if the inlet velocity is non-zero
if Vin>=0
model.component('mod1').physics('spf2').create('inl1', 'InletBoundary', 2);
model.component('mod1').physics('spf2').feature('inl1').selection.set([1]);
model.component('mod1').physics('spf2').create('out1', 'OutletBoundary', 2);
model.component('mod1').physics('spf2').feature('out1').selection.set([36]);
else
    model.component('mod1').physics('spf2').create('inl1', 'InletBoundary', 2);
model.component('mod1').physics('spf2').feature('inl1').selection.set([36]);
model.component('mod1').physics('spf2').create('out1', 'OutletBoundary', 2);
model.component('mod1').physics('spf2').feature('out1').selection.set([1]);
end
model.component('mod1').physics('spf2').create('wallbc2', 'WallBC', 2);
model.component('mod1').physics('spf2').feature('wallbc2').selection.named('geom1_csel3_bnd');

% Mesh building
model.component('mod1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('mod1').mesh('mesh1').create('cr1', 'CornerRefinement');
model.component('mod1').mesh('mesh1').create('ftet1', 'FreeTet');
model.component('mod1').mesh('mesh1').create('bl1', 'BndLayer');
model.component('mod1').mesh('mesh1').feature('ftri1').selection.named('geom1_csel3_bnd');
model.component('mod1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('mod1').mesh('mesh1').feature('cr1').selection.geom('geom1', 3);
model.component('mod1').mesh('mesh1').feature('cr1').selection.set([1]);
model.component('mod1').mesh('mesh1').feature('ftet1').create('size1', 'Size');
model.component('mod1').mesh('mesh1').feature('ftet1').feature('size1').selection.geom('geom1', 3);
model.component('mod1').mesh('mesh1').feature('ftet1').feature('size1').selection.set([1]);
model.component('mod1').mesh('mesh1').feature('bl1').selection.geom('geom1', 3);
model.component('mod1').mesh('mesh1').feature('bl1').selection.set([1]);
model.component('mod1').mesh('mesh1').feature('bl1').create('blp1', 'BndLayerProp');
model.component('mod1').mesh('mesh1').feature('bl1').feature('blp1').selection.set([2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35]);

% This is where the evaluated results will be kept
model.result.table('tbl1').comments('Global Evaluation 1 (U, V, W, avx, avy, avz)');

model.capeopen.label('Thermodynamics Package');

% This is where we implement force and torque-free swimming
% The swimming is 'torque-free' if we put the magnetic torque alongside other torques...
model.component('mod1').physics('ge').feature('ge1').set('name', {'U'; 'V'; 'W'; 'avy'; 'avz'});
model.component('mod1').physics('ge').feature('ge1').set('equation', {...
    'intop1(u_lm)'; ... % x direction forces (fluid)
    'intop1(v_lm-fwally*IFWALL)'; ... % y direction forces (contact and fluid)
    'intop1(w_lm-fwallz*IFWALL)+(Ftz+Fhz)*IFGRAV'; ... % z direction forces (contact, fluid and gravity)
    'intop1((z-z0)*u_lm) - intop1((x-x0)*(w_lm-fwallz*IFWALL)) - Tmy - Tgy*IFGRAV';... % y direction torques (fluid, contact, magnetic and gravity)
    'intop1((x-x0)*(v_lm-fwally*IFWALL)) - intop1((y-y0)*u_lm) - Tmz - Tgz*IFGRAV'}); % z direction torques (fluid, contact, magnetic and gravity)
model.component('mod1').physics('ge').feature('ge1').set('initialValueU', [0; 0; 0; 0; 0; 0]);
model.component('mod1').physics('ge').feature('ge1').set('initialValueUt', [0; 0; 0; 0; 0; 0]);
model.component('mod1').physics('ge').feature('ge1').set('description', {''; ''; ''; ''; ''; ''});
model.component('mod1').physics('spf2').prop('ShapeProperty').set('order_fluid', 1);
model.component('mod1').physics('spf2').prop('AdvancedSettingProperty').set('locCFL', 'spf2.localCFLvalue');
model.component('mod1').physics('spf2').feature('fp1').set('rho', 'rho');
model.component('mod1').physics('spf2').feature('fp1').set('mu', 'visc');
model.component('mod1').physics('spf2').feature('init1').set('nutilde_init', 'spf2.nutildeinit');
model.component('mod1').physics('spf2').feature('init1').set('G_init', 'spf2.G0');
model.component('mod1').physics('spf2').feature('init1').set('k_init', 'spf2.kinit');
model.component('mod1').physics('spf2').feature('init1').set('ep_init', 'spf2.epinit');
model.component('mod1').physics('spf2').feature('init1').set('om_init', 'spf2.omInit');
model.component('mod1').physics('spf2').feature('init1').set('yPlus_init', 'spf2.yPlusinit');
model.component('mod1').physics('spf2').feature('init1').set('uPlus_init', 'spf2.uPlusinit');
model.component('mod1').physics('spf2').feature('inl1').set('BoundaryCondition', 'LaminarInflow');
model.component('mod1').physics('spf2').feature('inl1').set('Uav', 'Vin');
model.component('mod1').physics('spf2').feature('wallbc2').set('utr', {'U+(z-z0)*avy-(y-y0)*avz'; 'V+(x-x0)*avz-(z-z0)*(avx)'; 'W+avx*(y-y0)-avy*(x-x0)'});
model.component('mod1').physics('spf2').feature('wallbc2').set('TranslationalVelocityOption', 'Manual');
model.component('mod1').physics('spf2').feature('wallbc2').set('weakConstraints', true);

model.component('mod1').mesh('mesh1').feature('size').set('hauto', 7);
model.component('mod1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('mod1').mesh('mesh1').feature('size').set('table', 'cfd');
model.component('mod1').mesh('mesh1').feature('size').set('hmax', 0.25);
model.component('mod1').mesh('mesh1').feature('size').set('hmin', 0.0717);
model.component('mod1').mesh('mesh1').feature('size').set('hnarrow', 0.5);
model.component('mod1').mesh('mesh1').feature('size').set('hgrad', 1.25);
model.component('mod1').mesh('mesh1').feature('ftri1').feature('size1').set('custom', 'on');
model.component('mod1').mesh('mesh1').feature('ftri1').feature('size1').set('hmax', '.25');
model.component('mod1').mesh('mesh1').feature('ftri1').feature('size1').set('hmaxactive', true);
model.component('mod1').mesh('mesh1').feature('ftri1').feature('size1').set('hmin', 0.27);
model.component('mod1').mesh('mesh1').feature('ftri1').feature('size1').set('hminactive', false);
model.component('mod1').mesh('mesh1').feature('cr1').selection('boundary').set([2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35]);
model.component('mod1').mesh('mesh1').feature('ftet1').feature('size1').set('hauto', 3);
model.component('mod1').mesh('mesh1').feature('bl1').set('sharpcorners', 'trim');
model.component('mod1').mesh('mesh1').feature('bl1').feature('blp1').set('blnlayers', 1);
model.component('mod1').mesh('mesh1').feature('bl1').feature('blp1').set('blhminfact', 10);
model.component('mod1').mesh('mesh1').run;

% Material assignment

model.frame('material1').sorder(1);

model.component('mod1').physics('spf2').feature('fp1').set('rho_mat', 'userdef');
model.component('mod1').physics('spf2').feature('fp1').set('mu_mat', 'userdef');

% Establishing the solver

model.study.create('std1');
model.study('std1').create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.result.dataset.create('surf1', 'Surface');
model.result.dataset.create('cpl1', 'CutPlane');
model.result.dataset.create('surf2', 'Surface');
model.result.dataset.create('cpl2', 'CutPlane');
model.result.dataset('surf1').selection.named('geom1_csel3_bnd');
model.result.dataset('surf2').selection.set([2]);
model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical('gev1').set('probetag', 'none');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('s1').feature('dDef').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('dDef').set('ooc', false);

end