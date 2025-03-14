% Part of the code package related to the publication:
%
% Caldag, H. O., & Yesilyurt, S. (2019). Trajectories of magnetically-actuated 
% helical swimmers in cylindricalchannels at low Reynolds numbers. Journal of 
% Fluids and Structures, 90, 164-176. https://doi.org/10.1016/j.jfluidstructs.2019.06.005
%
% This script, called from traj_helical.m, uses the model built in initialbuild.m, changes
% the updated parameters and runs the simulation for the current time step. The solutions are
% returned via V and out. V contains all computed velocities, forces and torques while out is the
% simulation model.

function [V,out] = sim_run(model,nlam,lam,Db,rt,Lh,rch,Lch,rb,R,ma,omega,pstime,...
    grav,Vin,rhohead,rhotail,freq,GRAV,WALL,SAVE,PLOT,run,zint,rot)

model.sol('sol1').clearSolutionData; % Clear solution data from previous time step
model.result.table('tbl1').clearTableData; % Clear table data from previous time step

model.param.set('xb', num2str(0));
model.param.set('yb', num2str(rb(2)));
model.param.set('zb', num2str(zint));
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
model.param.set('pstime',num2str(pstime));
model.param.set('psi', num2str(rot));
model.param.set('Vin', num2str(abs(Vin)));
model.param.set('B0', '3600/(rholiq*lscale*lscale*freq*freq)');
model.param.set('x0', 'xb');
model.param.set('y0', 'yb');
model.param.set('z0', 'zb');
model.param.set('Tmx', 'my*Bz-mz*By');
model.param.set('Tmy', 'mz*Bx-mx*Bz');
model.param.set('Tmz', 'mx*By-my*Bx');
model.param.set('Bx', '0');
model.param.set('By', 'B0*cos(omega*pstime)');
model.param.set('Bz', 'B0*sin(omega*pstime)'); 
model.param.set('x0com', '(rhotail*voltail*x0t+rhohead*volhead*x0b)/(rhotail*voltail+rhohead*volhead)');
model.param.set('Ftz', '-voltail*grav/lscale*tscale*tscale*(rhotail-rho)');
model.param.set('Fhz', '-volhead*grav/lscale*tscale*tscale*(rhohead-rho)');
model.param.set('Fhy', '0');
model.param.set('Fhx', '0');
model.param.set('Fty', '0');
model.param.set('Ftx', '0');
model.param.set('xbt', 'xb+hax*x0t');
model.param.set('ybt', 'yb+hay*x0t');
model.param.set('zbt', 'zb+haz*x0t');
model.param.set('xbh', 'xb+hax*x0b');
model.param.set('ybh', 'yb+hay*x0b');
model.param.set('zbh', 'zb+haz*x0b');
model.param.set('Tgx', '(ybh-y0)*Fhz - (zbh-z0)*Fhy + (ybt-y0)*Ftz - (zbt-z0)*Fty');
model.param.set('Tgy', '(zbh-z0)*Fhx - (xbh-x0)*Fhz + (zbt-z0)*Ftx - (xbt-x0)*Ftz');
model.param.set('Tgz', '(xbh-x0)*Fhy - (ybh-y0)*Fhx + (xbt-x0)*Fty - (ybt-y0)*Ftx');
model.param.set('thet', '-asin(haz)');
model.param.set('phi', 'atan2(hay/cos(thet),hax/cos(thet))');
% These parameters are updated for the new time step.
% Descriptions available at initial_build.m

% Rebuild the geometry
model.component('mod1').geom('geom1').run;

% Redefining variables for contact modelling
model.component('mod1').variable('var1').set('dwall', 'rch-rsw');
model.component('mod1').variable('var1').set('rsw', 'sqrt(y*y+z*z)');
model.component('mod1').variable('var1').set('fwally0', '-fwallc/dwall*velr*cos(theta)');
model.component('mod1').variable('var1').set('fwallz0', '-fwallc/dwall*velr*sin(theta)');
model.component('mod1').variable('var1').set('fwallc0', '6*pi*visc*rt*wallcheck');
model.component('mod1').variable('var1').set('theta', 'atan2(z,y)');
model.component('mod1').variable('var1').set('velr', 'v*cos(theta)+w*sin(theta)');
model.component('mod1').variable('var1').set('wallcheck', '(dwall<delw)');
model.component('mod1').variable('var1').set('fwally', 'fwallr*cos(theta)');
model.component('mod1').variable('var1').set('fwallz', 'fwallr*sin(theta)');
model.component('mod1').variable('var1').set('fwallr', '-max(0,Kspv)*wallcheck-Ksp*(delw-dwall)');
model.component('mod1').variable('var1').set('Tstressr', 'v_lm*cos(theta)+w_lm*sin(theta)');
model.component('mod1').variable('var1').set('Kspv', '-(intop1(Tstressr)/(intop1(wallcheck)+0.1))');
model.component('mod1').variable('var1').selection.named('geom1_csel3_bnd');

% Build mesh
model.component('mod1').mesh('mesh1').run;

% Material setting
model.frame('material1').sorder(1);

% Running the simulation
model.sol('sol1').runAll;

model.result.dataset('cpl1').set('quickx', -1);
model.result.dataset('cpl2').set('quickx', 1);
model.result.dataset('cpl2').set('spacevars', {'cpl1x' 'cpl1y'});
model.result.dataset('cpl2').set('normal', {'cpl1nx' 'cpl1ny' 'cpl1nz'});
model.result.numerical('gev1').set('table', 'tbl1');
model.result.numerical('gev1').setResult;

out = model;
% out will be returned to the main code 

% V contains all the data
V = zeros(45,1);
V(1) = mphglobal(out,'U');
V(2) = mphglobal(out,'V');
V(3) = mphglobal(out,'W');
V(4) = mphglobal(out,'avx');
V(5) = mphglobal(out,'avy');
V(6) = mphglobal(out,'avz');
V(7)=mphglobal(out,'intop1(u_lm)');
V(8)=mphglobal(out,'intop1(v_lm)');
V(9)=mphglobal(out,'intop1(w_lm)');
V(10)=mphglobal(out,'Ftz');
V(11)=mphglobal(out,'Fhz');
V(12)=mphglobal(out,'intop1((y-y0)*(w_lm-fwallz))');
V(13)=mphglobal(out,'intop1((z-z0)*(v_lm-fwally))');
V(14)=mphglobal(out,'Tmx');
V(15)=mphglobal(out,'Tgx');
V(16)=mphglobal(out,'intop1((z-z0)*u_lm)');
V(17)=mphglobal(out,'intop1((x-x0)*(w_lm-fwallz))');
V(18)=mphglobal(out,'Tmy');
V(19)=mphglobal(out,'Tgy');
V(20)=mphglobal(out,'intop1((x-x0)*(v_lm-fwally))');
V(21)=mphglobal(out,'intop1((y-y0)*u_lm)');
V(22)=mphglobal(out,'Tmz');
V(23)=mphglobal(out,'Tgz');
V(24)=mphglobal(out, 'intop2(u_lm)');
V(25)=mphglobal(out, 'intop2(v_lm)');
V(26)=mphglobal(out, 'intop2(w_lm)');
V(27)=mphglobal(out, 'intop3(u_lm)');
V(28)=mphglobal(out, 'intop3(v_lm)');
V(29)=mphglobal(out, 'intop3(w_lm)');
V(30)=mphglobal(out,'intop2((y-y0)*(w_lm))');
V(31)=mphglobal(out,'intop2((z-z0)*(v_lm))');
V(32)=mphglobal(out,'intop2((z-z0)*u_lm)');
V(33)=mphglobal(out,'intop2((x-x0)*(w_lm))');
V(34)=mphglobal(out,'intop2((x-x0)*(v_lm))');
V(35)=mphglobal(out,'intop2((y-y0)*u_lm)');
V(36)=mphglobal(out,'intop3((y-y0)*(w_lm))');
V(37)=mphglobal(out,'intop3((z-z0)*(v_lm))');
V(38)=mphglobal(out,'intop3((z-z0)*u_lm)');
V(39)=mphglobal(out,'intop3((x-x0)*(w_lm))');
V(40)=mphglobal(out,'intop3((x-x0)*(v_lm))');
V(41)=mphglobal(out,'intop3((y-y0)*u_lm)');
V(42)=mphglobal(out,'intop1(wallcheck)');
V(43)=mphglobal(out,'intop2(wallcheck)');
V(44)=mphglobal(out,'intop3(wallcheck)');
V(45)=mphglobal(out, 'intop1(fwallr)');

if SAVE % If the model is to be saved
    fname = sprintf('%s_t=%4.2f.mph',run,pstime); % Save it with the time information
    mphsave(model,fname);
end
if PLOT % This plots the current physical setup on MATLAB
    figure(1)
    model.result.create('pg1', 'PlotGroup3D');
    model.result('pg1').create('surf1', 'Surface');
    model.result('pg1').create('arws1', 'ArrowSurface');
    model.result('pg1').create('arws2', 'ArrowSurface');
    model.result('pg1').set('titletype', 'none');
    model.result('pg1').set('frametype', 'spatial');
    model.result('pg1').feature('surf1').set('data', 'surf1');
    model.result('pg1').feature('surf1').set('expr', 'p');
    model.result('pg1').feature('surf1').set('unit', 'Pa');
    model.result('pg1').feature('surf1').set('descr', 'Pressure');
    model.result('pg1').feature('surf1').set('titletype', 'none');
    model.result('pg1').feature('surf1').set('rangecoloractive', true);
    model.result('pg1').feature('surf1').set('rangecolormin', '-1e5');
    model.result('pg1').feature('surf1').set('rangecolormax', '1e5');
    model.result('pg1').feature('surf1').set('resolution', 'normal');
    model.result('pg1').feature('arws1').set('data', 'cpl1');
    model.result('pg1').feature('arws1').set('titletype', 'none');
    model.result('pg1').feature('arws1').set('scale', 0.2);
    model.result('pg1').feature('arws1').set('scaleactive', true);
    model.result('pg1').feature('arws1').set('arrowcount', 50);
    model.result('pg1').feature('arws2').set('data', 'cpl2');
    model.result('pg1').feature('arws2').set('titletype', 'none');
    model.result('pg1').feature('arws2').set('scale', 0.2);
    model.result('pg1').feature('arws2').set('scaleactive', true);
    model.result('pg1').feature('arws2').set('arrowcount', 50);
    model.result('pg1').set('view', 'view1');
    fname = sprintf('./%s/t=%4.2f.jpg',run,pstime);
    mphplot(model,'pg1');
    drawnow;pause(.01);
    saveas(gcf,fname);    % saving the frame, later can use these to generate an animation
end
