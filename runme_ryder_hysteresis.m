%% Project Description
% This project focuses on Ryder glacier, Greenland, and the forcings needed to
% regrow its ice tongue after an eventual future loss (here we assume simulations
% start in the year 2100). You should try changing the SMB anomaly, melt rates, 
% and calving parameters both individually and together to see which circumstances allow for regrowth
% of the ice tongue. The SMB anomaly is applied to the mean historical
% (1960-2014) forcing, so starting with an anomaly of 0 will let you know if reverting
% to a historical climate is enough to regrow the ice tongue.
% Reading starting point: https://www.nature.com/articles/s41467-022-29529-5

flowline = 'Exp/ryder_flowline.exp';

%% Toggles

smb_anomaly = 5; % Anomaly added to the mean historical SMB field in mm water equivalent
melt_rate = 25; %Default of 25 m/yr, change to whatever you would like
calving_threshold_grounded = 5e5; %Default 5e5, change to whatever you would like 
calving_threshold_floating = 200e3; %Default 200e3, change to whatever you would like

number_of_years = 2; %Length of simulation from 2100

run_name = 'ryder_SMB1_melt25_calving5e5_200e3';

%% Transient run - do not edit this section

load('./Model_Data/Ryder_Spinup.mat');
md = transientrestart(md);


 %Dont use damage model
md.damage.D=zeros(md.mesh.numberofvertices,1);
md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);

%Dont use thermal model
md.transient.isthermal=0;
md.thermal.spctemperature=NaN*ones(md.mesh.numberofvertices,1);

%Additional options
md.inversion.iscontrol=0;

disp('Reading and interpolating SMB data')

% md.smb.mass_balance = [];
% smbMAR = [];
% 
% for yy=1:(40*12) % Monthly data
%     smboutput = interpMAR_monthly(md.mesh.x,md.mesh.y,'SMB',yy, './Model_Data/MARv3.11.3-historical-combined.nc');
%     smbMAR = [smbMAR smboutput];
%     %progress = sprintf('Read %d timesteps out of %d',yy, nyrs_smb*12);
%     %disp(progress)
% end
% 
% pos=find(smbMAR==-9999);
% smbMAR(pos)=0.0;
%       
% md.smb.mass_balance = (smbMAR+smb_anomaly)/1000*12*(md.materials.rho_freshwater/md.materials.rho_ice); 
% md.smb.mass_balance = mean(md.smb.mass_balance,2);

%Read in pre-interpolated SMB
load('Model_Data/ryder_SMB.mat');
md.smb.mass_balance = smb_readin+smb_anomaly;

%Make sure bed is below base
pos=find(md.geometry.bed>md.geometry.base);
md.geometry.base(pos)=md.geometry.bed(pos);

%Recalculate surface
md.geometry.surface=md.geometry.base+md.geometry.thickness;

%Front and GL options
md.transient.ismovingfront=1;
md.transient.isgroundingline=1;
md.groundingline.migration='SubelementMigration';
md.groundingline.melt_interpolation='NoMeltOnPartiallyFloating';
md.groundingline.friction_interpolation='SubelementFriction1';	

isresetlevelset = 1;
md.levelset.spclevelset=NaN(md.mesh.numberofvertices,1);
%%% Set levelset options %%%
if isresetlevelset == 1
    md.levelset.spclevelset=NaN(md.mesh.numberofvertices,1);

    %Get mask from BedMachine
    M = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask','nearest','./Model_Data/BedMachineGreenland-2022-03-17.nc');
    pos=find(M<2);
    md.mask.ice_levelset(pos)=+1; %set levelset for no-ice vertices
    %remove 0 in ice_levelset (advection will fail if used)
    md.mask.ice_levelset(find(md.mask.ice_levelset==0))=-1;

    %make it a signed distance
    md.mask.ice_levelset = reinitializelevelset(md,md.mask.ice_levelset);

    % Reset levelset boundary conditions on domain boundary
    md.levelset.spclevelset(find(md.mesh.vertexonboundary)) = md.mask.ice_levelset(find(md.mesh.vertexonboundary));
else 
    %dont touch the spclevelset, just keep what is from the previous model and do nothing here
end

icelandspc=1;
if icelandspc == 1
    md.levelset.spclevelset=NaN(md.mesh.numberofvertices, 1);
    pos = find_iceLandBoundary(md, 1); %1=is2D
    md.levelset.spclevelset(pos)=-1;
end


%Basal Melt options
% Fixed melt
md.basalforcings=linearbasalforcings();
%md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.basalforcings.deepwater_melting_rate = melt_rate;
md.basalforcings.deepwater_elevation = -300;
md.basalforcings.upperwater_melting_rate = 0;
md.basalforcings.upperwater_elevation = -100;
md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);
md.basalforcings.geothermalflux=interpSeaRISE_new(md.mesh.x,md.mesh.y,'bheatflx');


 md.levelset.kill_icebergs=1;
%md.levelset.migration_max=10000; % -- maximum allowed migration rate (m/a)

    %Calving options
if md.transient.ismovingfront==1
    md.calving=calvingvonmises(); %activate von mises calving law

    %Stress threshold
    md.calving.stress_threshold_groundedice= calving_threshold_grounded; %default 1 MPa = 1e6 Pa
    md.calving.stress_threshold_floatingice= calving_threshold_floating; %default Petermann 300 kPa, default ISSM 150 kPa
    disp(['Calving sigma_max floatingice set to ' num2str(md.calving.stress_threshold_floatingice./1000)  ' kPa'])

    md.calving.min_thickness=50; %m, default NaN

    %Define calving rate and melt rate (only effective if ismovingfront==1)
    md.frontalforcings.meltingrate=melt_rate*ones(md.mesh.numberofvertices, 1); %only effective if front is grounded
end


%Timestepping options
md.timestepping.cycle_forcing = 1;
md.timestepping = timestepping();
md.timestepping.time_step = (1/12); % Adaptive test suggested 0.2 would even be ok
md.settings.output_frequency = 1; %montlhy
% 	md.settings.output_frequency=1; %1: every tstep; 5: every fifth tstep, etc (for debugging)
disp(['Setting fixed time step to ' num2str(md.timestepping.time_step) ' yrs'])
%md.timestepping=timesteppingadaptive();
%md.timestepping.time_step_max=5;
%md.timestepping.time_step_min=0.01;
md.timestepping.start_time = 2100;
md.timestepping.final_time = md.timestepping.start_time + number_of_years;



%Output options
md.transient.requested_outputs={'TotalSmb','SmbMassBalance',...
'IceVolume','IceVolumeAboveFloatation',...
'IceVolumeAboveFloatationScaled','GroundedAreaScaled',...
'MaskOceanLevelset','MaskIceLevelset',...
'FloatingAreaScaled','IceMass',...
'GroundedArea','FloatingArea','TotalFloatingBmb',...
'BasalforcingsFloatingiceMeltingRate',...
'TotalCalvingFluxLevelset',... %Gt/r
'GroundinglineMassFlux',... %Gt/yr
'CalvingMeltingrate','TotalCalvingMeltingFluxLevelset','IcefrontMassFluxLevelset',...
'TotalCalvingFluxLevelset','TotalGroundedBmb',...
'Calvingratex','Calvingratey','CalvingCalvingrate','SigmaVM'};

md.verbose=verbose('solution',true,'module',false,'convergence',false);

md.toolkits=toolkits();
md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
md.cluster=generic('name',oshostname,'np',4);

%Solve
md=solve(md,'Transient');

%savemodel(org,md);
save(['Outputs/' run_name],'md','-v7.3')

plotmodel(md, 'data', md.results.TransientSolution(1).Vel, 'data', md.results.TransientSolution(end).Vel, ...
    'mask', md.results.TransientSolution(1).MaskIceLevelset<0, 'mask#2-4', md.results.TransientSolution(end).MaskIceLevelset<0, ...
    'data', md.results.TransientSolution(end).Vel-md.results.TransientSolution(1).Vel,...
    'data', md.results.TransientSolution(end).Thickness-md.results.TransientSolution(1).Thickness,...
    'caxis#3-4', [-250 250], 'ncols', 4,'caxis#1-2', [0 max(md.results.TransientSolution(end).Vel)], ...
    'title','Initial Velocity (m/yr)' , 'title','Final Velocity (m/yr)' , 'title', 'End velocity - Starting velocity (m/yr)', ...
    'title', 'End thickness - Starting thickness (m)')








