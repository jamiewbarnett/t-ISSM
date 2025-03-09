%% Project Description
% This project focuses on Jakobshavn glacier, Greenland, and how it
% may respond to changes in future climate. You should run several
% simulations with various climatic forcings to look at the sensitivity to
% changes in SMB, frontal melt, and the von Mises calving stress threshold.
% Run some simulations where you just change one parameter at a time (e.g.
% keep the default melt and calving but change SMB), and some simulations
% where you change several parameters simultaneously.



flowline = 'Exp/Jakobshavn_flowline.exp';

%% Toggles

smb_scenario = 'ssp245'; % Alternatives: 'ssp245' or 'ssp585'
melt_rate = 4*365; %Default of 4 m/day, change to whatever you would like
calving_threshold = 1.5e6; %Default 1.5e6, change to whatever you would like 

number_of_years = 10; %Length of simulation from 2024

run_name = 'jakobshavn_ssp245_melt4_calving1e7';

%% Transient run - do not edit this section

load('./Model_Data/Jakobshavn_Spinup.mat');
md = transientrestart(md);


%Dont use damage model
md.damage.D=zeros(md.mesh.numberofvertices,1);
md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);

%Dont use thermal model
md.transient.isthermal=0;
md.thermal.spctemperature=NaN*ones(md.mesh.numberofvertices,1);

%Additional options
md.inversion.iscontrol=0;

%SMB interpolation
disp('Reading and interpolating SMB data')
md.smb.mass_balance = [];
smbMAR = [];

switch smb_scenario
    case{'ssp245'} 
        smb_file='./Model_Data/MARv3.11.3-ssp245-combined.nc';
    case{'ssp585'}
        smb_file='./Model_Data/MARv3.11.3-ssp585-combined.nc';
end

for yy=121:(((2024+number_of_years)-2024)+10)*12 % Monthly data
    smboutput = interpMAR_monthly(md.mesh.x,md.mesh.y,'SMB',yy, smb_file);
    smbMAR = [smbMAR smboutput];
end

pos=find(smbMAR==-9999);
smbMAR(pos)=0.0;
  
md.smb.mass_balance = [smbMAR/1000*12*(md.materials.rho_freshwater/md.materials.rho_ice); ...
[2024:1/12:(2024+number_of_years)-(1/12)]]; % Monthly transient forcing
mean_SMB = mean(smbMAR/1000*12*(md.materials.rho_freshwater/md.materials.rho_ice), 2);

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

% if icelandspc == 1
% md.levelset.spclevelset=NaN(md.mesh.numberofvertices, 1);
% pos = find_iceLandBoundary(md, 1); %1=is2D
% md.levelset.spclevelset(pos)=-1;
% end

md.levelset.kill_icebergs=1;
md.levelset.migration_max=10000; % -- maximum allowed migration rate (m/a)

%Basal Melt options
% Fixed melt
md.basalforcings=linearbasalforcings();
%md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.basalforcings.deepwater_melting_rate = melt_rate;
md.basalforcings.deepwater_elevation = -600;
md.basalforcings.upperwater_melting_rate = 0;
md.basalforcings.upperwater_elevation = -50;
md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);
md.basalforcings.geothermalflux=interpSeaRISE_new(md.mesh.x,md.mesh.y,'bheatflx');

%Define calving rate and melt rate (only effective if ismovingfront==1)


melt_year = [(0.2*365)/2 (0.2*365)/2 melt_rate/2 (0.2*365)/2];     %Winter melt 0.6 * 365 Rignot et al. 2016  [(0.6*365) (0.6*365) deep_melt/2 (0.6*365)]; 
melt_front = repmat(melt_year,md.mesh.numberofvertices,number_of_years);
melt_base = repmat(melt_year,1,number_of_years);
melt_time = [];

for i = 2024:(2024+number_of_years-1)
    melt_time = [melt_time ([0 (1/12)*5 (1/12)*8 (1/12)*11]+ i)];
end    

melt_front = [melt_front; melt_time];
md.frontalforcings.meltingrate=melt_front; %only effective if front is grounded



%Calving options
md.calving=calvingvonmises(); %activate von mises calving law
md.calving.stress_threshold_floatingice=400e3;
md.calving.stress_threshold_groundedice=calving_threshold;
md.calving.min_thickness=50; %m, default NaN



%Timestepping options
md.timestepping.cycle_forcing = 1;
md.timestepping = timestepping();
md.timestepping.time_step = (1/48); % Adaptive test suggested 0.02
md.settings.output_frequency = 4; %montlhy
% 	md.settings.output_frequency=1; %1: every tstep; 5: every fifth tstep, etc (for debugging)
% disp(['Setting fixed time step to ' num2str(md.timestepping.time_step) ' yrs'])
% md.timestepping=timesteppingadaptive();
% md.timestepping.time_step_max=5;
% md.timestepping.time_step_min=0.001;
md.timestepping.start_time = 2024;
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

%md = thinmodel(md,2024:1/12:final_year);

%savemodel(org,md);
save(['Outputs/' run_name],'md','-v7.3')

plotmodel(md, 'data', md.results.TransientSolution(1).Vel, 'data', md.results.TransientSolution(end).Vel, ...
    'mask', md.results.TransientSolution(1).MaskIceLevelset<0, 'mask#2-5', md.results.TransientSolution(end).MaskIceLevelset<0, ...
    'data', md.results.TransientSolution(end).Vel-md.results.TransientSolution(1).Vel,...
    'data', md.results.TransientSolution(end).Thickness-md.results.TransientSolution(1).Thickness,...
    'data', mean_SMB, 'caxis#3-4', [-250 250], 'ncols', 5,'caxis#1-2', [0 max(md.results.TransientSolution(end).Vel)], ...
    'title','Initial Velocity (m/yr)' , 'title','Final Velocity (m/yr)' , 'title', 'End velocity - Starting velocity (m/yr)', ...
    'title', 'End thickness - Starting thickness (m)', 'title', 'Mean SMB (mm/yr ice eq.)')
