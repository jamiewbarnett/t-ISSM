%% Project Description
% Here you can investigate the future of Petermann Glacier.
% The aim of the projection to to explore the future of the glacier with
% respect to the choice of calving law. You can choose Von Mises (VM),
% Eigen calving (EC), Crevasse Depth calving (CD) and a minimum thickness
% law (MT).

% In the toggles section below you can choose the calving law, the
% future SMB scenario and, if you're feeling adventurous, the ocean melt at
% the grounding line.

% I would begin by running short simulations, a few decades in length, and
% exploring the data analysis scripts to visualse the glacier's changes over
% time.

% If you feel crazy you can also change the calirated tuning parameters in
% the calving law. Read Choi et al (2018)/Wilner et al (2023) or ask for
% help!

% Possible questions to answer:
% When will the ice shelf disintegrate? Does this vary with calving law?
% What is the effect on mass loss if the shelf disintegrates?
% Is the glacier/shelf more sensitive to SMB or ocean melt changes? 

% Relevant Papers:
% Millan et al (2023) - https://www.nature.com/articles/s41467-023-42198-2
% Hill et al (2021) - https://doi.org/10.1017/jog.2020.97
% Choi et al (2018) - https://doi.org/10.5194/tc-12-3735-2018
% Wilner et al (2023) - https://doi.org/10.5194/tc-17-4889-2023

flowline = 'Exp/Petermann_flowline.exp';

%% %%%%%%%%%%%%% Toggles %%%%%%%%%%%%%%%%
% Parameters to play with for the transient simulations

%Transient
final_year = 2100; % Start year is 2025 

%%% Calving Law %%% 
calving_law = 'CD'; % VM | EIGEN | CD | MT

%%%% SMB %%%%
smb_scenario = 'ssp585'; % Low Emissions: ssp245 | High Emissions: ssp585

%%%% Melt at the Grounding Line %%%%
GL_melt = 55; % Default = 55 m/yr

%%%% Calving Calibrations %%%%
sigma_floating = 200e3; %Deafult 200e3
eigen = 60e8; %Default 60e8
waterH = 0; %Default 0
minT = 100; %Default 100

run_name = 'ice_go_byebye'; %Important to name your run well. 

%% %%%%%%%%%%%%%% Code to Run Simulation %%%%%%%%%%%%%%%

load(['Model_Data/Petermann_Spinup_Calving_' calving_law '.mat']); %Load spun-up model 

md = transientrestart(md); %Initialise model from spin up simulations

% Fix future SMB
disp('Reading and interpolating SMB data')

md.smb.mass_balance = [];
smbMAR = [];

switch smb_scenario
    case{'ssp245'} 
        %smb_file='./Model_Data/MARv3.11.3-ssp245-combined.nc';
        load('Model_Data/Petermann_SMB_ssp245.mat');
    case{'ssp585'}
        %smb_file='./Model_Data/MARv3.11.3-ssp585-combined.nc';
        load('Model_Data/Petermann_SMB_ssp585.mat');
end

md.smb.mass_balance = smb_readin;
% 
% smb_startyear = 2024; %Year to begin SMB series. Could be changed...
% 
% for yy=((smb_startyear-2014)*12)+1:((final_year-2014))*12 % Monthly data
%     smboutput = interpMAR_monthly(md.mesh.x,md.mesh.y,'SMB',yy, smb_file);
%     smbMAR = [smbMAR smboutput];
% end
% 
% pos=find(smbMAR==-9999);
% smbMAR(pos)=0.0;
%       
% md.smb.mass_balance = [smbMAR/1000*12*(md.materials.rho_freshwater/md.materials.rho_ice); ...
%     [2024:1/12:final_year-(1/12)]]; % Monthly transient forcing
% mean_SMB = mean(smbMAR/1000*12*(md.materials.rho_freshwater/md.materials.rho_ice), 2);

% Calibrated calving laws...
switch calving_law
    case{'VM'}
        md.calving=calvingvonmises();
        md.calving.stress_threshold_floatingice = sigma_floating;
        md.calving.min_thickness = 20;
    case{'EIGEN'}
        md.calving=calvinglevermann();
        md.calving.coeff = eigen*md.constants.yts *ones(md.mesh.numberofvertices,1);
    case{'CD'}
        md.calving=calvingcrevassedepth();
        md.calving.crevasse_opening_stress = 0;
        md.calving.water_height = waterH*ones(md.mesh.numberofvertices,1);
    case{'MT'}
        md.calving=calvingminthickness();
        md.calving.min_thickness = minT;
end

%Frontal melt, only works if glacier is grounded
md.frontalforcings.meltingrate = (GL_melt/2)*zeros(md.mesh.numberofvertices,1);   

% Set Grounding Line Melt
md.basalforcings.deepwater_melting_rate = GL_melt;

%Set intial and final timestep
md.timestepping = timestepping();
md.timestepping.start_time = 2025;
md.timestepping.final_time = final_year; 
md.timestepping.time_step = 1/48;
md.settings.output_frequency = 4; %Monthly time steps
md.timestepping.cycle_forcing = 1;

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
    'TotalCalvingFluxLevelset','TotalGroundedBmb'};

md.verbose=verbose('solution',true,'module',false,'convergence',false);

md.toolkits=toolkits();
md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
md.cluster=generic('name',oshostname,'np',4);

%Solve
md=solve(md,'Transient');

%md = thinmodel(md,2025:1/12:final_year); %Produce a montly output
md.miscellaneous.name = run_name;
save(['Outputs/Petermann_' run_name],'md','-v7.3');









