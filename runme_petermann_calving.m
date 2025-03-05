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

% Possible questions to answer:
% When will the ice shelf disintegrate? Does the vary with calving law?
% What is the effect on mass loss if the shelf disintegrates?
% Is the glacier/shelf more sensitive to SMB or ocean melt changes?

% Relevant Papers:

%% %%%%%%%%%%%%% Toggles %%%%%%%%%%%%%%%%
% Parameters to play with for the transient simulations (step 5)

%Transient
final_year = 2100; % Start year is 2025 

%%% Calving Law %%% 
calving_law = 'VM'; % VM | EIGEN | CD | MT

%%%% SMB %%%%
smb_scenario = 'ssp585'; % Low Emissions: ssp245 | High Emissions: ssp585

%%%% Melt at the Grounding Line %%%%
GL_melt = 55; % Default = 55 m/yr

run_name = 'ice_go_byebye'; %Important to name your run well. 

%% %%%%%%%%%%%%%% Code to Run Simulation %%%%%%%%%%%%%%%

load(['Petermann_Spinup_Calving_' calving_law '.mat']); %Load spun-up model 

md = transientrestart(md); %Initialise model from spin up simulations

% Fix future SMB
disp('Reading and interpolating SMB data')

md.smb.mass_balance = [];
smbMAR = [];

switch smb_scenario %Select SMB scenario
    case{'ssp245'} 
        smb_file='./Model_Data/MARv3.11.3-ssp245-combined.nc';
    case{'ssp585'}
        smb_file='./Model_Data/MARv3.11.3-ssp585-combined.nc';
end

for yy=133:((final_year-2025)+11)*12 % Monthly data
    smboutput = interpMAR_monthly(md.mesh.x,md.mesh.y,'SMB',yy, smb_file);
    smbMAR = [smbMAR smboutput];
end

pos=find(smbMAR==-9999);
smbMAR(pos)=0.0;
      
md.smb.mass_balance = [smbMAR/1000*12*(md.materials.rho_freshwater/md.materials.rho_ice); ...
    [2025:1/12:final_year-(1/12)]]; % Monthly transient forcing
mean_SMB = mean(smbMAR/1000*12*(md.materials.rho_freshwater/md.materials.rho_ice), 2);

% Set Grounding Line Melt
md.basalforcings.deepwater_melting_rate = GL_melt;

%Set intial and final timestep
md.timestepping.start_time = 2025;
md.timestepping.final_time = final_year; 
md.settings.output_frequency = 10;


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

md = thinmodel(md,2025:1:final_year); %Produce a yearly output
save(['Outputs/Petermann_' run_name],'md','-v7.3');









