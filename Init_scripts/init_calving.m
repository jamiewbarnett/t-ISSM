steps = [6];

% Todo:
% Spun up .mat files
% New shortened runme
% Questions/specific projects
% Check timesteps
% Processing Scripts into lecture
% Check with old binaries

%% %%%%%%%%%%%%% Glacier Selection %%%%%%%%%%%%%%

% Type the glacier you want to model below

glacier = '79'; %'79', 'Helheim', 'Kangerlussuaq' etc...

% Find correct exp and flowline files
switch glacier
    case{'79'} %Jamie
        exp_file = './Exp/79.exp';
        flowline_file = './Exp/79_flowline.exp';
        hmin = 750;
        hmax = 20000;
        fjordmesh = 1000;
        seasonalmelt = 0;
        deep_melt = 45;
        deep_depth = -600;
        upper_melt = 0;
        upper_depth = -100;
        sigma_floating = 200e3;
        eigen = 2e8;
        waterH = 0;
        minT = 90;
        icelandspc = 0;

    case{'Petermann'}%Felis
        exp_file = './Exp/petermann.exp';
        flowline_file = './Exp/petermann_flowline.exp';
        hmin = 750;
        hmax = 10000;
        fjordmesh = 750;
        sigma_floating = 200e3;
        eigen = 60e8;
        waterH = 0;
        minT = 100;
        deep_melt = 55;
        deep_depth = -700;
        upper_melt = 0;
        upper_depth = -75;
        icelandspc = 0;
end


parameterize_file = './Greenland.par';

%% %%%%%%%%%%%%% Transient Toggles %%%%%%%%%%%%%%%%
% Parameters to play with for the transient simulations (step 5)

%Transient
final_year = 2100; % Start year is 2024 and max possible final year 2100


%%%% SMB %%%%
smb_scenario = ['ssp245']; 

%%%% Calving %%%%%
calving_law = ['MT'];
 
%%%% Submarine Melt %%%%
melt_transient = []; %m/yr
melt_transient_time = [];
 
%%%% Calving %%%%
grounded_transient_sigmaMax =  [];
grounded_transient_time = [];% Times to apply the change in sigma max
floating_transient_sigmaMax = [];
floating_transient_time =  [];% Times to apply the change in sigma max

grounded_transient_sigmaMax = [];
grounded_transient_time = [];% Times to apply the change in sigma max
floating_transient_sigmaMax = [];
floating_transient_time = [];

%%%% Model name %%%%
ModelName = 'init'; %set your transient run name here
org = organizer('repository','Outputs','prefix',[glacier ModelName],'steps',steps);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%% Step 1: Meshing %%%%%%%%%%%%%%%
if perform(org,'Mesh')

    md = model;

    %Create initial mesh at resolution of 500m
    md=bamg(md,'domain',exp_file,'hmax',500);


    disp('Refining mesh to velocities')

    %Build grid & Velocity data
    x = double(ncread("Model_Data/vx.nc",'x'));
    y = double(ncread("Model_Data/vx.nc",'y'));
    vx	= double(ncread("Model_Data/vx.nc",'Band1'));
    vy	= double(ncread("Model_Data/vy.nc",'Band1'));
    vobs = sqrt(vx.^2+vy.^2);

    %Clip data to mesh domain
    posx  = find(x<=max(md.mesh.x) & x>=min(md.mesh.x));
    x_dom = x(posx);
    posy  = find(y>min(md.mesh.y) & y<=max(md.mesh.y));
    y_dom = (y(posy));

    vx_dom  = (vx(posx,posy)');
    vy_dom  = (vy(posx,posy)');
    vel_dom = sqrt(vx_dom.^2+vy_dom.^2);


    %Interpolate Velocities to Mesh
    md.inversion.vx_obs = InterpFromGridToMesh(x_dom,y_dom,(vx_dom),md.mesh.x,md.mesh.y,0);
    md.inversion.vy_obs = InterpFromGridToMesh(x_dom,y_dom,(vy_dom),md.mesh.x,md.mesh.y,0);
    md.inversion.vel_obs  = sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2);  

    %Define refinement everywhere where vel is fast (> vel_threshold)
    hmaxVertices=NaN*ones(md.mesh.numberofvertices,1);
    in=find(md.inversion.vel_obs>500); %find elemens where observed vel is fast >500m/yr
    hmaxVertices(in)=hmin;

    disp('Refining fjord area')

    %Refine areas in the fjord (ocean) to hmin
    md.mask.ocean_levelset = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask','nearest','./Model_Data/BedMachineGreenland-2022-03-17.nc');
    in = find(md.mask.ocean_levelset == 0);
    hmaxVertices(in) = fjordmesh;
    
    %Remesh again, adding refined region
    md=bamg(md,'hmin',hmin,'hmax',hmax,'field',...	
	    md.inversion.vel_obs,...
	    'gradation',1.1,'err',25,'anisomax',1.,'hmaxVertices',hmaxVertices);

    %Re-interpolate Velocities and Plot mesh
    md.inversion.vx_obs = InterpFromGridToMesh(x_dom,y_dom,(vx_dom),md.mesh.x,md.mesh.y,0);
    md.inversion.vy_obs = InterpFromGridToMesh(x_dom,y_dom,(vy_dom),md.mesh.x,md.mesh.y,0);
    md.inversion.vel_obs  = sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2);  
    plotmodel(md,...
        'data',md.inversion.vel_obs,'title','Velocity', ...
        'data','mesh', 'title', 'Mesh')

    %savemodel(org,md);
    save(['Outputs/' char(glacier) '_Mesh'],'md','-v7.3'); %char(org.steps(end).string)

end

%% %%%%%%%%%%%%% Step 2: Parameterize %%%%%%%%%%%%
if perform(org,'Parameterization')

    %md = loadmodel(org,'Mesh');
    load(['Outputs/' char(glacier) '_Mesh'])

    md = parameterize(md,parameterize_file);

    plotmodel(md,'data',md.mask.ocean_levelset, 'title', 'Ice/Ocean Mask', ...
        'data', md.geometry.thickness, 'title', 'Thickness (m)', ...
        'data', md.geometry.bed, 'title', 'Bedrock (m)', ...
        'data', 'BC','caxis#1',[-1,1])

    %savemodel(org,md);
    save(['Outputs/' char(glacier) '_Parameterization'],'md','-v7.3');

end

%% %%%%%%%%%%%%% Step 3: Rheology %%%%%%%%%%%

if perform(org,'Rheology_Inv')
    load(['Outputs/' char(glacier) '_Parameterization'])
    
    disp('Inverting for rheology_B')
    %Activate m1qn3-type inversion
    md.inversion=m1qn3inversion();
    
    %Set levelset options
    md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
    md.levelset.stabilization=1;
        %make it a signed distance
	md.mask.ice_levelset = reinitializelevelset(md,md.mask.ice_levelset);

    %Get mask from BedMachine
    M = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask','nearest','./Model_Data/BedMachineGreenland-2022-03-17.nc');
    pos=find(M<2);
    md.mask.ice_levelset(pos)=+1; %set levelset for no-ice vertices
    %remove 0 in ice_levelset (advection will fail if used)
    md.mask.ice_levelset(find(md.mask.ice_levelset==0))=-1;

    switch glacier
        case{'Petermann'}%        

            %Cost functions
            md.inversion.cost_functions=[101 103 502];
            md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
            md.inversion.cost_functions_coefficients(:,1)=1000; 
            md.inversion.cost_functions_coefficients(:,2)=1;  
            md.inversion.cost_functions_coefficients(:,3)=1.e-18; 

%             pos = find(md.mesh.x>-2.86e5 & md.mesh.x<-2.82e5 & md.mesh.y>-9.65e5 & md.mesh.y<-9.3e5);
%             md.stressbalance.spcvx(pos) = md.initialization.vx(pos);
%             md.stressbalance.spcvy(pos) = md.initialization.vy(pos);

%             in = ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,'Exp/Petermann_west_tongue.exp','node',1);
%             pos = find(in == 1 & md.geometry.bed > 0);
%             md.mask.ice_levelset(pos) = 1;
% 
%             %Fix the side wall of ice tongue...
%             in = ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,'Exp/HenningPetermannDomainWestFjordWall.exp','node',1);
%             pos = find(in == 1 & md.mask.ice_levelset<0);
% %             md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);
%             md.stressbalance.spcvx(pos) = md.initialization.vx(pos);
%             md.stressbalance.spcvy(pos) = md.initialization.vy(pos);

            %Set cost function weight to 0 by the west fjord
            %md.inversion.cost_functions_coefficients(pos,:) = 0;

        case{'79','79Z'}

            %Cost functions
            md.inversion.cost_functions=[101 103 502];
            md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
            md.inversion.cost_functions_coefficients(:,1)=1000; 
            md.inversion.cost_functions_coefficients(:,2)=1;  
            md.inversion.cost_functions_coefficients(:,3)=1.e-18; 


    end 

    %Set BCs along domain boundaries
    md.masstransport.spcthickness = NaN(md.mesh.numberofvertices,1);
    idxB=find(md.mesh.vertexonboundary==1); %vertices on boundary
    md.masstransport.spcthickness(idxB)=md.geometry.thickness(idxB);
    
    %Find vertices on boundary and ice
    pos=find(md.mesh.vertexonboundary & M>=2);
    
    %Set Dirichlet
    md.stressbalance.spcvx(pos)=md.initialization.vx(pos); % west/east flow
    md.stressbalance.spcvy(pos)=md.initialization.vy(pos); % north/south flow

    %Set Dirichlet for non ice vertices
    pos = find(M==1 | ~isnan(md.stressbalance.spcvx));
    md.stressbalance.spcvx(pos)=md.initialization.vx(pos);
    md.stressbalance.spcvy(pos)=md.initialization.vy(pos);

    %set flowlaw
    md = setflowequation(md,'SSA','all');

    %Get velocities for the inversion
    md.inversion.vx_obs=md.initialization.vx;
    md.inversion.vy_obs=md.initialization.vy;
    md.inversion.vz_obs=md.initialization.vz;
    md.inversion.vel_obs=md.initialization.vel;

    %Control general
    md.inversion.iscontrol=1;
    md.inversion.maxsteps=50;
    md.inversion.maxiter=50;
    md.inversion.dxmin=1.0e-6; %default 0.1
    md.inversion.gttol=1.0e-6;
    md.verbose=verbose('control',true);


    %Controls
    md.inversion.control_parameters={'MaterialsRheologyBbar'};
    md.inversion.min_parameters=md.materials.rheology_B;
    md.inversion.max_parameters=md.materials.rheology_B;
    
    %Set min and max parameters for floating ice
    pos=find(md.mask.ocean_levelset<0);
    md.inversion.min_parameters(pos)=5e6;
    md.inversion.max_parameters(pos)=cuffey(200);
    
    %Additional parameters
    md.stressbalance.restol=0.01;
    md.stressbalance.reltol=0.1;
    md.stressbalance.abstol=NaN;

    %Fill in blanks in velocity data
    pos = find(isnan(md.inversion.vx_obs) | isnan(md.inversion.vy_obs) | md.inversion.vel_obs == 0);
    md.inversion.vx_obs(pos) = 0;
    md.inversion.vy_obs(pos) = 0;
    md.inversion.vel_obs(pos) = 0;
    md.inversion.cost_functions_coefficients(pos,[1 2 3]) = 0;

    pos = find(md.inversion.vel_obs>0 & md.mask.ice_levelset>0);
    md.inversion.vx_obs(pos) = 0;
    md.inversion.vy_obs(pos) = 0;
    md.inversion.vel_obs(pos) = 0;
    md.inversion.cost_functions_coefficients(pos,[1 2 3]) = 0;
    

    % Solve
    md.cluster=generic('name',oshostname,'np',4);


    switch glacier
        case{'Petermann'}
            mds=extract(md,md.mask.ocean_levelset<0); %only do inversion for floating ice

             %Cost function set to 0 in areas with no velocity
             pos = find(mds.inversion.vel_obs==0);
             mds.inversion.cost_functions_coefficients(pos,[1 2]) = 0;
             % Weigh abs velocity as more important
             mds.inversion.cost_functions_coefficients(:,1) = 1*mds.inversion.cost_functions_coefficients(:,1);
             % Reduce minimum allowed rheology_B to better capture shear margins
             mds.inversion.min_parameters(:)=0.1*cuffey(273.15);
             %Retreat ice levelset by 1, to prevent artifacts at the ice front
             pos = find(max(mds.mask.ice_levelset(mds.mesh.elements),[],2)>0);
             mds.mask.ice_levelset(mds.mesh.elements(pos,:)) = +1;

        case{'Ryder'}
            mds=extract(md,md.mask.ocean_levelset<0); %only do inversion for floating ice
             %Cost function set to 0 in areas with no velocity
             pos = find(mds.inversion.vel_obs==0);
             mds.inversion.cost_functions_coefficients(pos,[1 2]) = 0;
             % Weigh abs velocity as more important
             mds.inversion.cost_functions_coefficients(:,1) = 5*mds.inversion.cost_functions_coefficients(:,1);
             % Reduce minimum allowed rheology_B to better capture shear margins
             mds.inversion.min_parameters(:)=5e6;
             %Retreat ice levelset by 1, to prevent artifacts at the ice front
             pos = find(max(mds.mask.ice_levelset(mds.mesh.elements),[],2)>0);
             mds.mask.ice_levelset(mds.mesh.elements(pos,:)) = +1;

        case{'Storstrommen'}
            in = ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,'Exp/storstrommen_front.exp','node',1);
%           pos = find(in == 0 & mds.mask.ice_levelset<0);
%             mds.inversion.cost_functions_coefficients(pos,[1 2 3]) = 0;
            mds=extract(md,(md.mask.ocean_levelset<0 & in == 1)); %only do inversion for floating ice

        case{'79','79Z'}
            mds=extract(md,md.mask.ocean_levelset<0); %only do inversion for floating ice
            
            %Cost function set to 0 in areas with no velocity
            pos = find(mds.inversion.vel_obs==0);
            mds.inversion.cost_functions_coefficients(pos,[1 2]) = 0;
            % Weigh abs velocity as more important
            mds.inversion.cost_functions_coefficients(:,1) = 5*mds.inversion.cost_functions_coefficients(:,1);
            % Reduce minimum allowed rheology_B to better capture shear margins
            mds.inversion.min_parameters(:)=5e6;
            %Retreat ice levelset by 1, to prevent artifacts at the ice front
            pos = find(max(mds.mask.ice_levelset(mds.mesh.elements),[],2)>0);
            mds.mask.ice_levelset(mds.mesh.elements(pos,:)) = +1;

        case{'Steensby'}
            mds=extract(md,md.mask.ocean_levelset<0); %only do inversion for floating ice

             %Cost function set to 0 in areas with no velocity
             pos = find(mds.inversion.vel_obs==0);
             mds.inversion.cost_functions_coefficients(pos,[1 2]) = 0;
             % Weigh abs velocity as more important
             mds.inversion.cost_functions_coefficients(:,1) = 100*mds.inversion.cost_functions_coefficients(:,1);
             % Reduce minimum allowed rheology_B to better capture shear margins
             mds.inversion.min_parameters(:)=5e6;
             %Retreat ice levelset by 1, to prevent artifacts at the ice front
             pos = find(max(mds.mask.ice_levelset(mds.mesh.elements),[],2)>0);
             mds.mask.ice_levelset(mds.mesh.elements(pos,:)) = +1;
    end

    mds.settings.solver_residue_threshold = 1e-04;

    mds=solve(mds,'Stressbalance');


    % Update model rheology_B accordingly
    md.materials.rheology_B(mds.mesh.extractedvertices)=mds.results.StressbalanceSolution.MaterialsRheologyBbar;
    
    
    %md.materials.rheology_B = md.results.StressbalanceSolution.MaterialsRheologyBbar;
    plotmodel(mds,'data',mds.results.StressbalanceSolution.MaterialsRheologyBbar,'data',mds.results.StressbalanceSolution.Vel,'data',mds.inversion.vel_obs,'data',mds.results.StressbalanceSolution.Vel - mds.initialization.vel,'ncols',4);


%     md.mask.ice_levelset = reinitializelevelset(md,md.mask.ice_levelset);
%     md.mask.ocean_levelset = reinitializelevelset(md,md.mask.ocean_levelset);
% 
% 
    flags = zeros(md.mesh.numberofvertices,1);
	flags(md.mask.ocean_levelset<0 & md.mask.ice_levelset>-3000 & md.mask.ice_levelset<0 & md.materials.rheology_B<1.4e8) = 1;
	pos1 = find(flags);
    flags = zeros(md.mesh.numberofvertices,1);
    flags(md.mask.ocean_levelset<0 & md.mask.ocean_levelset>-3000 & md.mask.ice_levelset>0) = 1;
    pos2 = find(flags);
	md.materials.rheology_B(pos2) = griddata(md.mesh.x(pos1),md.mesh.y(pos1), md.materials.rheology_B(pos1), md.mesh.x(pos2),md.mesh.y(pos2), 'nearest');

    save(['Outputs/' char(glacier) '_Rheology'],'md','-v7.3');

end


%% %%%%%%%%%%%%% Step 4: Stressbalance %%%%%%%%%%%
if perform(org,'Stressbalance')

    %md = loadmodel(org,'Parameterization');
    load(['Outputs/' char(glacier) '_Rheology'])
    
    md = setflowequation(md,'SSA','all');

    %Set levelset options
    md.levelset.stabilization=1;

    %Get mask from BedMachine
    M = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask','nearest','./Model_Data/BedMachineGreenland-2022-03-17.nc');
    pos=find(M<2);
    md.mask.ice_levelset(pos)=+1; %set levelset for no-ice vertices
    %remove 0 in ice_levelset (advection will fail if used)
    md.mask.ice_levelset(find(md.mask.ice_levelset==0))=-1;
    
    %make it a signed distance
	md.mask.ice_levelset = reinitializelevelset(md,md.mask.ice_levelset);

    %Set spcthickness to NaN
    md.masstransport.spcthickness=NaN(md.mesh.numberofvertices,1);
    
    %Set Ice Boundary conditions
    pos=find(md.mesh.vertexonboundary & M>=2);
    md.stressbalance.spcvx(pos)= md.inversion.vx_obs(pos);
    md.stressbalance.spcvy(pos)= md.inversion.vy_obs(pos);
    md.masstransport.spcthickness(pos) = md.geometry.thickness(pos);

    %Set Dirichlet for Non-ice vertices
    pos = find(M==1 | ~isnan(md.stressbalance.spcvx));
    md.stressbalance.spcvx(pos)=md.inversion.vx_obs(pos);
    md.stressbalance.spcvy(pos)=md.inversion.vy_obs(pos);
   
    %%% PERFORM INVERSION %%%
    disp('Inverting for Friction')

    %Activate m1qn3-type inversion
    md.inversion=m1qn3inversion();

    %Get velocities for the inversion
    md.inversion.vx_obs=md.initialization.vx;
    md.inversion.vy_obs=md.initialization.vy;
    md.inversion.vz_obs=md.initialization.vz;
    md.inversion.vel_obs=md.initialization.vel;

    %Control general
    md.inversion.iscontrol=1;
    md.inversion.maxsteps=80;
    md.inversion.maxiter=80;

    %Cost functions
    md.inversion.cost_functions=[101 103 501];
    md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
    md.inversion.cost_functions_coefficients(:,1)=1000; 
    md.inversion.cost_functions_coefficients(:,2)=1;  
    md.inversion.cost_functions_coefficients(:,3)=1e-17; 

    %Controls
    md.inversion.control_parameters={'FrictionCoefficient'};
    md.inversion.min_parameters=10*ones(md.mesh.numberofvertices,1); %default 10
    md.inversion.max_parameters=200*ones(md.mesh.numberofvertices,1);

    %Additional parameters
    md.stressbalance.restol=0.01;
    md.stressbalance.reltol=0.1;
    md.stressbalance.abstol=NaN;

    %%% INVERSION END %%%
   

    %No Damage model
    md.damage.D = zeros(md.mesh.numberofvertices,1);
    md.damage.spcdamage = NaN*ones(md.mesh.numberofvertices,1);

    %Dont use thermal model
    md.thermal.spctemperature=NaN*ones(md.mesh.numberofvertices,1);

    %    md.toolkits=toolkits;
    md.toolkits = toolkits();
    md.verbose=verbose('solution',false,'module',false,'convergence',false);


    md.cluster=generic('name',oshostname,'np',4);
	%Solve for inverted friction
	md=solve(md,'sb','np',4);

    %Update model friction from stressbalance solution
    md.friction.coefficient = md.results.StressbalanceSolution.FrictionCoefficient;


    %Set friction for non-ice regions
    pos = find(md.mask.ocean_levelset<0 | md.mask.ice_levelset>1 | (md.friction.coefficient>150 & md.mask.ocean_levelset<1000 & md.geometry.bed <0));
    md.friction.coefficient(pos)=  120*((min(max(0,md.geometry.bed(pos)+800),max(md.geometry.bed))/max(md.geometry.bed))); %Akesson (2018)

    md.friction.p = ones(md.mesh.numberofelements,1);
    md.friction.q = ones(md.mesh.numberofelements,1);


    %Make sure base is not below bed
	pos1=find(md.geometry.bed>md.geometry.base);
	md.geometry.base(pos1)=md.geometry.bed(pos1);

	%Recalculate surface
	md.geometry.surface=md.geometry.base+md.geometry.thickness;

	%Make sure bed = base for grounded ice
	pos1=find(md.mask.ocean_levelset>0);
	md.geometry.base(pos1)=md.geometry.bed(pos1);

	%Recalculate surface
	md.geometry.surface=md.geometry.base+md.geometry.thickness;

    M = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask','nearest','./Model_Data/BedMachineGreenland-2022-03-17.nc');

    %Find vertices on boundary and within exp files
    pos=find(md.mesh.vertexonboundary & M>=1);
    
    %Set Dirichlet
    md.stressbalance.spcvx(pos)= md.inversion.vx_obs(pos);
    md.stressbalance.spcvy(pos)= md.inversion.vy_obs(pos);

    %REMOVE dirichlet conditions on non-ice vertices
%     pos = find(M==1 & md.mesh.vertexonboundary==0);
%     md.stressbalance.spcvx(pos)=NaN;
%     md.stressbalance.spcvy(pos)=NaN;

	%Test plot
	plotmodel(md,...
 		'data',md.inversion.vel_obs,'title','Observed velocity (m/yr)',...
 		'data',md.results.StressbalanceSolution.Vel,'title','Modeled Velocity (m/yr)',...
 		'data',md.geometry.bed,'title','Bed elevation (m)',...
 		'data',md.friction.coefficient,'title','Friction Coefficient',...
 		'colorbar#all','on',...
 		'caxis#1-2',([1.5,round(max(md.inversion.vel_obs),-2)]),...
 		'log#1-2',10,'ncols',4);

    %savemodel(org,md);
    save(['Outputs/' char(glacier) '_Stressbalance'],'md','-v7.3');

end

%% %%%%%%%%%%%%% Step 5: Spin_Up Melt %%%%%%%%%%%%%
if perform(org,'Spin_Up_Melt')

    load(['Outputs/' char(glacier) '_Stressbalance'])

    %Intial Velocities
    md.initialization.vx = md.results.StressbalanceSolution.Vx;
    md.initialization.vy = md.results.StressbalanceSolution.Vy;
    md.initialization.vel = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2);

    %Dont use damage model
	md.damage.D=zeros(md.mesh.numberofvertices,1);
	md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);

	%Dont use thermal model
	md.transient.isthermal=0;
	md.thermal.spctemperature=NaN*ones(md.mesh.numberofvertices,1);

	%Additional options
	md.inversion.iscontrol=0;
    
    disp('Reading and interpolating SMB data')
    
    md.smb.mass_balance = [];
    smbMAR = [];

    for yy=(648-(10*12)):648 % Monthly data 2004 2014
        smboutput = interpMAR_monthly(md.mesh.x,md.mesh.y,'SMB',yy, './Model_Data/MARv3.11.3-historical-combined.nc');
        smbMAR = [smbMAR smboutput];
        %progress = sprintf('Read %d timesteps out of %d',yy, nyrs_smb*12);
        %disp(progress)
    end
   
    pos=find(smbMAR==-9999);
    smbMAR(pos)=0.0;
          
    md.smb.mass_balance = smbMAR/1000*12*(md.materials.rho_freshwater/md.materials.rho_ice); 
    md.smb.mass_balance = mean(md.smb.mass_balance,2);


    %Make sure bed is below base
    pos=find(md.geometry.bed>md.geometry.base);
    md.geometry.base(pos)=md.geometry.bed(pos);

    %Recalculate surface
    md.geometry.surface=md.geometry.base+md.geometry.thickness;
    
    %Front and GL options
    md.transient.ismovingfront=0;
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

    if icelandspc == 1
        md.levelset.spclevelset=NaN(md.mesh.numberofvertices, 1);
        pos = find_iceLandBoundary(md, 1); %1=is2D
        md.levelset.spclevelset(pos)=-1;
    end


    %Basal Melt options
    % Fixed melt
    md.basalforcings=linearbasalforcings();
	%md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.deepwater_melting_rate = deep_melt;
    md.basalforcings.deepwater_elevation = deep_depth;
    md.basalforcings.upperwater_melting_rate = upper_melt;
    md.basalforcings.upperwater_elevation = upper_depth;
    md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);
	md.basalforcings.geothermalflux=interpSeaRISE_new(md.mesh.x,md.mesh.y,'bheatflx');


     md.levelset.kill_icebergs=1;
    %md.levelset.migration_max=10000; % -- maximum allowed migration rate (m/a)

    %Timestepping options
    md.timestepping = timesteppingadaptive;
    md.timestepping.time_step_min = (1/365)*0.72; 
    md.timestepping.time_step_max = (1/365)*7.2; 
    md.settings.output_frequency = 25;
    md.timestepping.start_time = 0;
    md.timestepping.final_time= 10;

    %Output options
    md.transient.requested_outputs={'TotalSmb','SmbMassBalance',...
        'IceVolume','BasalforcingsFloatingiceMeltingRate','MaskIceLevelset','MaskOceanLevelset'};
    
    md.verbose=verbose('solution',true,'module',false,'convergence',false);

    md.toolkits=toolkits();
    md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
	md.cluster=generic('name',oshostname,'np',4);
    md.settings.solver_residue_threshold=1e-4;

	%Solve
	md=solve(md,'Transient');

    %savemodel(org,md);
    save(['Outputs/' char(glacier) '_Spinup_Melt'],'md','-v7.3');

    plotmodel(md, 'data', md.inversion.vel_obs, 'data', md.results.TransientSolution(end).Vel, ...
        'mask', md.mask.ice_levelset<0, 'mask#2-4', md.results.TransientSolution(end).MaskIceLevelset<0, ...
        'caxis#1-2', [0 max(max(md.inversion.vel_obs),max(md.results.TransientSolution(end).Vel))], 'data', md.results.TransientSolution(end).Thickness-md.results.TransientSolution(1).Thickness,...
        'data', md.results.TransientSolution(end).MaskOceanLevelset, 'caxis#3', [-250 250] , 'caxis#4', [-1 1], 'ncols', 4,...
        'title','Observed Velocity (m/yr)' , 'title','Modelled Velocity (m/yr)' , 'title', 'End thickness - Starting thickness (m)' , 'title', 'Ocean mask' )
       

end

%% %%%%%%%%%%%%% Step 6: Spin Up Calving %%%%%%%%%%%%%

if perform(org,'Spin_Up_Calving')

    load(['Outputs/' char(glacier) '_Spinup_Melt'])

    %set some things
    md.miscellaneous.name = [char(glacier) '_' char(ModelName)];

    %Intial Velocities
    md = transientrestart(md);
    
    %Front and GL options
    md.transient.ismovingfront=1;
    md.transient.isgroundingline=1;
    md.groundingline.migration='SubelementMigration';
    md.groundingline.melt_interpolation='NoMeltOnPartiallyFloating';
    md.groundingline.friction_interpolation='SubelementFriction1';	

  
    md.levelset.kill_icebergs=1;
    %md.levelset.migration_max=10000; % -- maximum allowed migration rate (m/a)

    switch calving_law
        case{'VM'}
            md.calving=calvingvonmises();
            md.calving.stress_threshold_floatingice = sigma_floating;
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

    %Set frontal melt to 0
    md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);    

   %Timestepping options
    md.timestepping = timesteppingadaptive;
    md.timestepping.time_step_min = (1/365)*0.72; 
    md.timestepping.time_step_max = (1/365)*7.2; 
    md.settings.output_frequency = 25;
    md.timestepping.start_time = 0;
    md.timestepping.final_time = 25;

    %Output options
    md.transient.requested_outputs={'TotalSmb','SmbMassBalance',...
        'IceVolume','BasalforcingsFloatingiceMeltingRate','MaskIceLevelset','MaskOceanLevelset'};
    
    md.verbose=verbose('solution',true,'module',false,'convergence',false);

    md.toolkits=toolkits();
    md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
	md.cluster=generic('name',oshostname,'np',4);
    md.settings.solver_residue_threshold=1e-4;

	%Solve
	md=solve(md,'Transient');

    %savemodel(org,md);
    save(['Outputs/' char(glacier) '_Spinup_Calving_' char(calving_law)],'md','-v7.3');

    plotmodel(md, 'data', md.inversion.vel_obs, 'data', md.results.TransientSolution(end).Vel, ...
        'mask', md.mask.ice_levelset<0, 'mask#2-4', md.results.TransientSolution(end).MaskIceLevelset<0, ...
        'caxis#1-2', [0 max(max(md.inversion.vel_obs),max(md.results.TransientSolution(end).Vel))], 'data', md.results.TransientSolution(end).Thickness-md.results.TransientSolution(1).Thickness,...
        'data', md.results.TransientSolution(end).MaskOceanLevelset, 'caxis#3', [-250 250] , 'caxis#4', [-1 1], 'ncols', 4,...
        'title','Observed Velocity (m/yr)' , 'title','Modelled Velocity (m/yr)' , 'title', 'End thickness - Starting thickness (m)' , 'title', 'Ocean mask' )
end
















