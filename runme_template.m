steps = [1:3];

%%%%%%%%%%%%% Make mesh domain/exp file %%%%%%%%%%%%%%

exp_file = './Exp/79.exp';
parameterize_file = './Greenland.par';

%% %%%%%%%%%%%%% Toggles and things %%%%%%%%%%%%%%

%Mesh
hmin = 1000;
hmax = 20000;


%Transient
nyrs = 1;

%Timestepping
timestep = 0.05;
outfreq = 1/timestep; % Annual output

%SMB Years
nyrs_smb = 2100-2099; % End and start year of dataset


%%%% Model name %%%%
ModelName = 'NEWTEST';
org = organizer('repository','Models','prefix',[ModelName num2str(nyrs) 'years'],'steps',steps);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
%% %%%%%%%%%%%%%%% Step 1: Meshing %%%%%%%%%%%%%%%

if perform(org,'Mesh')

    md = model;

    %Create initial mesh at resolution of 500m
    md=bamg(md,'domain',exp_file,'hmax',500);


    disp('Refining mesh to bedmachine')
    %Interpolate bedrock data from bedmachine onto mesh
    md.geometry.bed = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'bed','nearest','./Model_Data/BedMachineGreenland-2022-03-17.nc');
    

    hmaxVertices = NaN*ones(md.mesh.numberofvertices,1);
    in = ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,exp_file,'node',1);
    hmaxVertices(find(in)) = hmax;

    %Refine mesh with bedmachine data
    md = bamg(md,'hmin', hmin, 'hmax', hmax, ...
        'field', md.geometry.bed, 'gradation', 1.1, 'err', 75, ...
        'anisomax', 1., 'hmaxVertices', hmaxVertices);

    %Interpolate bedrock data again onto new mesh
    md.geometry.bed = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'bed','nearest','./Model_Data/BedMachineGreenland-2022-03-17.nc');

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
    in=find(md.inversion.vel_obs>300); %find elemens where observed vel is fast >300m/yr
    hmaxVertices(in)=hmin;

    disp('Refining fjord area')

    %Refine areas in the fjord (ocean) to hmin
    md.mask.ocean_levelset = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask','nearest','./Model_Data/BedMachineGreenland-2022-03-17.nc');
    in = find(md.mask.ocean_levelset == 0);
    hmaxVertices(in) = hmin;
    
    %Remesh again, adding refined region
    md=bamg(md,'hmin',hmin,'hmax',hmax,'field',...	
	    md.geometry.bed,...
	    'gradation',1.1,'err',50,'anisomax',1.,'hmaxVertices',hmaxVertices);

    %Re-interpolate Velocities and Plot mesh
    md.inversion.vx_obs = InterpFromGridToMesh(x_dom,y_dom,(vx_dom),md.mesh.x,md.mesh.y,0);
    md.inversion.vy_obs = InterpFromGridToMesh(x_dom,y_dom,(vy_dom),md.mesh.x,md.mesh.y,0);
    md.inversion.vel_obs  = sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2);  
    plotmodel(md,...
        'data',md.inversion.vel_obs,'title','Velocity', ...
        'data','mesh', 'title', 'Mesh')

    savemodel(org,md);

end



%%
%% %%%%%%%%%%%%% Step 2: Parameterize %%%%%%%%%%%%
if perform(org,'Parameterization')
    md = loadmodel(org,'Mesh');

    md = parameterize(md,parameterize_file);

%     plotmodel(md,'data',md.mask.ocean_levelset, 'title', 'Ice/Ocean Mask', ...
%         'data', md.geometry.thickness, 'title', 'Thickness', ...
%         'data', md.geometry.bed, 'title', 'Bedrock', ...
%         'data', 'BC','caxis#1',[-1,1])

    savemodel(org,md);

end
%%
%% %%%%%%%%%%%%% Step 3: Stressbalance %%%%%%%%%%%

if perform(org,'Stressbalance')
    md = loadmodel(org,'Parameterization');
    
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
    md.stressbalance.spcvx(pos)=0; %no west/east flow
    md.stressbalance.spcvy(pos)=0; %no north/south flow


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
    md.inversion.maxsteps=40;
    md.inversion.maxiter=40;

    %Cost functions
    md.inversion.cost_functions=[101 103 501];
    md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
    md.inversion.cost_functions_coefficients(:,1)=2000; 
    md.inversion.cost_functions_coefficients(:,2)=1;  
    md.inversion.cost_functions_coefficients(:,3)=6e-7; 

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
    pos = find(md.mask.ocean_levelset<0 | md.mask.ice_levelset>1);
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

	%Test plot
	plotmodel(md,...
 		'data',md.inversion.vel_obs,'title','Observed velocity',...
 		'data',md.results.StressbalanceSolution.Vel,'title','Modeled Velocity',...
 		'data',md.geometry.bed,'title','Bed elevation',...
 		'data',md.friction.coefficient,'title','Friction Coefficient',...
 		'colorbar#all','on','colorbartitle#1-2','(m/yr)',...
 		'caxis#1-2',([1.5,1500]),...
 		'colorbartitle#3','(m)', 'log#1-2',10);

    savemodel(org,md);
end

%%
%% %%%%%%%%%%%%%%% Step 4: Transient %%%%%%%%%%%%%

if perform(org,'Transient')

    md = loadmodel(org,'Stressbalance');

    if isrestart == 1
        %Initialize geometry and update mask from previous transient results
        md = transientrestart(md);
    end

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

    for yy=1:(nyrs_smb*12) % Monthly data
        smboutput = interpMAR_monthly(md.mesh.x,md.mesh.y,'SMB',yy, './Model_Data/MARv3.11.3-ssp585-combined.nc');
        smbMAR = [smbMAR smboutput];
        %progress = sprintf('Read %d timesteps out of %d',yy, nyrs_smb*12);
        %disp(progress)
    end
   
    for zz=1:length(smbMAR)
        for vv=1:width(smbMAR)
            if smbMAR(zz, vv) == -9999 % Remove instances of fill/missing value
                smbMAR(zz, vv) = 0.0;
            else
                smbMAR(zz, vv) = smbMAR(zz, vv);
            end
        end
    end
          
    md.smb.mass_balance = ...
           [smbMAR/1000*12*(md.materials.rho_freshwater/md.materials.rho_ice); ...
           [0:1/12:(nyrs_smb)-(1/12)]]; % 1/12 used with monthly input. Change both instances to e.g. 1 for yearly data 



    %Make sure we're using correct flow equation
    if isHO == 1
	    md = setflowequation(md,'HO','all');
    else
	    md = setflowequation(md,'SSA','all');
    end
    
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
        if isHO == 1
		    %do nothing
	    else
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
        end
    else 
	    %dont touch the spclevelset, just keep what is from the previous model and do nothing here
    end

    md.levelset.kill_icebergs=1;
    %md.levelset.migration_max=10000; % -- maximum allowed migration rate (m/a)

        %Calving options
    if md.transient.ismovingfront==1
	    md.calving=calvingvonmises(); %activate von mises calving law
    
	    %Stress threshold
	    md.calving.stress_threshold_groundedice=1e6; %default 1 MPa = 1e6 Pa
    
	    md.calving.stress_threshold_floatingice=300e3; %default Petermann 300 kPa, default ISSM 150 kPa
	    disp(['Calving sigma_max floatingice set to ' num2str(md.calving.stress_threshold_floatingice./1000)  ' kPa'])
    
	    md.calving.min_thickness=50; %m, default NaN
    
	    %Define calving rate and melt rate (only effective if ismovingfront==1)
	    md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices,1); %only effective if front is grounded
    end

    %Basal Melt options
    % Fixed melt
    md.basalforcings=basalforcings();
	%md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate=10*ones(md.mesh.numberofvertices,1);
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);


    %Timestepping options

    md.timestepping = timestepping();
    md.timestepping.time_step = timestep;
    md.settings.output_frequency = outfreq; %yearly
% 	md.settings.output_frequency=1; %1: every tstep; 5: every fifth tstep, etc (for debugging)
    disp(['Setting fixed time step to ' num2str(md.timestepping.time_step) ' yrs'])




    md.timestepping.final_time=md.timestepping.start_time+nyrs;

    %Output options
    if md.transient.ismovingfront == 1
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
							'Calvingratex','Calvingratey','CalvingCalvingrate'};
	

    else
	    md.transient.requested_outputs={'TotalSmb','SmbMassBalance',...
						    'IceVolume','IceVolumeAboveFloatation',...
						    'IceVolumeAboveFloatationScaled','GroundedAreaScaled',...
						    'FloatingAreaScaled','IceMass',...
						    'GroundedArea','FloatingArea','TotalFloatingBmb',...
						    'BasalforcingsFloatingiceMeltingRate',...
						    'GroundinglineMassFlux'}; %Gt/yr
    end
    
    % 					'GroundinglineHeight',...
    % 							'StrainRatexx','StrainRateyy','StrainRatexy',...
    % 							'StrainRateparallel','StrainRateperpendicular'};
    
	    md.verbose=verbose('solution',true,'module',false,'convergence',false);

        md.toolkits=toolkits();
		md.cluster=generic('name',oshostname,'np',4);

		%Solve
		md=solve(md,'Transient');

        savemodel(org,md);

end







