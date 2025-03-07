function export_csv = export_csv(md,flowline)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

data = [];
vel_avg = [];
vel_max = [];
distance = [];
icevolume = [];
time = [];
calving = [];
frontal_melt = [];
shelf_melt = [];
discharge = [];
smb = [];


if md.mesh.dimension == 2
    resolution = 100;
else
    resolution = [100 100];
end

%Calculate the front of the glacier along the flowline
[elementsL,xL,yL,zL,sL,hL]=SectionValues(md,md.results.TransientSolution(1).Surface,flowline,resolution);
for i = [1:length(sL)]
    if hL(i) > 1
        distance0 = sL(i); %Intital start value
        break
    end
end


for i = 1:length(md.results.TransientSolution)
    ice_area = find(md.results.TransientSolution(i).MaskIceLevelset<0);

    floating = find(md.results.TransientSolution(i).MaskIceLevelset<0 & md.results.TransientSolution(i).MaskOceanLevelset<0);
    time = [time md.results.TransientSolution(i).time];
    vel_avg = [vel_avg mean(md.results.TransientSolution(i).Vel(ice_area))];
    vel_max = [vel_max max(md.results.TransientSolution(i).Vel)];
    icevolume = [icevolume md.results.TransientSolution(i).IceVolume];
    %calving = [calving mean(md.results.TransientSolution(i).CalvingCalvingrate)];
    frontal_melt = [frontal_melt mean(md.results.TransientSolution(i).CalvingMeltingrate)];
    shelf_melt = [shelf_melt mean(md.results.TransientSolution(i).BasalforcingsFloatingiceMeltingRate(floating))];
    discharge = [discharge md.results.TransientSolution(i).GroundinglineMassFlux];
    smb = [smb md.results.TransientSolution(i).TotalSmb];


    [elementsL,xL,yL,zL,sL,hL]=SectionValues(md,md.results.TransientSolution(i).Surface,flowline,resolution);
    for i = [1:length(sL)]
        if hL(i) > 1
            distance = [distance (distance0-sL(i))];
            break
        end
    end
end

%% Sea level rise stuff

% Sea level rise contribution 
rho_ice = 917; % (kg/m^3) density of ice
rho_sw = 1030; % (kg/m^3) density of seawater
rho_fw = 1000; % (kg/m^3) density of freshwater
SL_potential=[];
SL_contrib=[];
grid_size=100;

if md.mesh.dimension == 3
    index=md.mesh.elements2d;   
    x=md.mesh.x2d;
    y=md.mesh.y2d;
else
    index=md.mesh.elements;   
    x=md.mesh.x;
    y=md.mesh.y;
end

%Not including impact of ocean freshening and/or firn density corrections
for i=1:(length(md.results.TransientSolution))
    base=md.results.TransientSolution(i).Base;
    surf=md.results.TransientSolution(i).Surface;
    % Calculate the freeboard height:
    fb = -base*rho_sw/rho_ice + base;
    fb(md.geometry.bed>=0) = md.geometry.bed(md.geometry.bed>=0);
    slr_thickness=surf-fb;
    slr_thickness(md.mask.ocean_levelset<0)=NaN;
    % Fix negative values around steep topography:
    slr_thickness(slr_thickness<0) = 0;
    % Project onto a fixed grid
    if md.mesh.dimension==3; slr_thickness=project2d(md, slr_thickness, 7); end
    slr_thickness_gridded=InterpFromMeshToGrid(index,x,y,slr_thickness,min(x):500:max(x),min(y):500:max(y),NaN);
    % Total ice volume available for SLR, in cubic meters:
    TotalSLRVolume=sum(slr_thickness_gridded(:)*(500*500), 'omitnan');
    % Total ice mass available for SLR (kg):
    TotalSLRmass = TotalSLRVolume*rho_ice;
    % Convert mass to Gt, and then to impact on SLR in milimeters:
    SLR_mm_per_Gt = (0.001/362)*1000;  % (mm/Gt) sea level potential per Gt water (or ice) from Meier et al. 2007
    SL_density_potential=md.results.TransientSolution(i).IceMass*SLR_mm_per_Gt*1e-12*((rho_sw-rho_fw)/rho_sw);
    SL_potential = [SL_potential ((TotalSLRmass*1e-12)*SLR_mm_per_Gt)+SL_density_potential];
    if i==1
        SL_contrib=[0];
    else
        sl_MassChange=(md.results.TransientSolution(i-1).IceMass-md.results.TransientSolution(i).IceMass)*SLR_mm_per_Gt*1e-12;
	    sl_density=sl_MassChange*((rho_sw-rho_fw)/rho_sw);
	    SL_contrib=[SL_contrib ((SL_contrib(i-1)+SL_potential(i-1)-SL_potential(i))+sl_density)];
	    %SL_contrib=[SL_contrib (SL_contrib(i-1)+SL_potential(i-1)-SL_potential(i))];
    end
end




headings = ["Time" "Vel_max (m/yr)" "Vel_avg (m/yr)" "Terminus_position_change (m)" "Ice Volume (Gt)" "SL Contribution (mm)" "SL Potential (mm)" "Grounded Melt (m/yr)" "Floating Melt (m/yr)" "Dischagre (Gt/yr)" "SMB (Gt/yr"];

data = [headings ; [time' vel_max' vel_avg' distance' icevolume' SL_contrib' SL_potential' frontal_melt' shelf_melt' discharge' smb']];

writematrix(data, ['./Outputs/' char(md.miscellaneous.name) '.csv'])

fprintf('File saved as ./Outputs/%s.csv\n',md.miscellaneous.name);

