
% Plotting code for ISSM simulations

% Plotting line graphs of average velocity/ SMB etc
vel = [];
smb = [];
calving = [];
frontal_melt = [];
shelf_melt = [];
volume = [];
volume_float = [];
gl_flux=[];

output_steps = length(md.results.TransientSolution)-1; % We want to start from 0
front_area=find(md.mesh.y>(-1*10e5) & md.mask.ice_levelset<0);

for i=1:(output_steps+1)
    vel = [vel mean(md.results.TransientSolution(i).Vel(front_area))];
    %smb = [smb md.results.TransientSolution(i).TotalSmb];
    calving = [calving mean(md.results.TransientSolution(i).CalvingCalvingrate)];
    frontal_melt = [frontal_melt mean(md.results.TransientSolution(i).CalvingMeltingrate)];
    shelf_melt = [shelf_melt mean(md.results.TransientSolution(i).BasalforcingsFloatingiceMeltingRate)];
    volume = [volume md.results.TransientSolution(i).IceVolume];
    %volume_float = [volume_float md.results.TransientSolution(i).FloatingArea];
    %gl_flux=[gl_flux md.results.TransientSolution(i).GroundinglineMassFlux];
end

figure()
subplot(5,1,1);
plot([0:output_steps], vel, 'color', 'b', 'linewidth', 2);
hold on;
plot([0:output_steps], calving, 'color', 'g', 'linewidth', 2);
title('Calving rate (green) and mean frontal velocity (blue) (m/yr)');
xlabel('Simulation years');

subplot(5,1,2);
plot([0:output_steps], volume, 'color', 'r', 'linewidth', 2)
title('Total ice volume (m3)');
xlabel('Simulation years');

subplot(5,1,3);
plot([0:output_steps], frontal_melt, 'color', 'k', 'linewidth', 2);
hold on;
plot([0:output_steps], shelf_melt, 'color', 'm', 'linewidth', 2); 
title('Mean melt at grounded fronts (black) and under floating ice (pink) (m/yr)');
xlabel('Simulation years');

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
for i=1:(output_steps+1)
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

subplot(5,1,4);
plot([2:output_steps], SL_potential(2:output_steps), 'color', [0.3010 0.7450 0.9330], 'linewidth', 2); 
title('Sea level potential locked up in glacier year-on-year (mm)');
xlabel('Simulation years');
subplot(5,1,5);
plot([2:output_steps], SL_contrib(2:output_steps), 'color', [0.9290 0.6940 0.1250], 'linewidth', 2); 
title('Cumulative Sea level contribution (mm)');
xlabel('Simulation years');




