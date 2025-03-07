function output = plot_output(md,flowline)
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
distance = [];
time=[];
time2=[];
discharge = [];
smb=[];
icearea = [];


output_steps = length(md.results.TransientSolution)-1; % We want to start from 0
%front_area=find(md.mesh.y>(-1*10e5) & md.mask.ice_levelset<0);

%Set reslution of interpolation for flowline
if md.mesh.dimension == 2
    resolution = 1000;
else
    resolution = [1000 1000];
end

%Calculate the front of the glacier along the flowline
[elementsL,xL,yL,zL,sL,hL]=SectionValues(md,md.results.TransientSolution(1).Surface,flowline,resolution);
for i = [1:length(sL)]
    if hL(i) > 1
        distance0 = sL(i); %Intital start value
        break
    end
end



for i=1:(output_steps+1)
    time = [time md.results.TransientSolution(i).time];
    vel = [vel max(md.results.TransientSolution(i).Vel)];
    %smb = [smb md.results.TransientSolution(i).TotalSmb];
    floating = find(md.results.TransientSolution(i).MaskIceLevelset<0 & md.results.TransientSolution(i).MaskOceanLevelset<0);
    %calving = [calving mean(md.results.TransientSolution(i).CalvingCalvingrate)];
    frontal_melt = [frontal_melt mean(md.results.TransientSolution(i).CalvingMeltingrate)];
    shelf_melt = [shelf_melt mean(md.results.TransientSolution(i).BasalforcingsFloatingiceMeltingRate(floating))];
    volume = [volume md.results.TransientSolution(i).IceVolume];
    %volume_float = [volume_float md.results.TransientSolution(i).FloatingArea];
    %gl_flux=[gl_flux md.results.TransientSolution(i).GroundinglineMassFlux];

    [elementsL,xL,yL,zL,sL,hL]=SectionValues(md,md.results.TransientSolution(i).Surface,flowline,resolution);
    for i = [1:length(sL)]
        if hL(i) > 1
            distance = [distance (distance0-sL(i))];
            break
        end
    end
end

for i = 1:12:length(md.results.TransientSolution)-3
    if i ~= length(md.results.TransientSolution)
        year_smb = (md.results.TransientSolution(i:i+11).TotalSmb) ;
        year_smb = sum(year_smb);
        year_discharge = md.results.TransientSolution(i:i+11).GroundinglineMassFlux;
        year_discharge = sum(year_discharge);
        smb = [smb year_smb];
        discharge = [discharge year_discharge];
        time2 = [time2 md.results.TransientSolution(i).time];
    end
end

% if time(1)>2024
%     time(1) = 2024;
%     time2(1) = 2024;
% end


figure()
ax1 = subplot(5,1,1);
plot(time, vel, 'color', 'b', 'linewidth', 2);
xlim([time(1) time(end)]);
hold on;
yyaxis right
%plot(time, calving, 'color', 'g', 'linewidth', 2);
title('Max velocity (blue) (m/yr)'); %and Calving rate (green)
xlabel('Simulation years');

ax2 = subplot(5,1,2);
plot(time, distance, 'color', "#0072BD", 'linewidth', 2)
xlim([time(1) time(end)]);
hold on 
yyaxis right;
plot(time, volume, 'color', 'r', 'linewidth', 2)
title('Change in terminus position (blue) (m) and total ice volume (red) (m3)');
xlabel('Simulation years');

ax3 = subplot(5,1,3);
plot(time, frontal_melt, 'color', 'k', 'linewidth', 2);
xlim([time(1) time(end)]);
hold on;
plot(time, shelf_melt, 'color', 'm', 'linewidth', 2); 
title('Mean melt at grounded fronts (black) and under floating ice (pink) (m/yr)');
xlabel('Simulation years');


ax4 = subplot(5,1,4);
plot(time2+1,smb, 'color', "#D95319", 'linewidth' , 2);
hold on
plot(time2+1,discharge, 'color', "#77AC30", 'linewidth', 2);
title('Annual SMB (orange) (Gt/yr) and Annual Discharge (green) (Gt/yr)');
xlim([time(1) time(end)]);



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

ax5 = subplot(5,1,5);
plot(time, SL_potential, 'color', [0.3010 0.7450 0.9330], 'linewidth', 2); 
xlim([time(1) time(end)]);
title('Sea level potential locked up in glacier year-on-year (mm) and Cumulative Sea level contribution (mm)');
xlabel('Simulation years');
yyaxis right;
plot(time, SL_contrib, 'color', [0.9290 0.6940 0.1250], 'linewidth', 2); 
xlim([time(1) time(end)]);


linkaxes([ax1,ax2, ax3, ax4, ax5], 'x')

