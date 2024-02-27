function output = plot_outlines(md,flowline)
%PLOT_OUTLINES - used to plot the shape/velocities of a glacier along a
%defined flow line
%
%   usage = plot_outlines(md,'EXPFILE')

% if no flowline exp given, default to ryder.
% if exist(varargin,'flowline')
%     flowline = varargin{'flowline'};
% else
%     flowline = './Exp/ryder_2d.exp';
% end



if nargin < 2;
    flowline = './Exp/ryder_2d.exp';
else
    flowline = flowline;
end

if md.mesh.dimension()==2
    resolution = 100;
else
    resolution = [1000 1000];
end

%Resolution of transient run
res = (md.results.TransientSolution(end).time / length(md.results.TransientSolution)-1);
%Number of timesteps
ts = length(md.results.TransientSolution);

maxsurface = 0;

%starting slice
b = 1;
%Plot every x timestep
j = 12; 


figure();
%Make bottom plot for 2d view
subplot(5,5,[11:15 20:25])

fprintf('-- Interpolation glacier surface \n')

%ts = 60;

%Loop through every 12th timestep
for i = [b:j:ts]
        
    fprintf('  -- Timestep %i of %i... \n',i,ts);

     %load in surface and base of the glacier along flowline at 100m
     %resolution
     [elements,x,y,z,s,surf]=SectionValues(md,md.results.TransientSolution(i).Surface,flowline,resolution);
     [elements,x,y,z,s,base]=SectionValues(md,md.results.TransientSolution(i).Base,flowline,resolution);
     
     %Loop through the surface data and remove the very thin ice infront of
     %the glacier
    for m = [1:length(surf)]
        if surf(m)<1
            surf(m) = NaN;
        else
            surf(m)=0;
            limit = m;
            break
        end
    end

    %flip the x and y of the surface values
    y_line = flip(surf(limit:end));
    x_line = flip(s(limit:end));

    %Loop through the base data to remove thin ice
    for m = [1:length(base)]
        if base(m)>(-2);
            base(m) = NaN;
        else
            base(m)=0;
            limit = m;
            break
        end
    end

    %Concatenate the surface and base data so they appear joined up in the
    %final plot
    y_line = [y_line; base(limit:end)];
    x_line = [x_line; s(limit:end)];
    
    if i == 1 || i == ts
        plot(x_line,y_line,'LineWidth',3,'DisplayName',string(i))
    else
        plot(x_line,y_line,'DisplayName',string(i))
    end

    maxsurface = max(maxsurface,max(y_line));

    hold on
  
end


xlim([0 s(end)]);
%ylim([-1000 2000]);
xticklabels(xticks/1000);
%Create cmap based on the number of timesteps
cmap = parula((length([b:j:ts]))+1);
set(gca,'colororder',cmap); ax = gca; ax.FontSize = 14; 
ylabel('Elevation (m)', 'FontSize', 14)
xlabel('Distance along flowline (km)', 'FontSize', 14)
c = colorbar; c.Location = "southoutside";
c.Ticks = [0, 0.5, 1];




if b==1
    if md.results.TransientSolution(1).time > 2024
        c.TickLabels = {'2024', string(((md.results.TransientSolution(ts).time)/2)+2024), string(md.results.TransientSolution(ts).time+2024)};
    else
        c.TickLabels = {'0', string(((md.results.TransientSolution(ts).time)/2)), string(md.results.TransientSolution(ts).time)};
    end
else
    c.TickLabels = {string(md.results.TransientSolution(b).time), string((md.results.TransientSolution(b).time+md.results.TransientSolution(ts).time)/2), string(md.results.TransientSolution(ts).time)};
end

c.Label.String = 'Simulation Years'; c.Label.FontSize = 16; c.FontSize = 14;


%Get bedrockdata
[elemets,x,y,z,s,bed] = SectionValues(md,md.geometry.bed,flowline,resolution);
plot(s,bed,'black', LineWidth=1)
fill_x = [s; flip(s)];
fill_y = [bed; (min(bed))*ones(length(bed),1)];
fill(fill_x,fill_y, [100 58 38]./255, 'FaceAlpha', 0.5)

set ( gca, 'xdir', 'reverse' )
%xt = get(gca, 'XTick');
%set(gca, 'XTickLabel', str([0:2:20]));

maxsurface =ceil(maxsurface/500)*500;
ylim([min(bed) maxsurface])


%Make top plot for velocity profile
subplot(5,5,[1:5 6:10])

fprintf('-- Interpolation velocities \n')

%Loop through every timestep
for i = [b:j:ts]

    fprintf('   -- Timestep %i of %i... \n',i,ts);

    %load in velocity data
    [elemets,x,y,z,s,vel] = SectionValues(md,md.results.TransientSolution(i).Vel,flowline,resolution);
    
    maxvel = 0;
    %Assumining the max velocity is the front of the glacier, remove any
    %slow moving ice
    for i = [1:length(vel)]
        if vel(i) == 0 ||  (vel(i)<(vel(i+1)*0.97))    %vel(i) < max(vel)
            vel(i) = NaN;
        else
            maxvel = max(maxvel,max(vel));
            break
        end
    end
    plot(s,vel)
    hold on

end

maxvel = ceil(maxvel/1000)*1000;

ylim([0 maxvel])
xlim([0 s(end)]);
xticklabels(xticks/1000);
ax = gca; ax.FontSize = 14; 
ylabel('Velocity (m/yr)','FontSize',14)
set(gca,'colororder',cmap);


set ( gca, 'xdir', 'reverse' )
%xt = get(gca, 'XTick');
%set(gca, 'XTickLabel', fliplr(xt/1000))




