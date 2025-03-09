function output = plot_outlines(md,flowline,observations)
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
    resolution = 500;
else
    resolution = [1000 1000];
end

%Resolution of transient run
res = (md.results.TransientSolution(end).time / length(md.results.TransientSolution)-1);
%Number of timesteps
ts = length(md.results.TransientSolution);

ice_array = [];

maxsurface = 0;

%starting slice
b = 1;
%Plot monthly timestep
j = 12;  

figure();
%Make bottom plot for 2d view
subplot(5,5,[11:15 20:25])

fprintf('-- Interpolation glacier surface \n')

%ts = 168;

%Loop through every 12th timestep
for i = [b:j:ts]
        
    fprintf('  -- Timestep %i of %i... \n',i,ts);

     %load in surface and base of the glacier along flowline at 100m
     %resolution
     [elements,x,y,z,s,surf]=SectionValues(md,md.results.TransientSolution(i).Surface,flowline,resolution);
     [elements,x,y,z,s,base]=SectionValues(md,md.results.TransientSolution(i).Base,flowline,resolution);
     [elements,x,y,z,s,ice]=SectionValues(md,md.results.TransientSolution(i).MaskIceLevelset,flowline,resolution);

     ice_array = [ice_array ice]; %save positions of ice
     pos = find(ice<0); %find where there is ice
     

    %flip the x and y of the surface values
    y_line = flip(surf(pos));
    x_line = flip(s(pos));

    %Concatenate the surface and base data so they appear joined up in the
    %final plot
    y_line = [y_line; base(pos)];
    x_line = [x_line; s(pos)];
    
    if i == 1 || i == ts
        plot(x_line,y_line,'LineWidth',2,'DisplayName',string(i))
    else
        plot(x_line,y_line,'DisplayName',string(i))
    end

    disp(max(y_line));
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



c.TickLabels = {string(md.results.TransientSolution(b).time), string((md.results.TransientSolution(b).time+md.results.TransientSolution(ts).time)/2), string(md.results.TransientSolution(ts).time)};

c.Label.String = 'Simulation Years'; c.Label.FontSize = 16; c.FontSize = 14;



%Get bedrockdata
[elemets,x,y,z,s,bed] = SectionValues(md,md.geometry.bed,flowline,resolution);
minbed = floor(min(bed)/500)*500;
plot(s,bed,'black', LineWidth=1)
fill_x = [s; flip(s)];
fill_y = [bed; (minbed)*ones(length(bed),1)];
fill(fill_x,fill_y, [100 58 38]./255, 'FaceAlpha', 0.5)

set ( gca, 'xdir', 'reverse' )
%xt = get(gca, 'XTick');
%set(gca, 'XTickLabel', str([0:2:20]));

maxsurface =ceil(maxsurface/500)*500;
ylim([minbed maxsurface])


%Make top plot for velocity profile
subplot(5,5,[1:5 6:10])
maxvel = 0;

fprintf('-- Interpolation velocities \n')
count = 1;
%Loop through every timestep
for i = [b:j:ts]
    fprintf('   -- Timestep %i of %i... \n',i,ts);

    %load in velocity data
    [elemets,x,y,z,s,vel] = SectionValues(md,md.results.TransientSolution(i).Vel,flowline,resolution);
      
    pos = find(ice_array(:,count)<0);
    
    plot(s(pos),vel(pos))
    hold on

    maxvel = max(maxvel,max(vel));
    
    count = count+1;
end




%Adapt the scale depending on if the glacier is fast (>3000) or slow
%(<3000)
if maxvel > 3000
    maxvel = ceil(maxvel/1000)*1000;
else
    maxvel = ceil(maxvel/100)*100;
end


ylim([0 maxvel])
xlim([0 s(end)]);
xticklabels(xticks/1000);
ax = gca; ax.FontSize = 14; 
ylabel('Velocity (m/yr)','FontSize',14)
set(gca,'colororder',cmap);

set ( gca, 'xdir', 'reverse' )
%xt = get(gca, 'XTick');
%set(gca, 'XTickLabel', fliplr(xt/1000))




