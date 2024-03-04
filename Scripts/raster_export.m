function raster_export(md, variable, timesteps, savefile)
%%%%%% %Interpolate ISSM output from model mesh to netCDF grid

%% Example usage for creating file 'test_export.nc' including surface results for 5 timesteps:
%% raster_export(md, 1, 'Surface', [1:5], 7, 'ssp585_surface_2100')

% Layer input variable is unused if the model is 2D

isverbose=1;
data_loc='md.results.TransientSolution';

index=md.mesh.elements;
X=md.mesh.x;
Y=md.mesh.y;
message = strcat('   -- raster Export: Model is 2D, layer argument being ignored');
if isverbose; disp(message); end



% Create grid
xlim = [min(X) max(X)];
ylim = [min(Y) max(Y)];
gridcell = 500;
x_g = xlim(1):gridcell:xlim(2);
y_g = ylim(1):gridcell:ylim(2);

message = strcat('   -- raster Export:',{' '}, 'Reading timestep',{' '},string(timesteps(1)),{' '},'to timestep',{' '},string(timesteps(end)),{' '},'of variable', {' '},variable);
if isverbose, disp(message); end
data = [];
for i=1:length(timesteps)
    name = char(strcat(string(data_loc), '(',string(timesteps(i)), ')', '.', variable));
    var = eval(name);
    tmp = InterpFromMeshToGrid(index,X,Y,var,x_g(:),y_g(:),NaN);
    tmp=tmp';
    data(:, :, i) = [tmp];
end


if isverbose, disp('   -- raster Export: Writing to netCDF file'); end
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
ncid=netcdf.create(strcat('Outputs/', savefile, '.nc'),mode);

netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'Projection', "Polar Stereographic North (70N, 45W)");
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'Conventions', "CF-1.7");
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'proj4', "+init=epsg:3413");
x_dimID = netcdf.defDim(ncid,'x',length(x_g));
y_dimID = netcdf.defDim(ncid,'y',length(y_g));
time_dimID = netcdf.defDim(ncid,'time',length(timesteps(1):timesteps(end)));
data_var_id = netcdf.defVar(ncid,variable,'NC_FLOAT',[x_dimID,y_dimID,time_dimID]);
x_var_id = netcdf.defVar(ncid,'x','DOUBLE',x_dimID);
y_var_id = netcdf.defVar(ncid,'y','DOUBLE',y_dimID);
time_var_id = netcdf.defVar(ncid,'time','DOUBLE',time_dimID);
% Set standard_name and units attributes for x and y variables
netcdf.putAtt(ncid, x_var_id, 'standard_name', 'projection_x_coordinate');
netcdf.putAtt(ncid, x_var_id, 'units', 'm');
netcdf.putAtt(ncid, y_var_id, 'standard_name', 'projection_y_coordinate');
netcdf.putAtt(ncid, y_var_id, 'units', 'm');
netcdf.putAtt(ncid, time_var_id, 'standard_name', 'time');
netcdf.putAtt(ncid, time_var_id, 'units', 'years');
netcdf.endDef(ncid);

netcdf.putVar(ncid,x_var_id,x_g);
netcdf.putVar(ncid,y_var_id,y_g);
netcdf.putVar(ncid,time_var_id,timesteps(1):timesteps(end));
netcdf.putVar(ncid,data_var_id,data);

netcdf.close(ncid);

% disp('   -- raster Export: Creating GeoTiff using gdal')
% system(sprintf('gdal_translate -a_srs EPSG:3413 NETCDF:%s.nc:%s %s.tif', savefile, variable, savefile));
% disp('   -- raster Export: All done')



















