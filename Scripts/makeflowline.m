function mesh = makemesh(data,filename)


mesh_data = readtable(data,"Delimiter",',');
mesh = struct('x',mesh_data.X,'y',mesh_data.Y);

expwrite(mesh,filename)

expdisp(filename)
%Make mesh