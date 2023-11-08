function mesh = makemesh(data,filename)


mesh_data = readtable(data,"Delimiter",';');
mesh = struct('x',mesh_data.X,'y',mesh_data.Y);

if (mesh.x(1) ~= mesh.x(end)) && (mesh.y(1) ~= mesh.y(end))
    mesh.x = [mesh.x; mesh.x(1)];
    mesh.y = [mesh.y; mesh.y(1)];
end

expwrite(mesh,filename)

expdisp(filename)
%Make mesh