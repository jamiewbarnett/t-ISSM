%% Script for delineating the ice-land boundary 

function output = find_iceLandBoundary(md, is2D)

% Read in the mask 
M = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask', 'nearest',...
    '../Model_Data/BedMachineGreenland-v5.nc');
% 0 is water, 1 is land, 2 is grounded ice, 3 is floating ice

if is2D==1
    nodes_per_element=width(md.mesh.elements);
    landnodes=(M==1);
    landelement=sum(landnodes(md.mesh.elements), 2);
    icennodes=md.mask.ice_levelset < 0;
    iceelement=sum(icennodes(md.mesh.elements), 2);
 else
    nodes_per_element=width(md.mesh.elements2d);
    landnodes=(M==1);
    landelement=sum(landnodes(md.mesh.elements2d), 2);
    icesurf=project2d(md, md.mask.ice_levelset, md.mesh.numberoflayers); % Check if ice at surf
    icenodes=icesurf < 0;
    iceelement=sum(icenodes(md.mesh.elements2d), 2);
 end
               
  land_ice_boundary = zeros(length(iceelement), 1);
  for element=1:length(iceelement) 
     if iceelement(element) == nodes_per_element || landelement(element) == 0 || landelement(element) == nodes_per_element
         land_ice_boundary(element) = 0;
     else
         land_ice_boundary(element) = 1;
     end
  end

  if is2D==0
    land_ice_nodes=[];
    for element=1:length(land_ice_boundary)
        if land_ice_boundary(element)==0
            continue
        else 
        land_ice_nodes=[land_ice_nodes md.mesh.elements2d(element, :)];
        end
     end
  else
    land_ice_nodes=[];
    for element=1:length(land_ice_boundary)
        if land_ice_boundary(element)==0
            continue
        else 
        land_ice_nodes=[land_ice_nodes md.mesh.elements(element, :)];
        end
    end
  end


land_ice_nodes=unique(land_ice_nodes(:));
output=land_ice_nodes;

