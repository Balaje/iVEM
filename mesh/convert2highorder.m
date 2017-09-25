function V = convert2highorder(mesh, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that converts the input mesh to an higher order problem by
% adding new degrees of freedom
%       Input :
%               mesh : Input mesh structure containing
%                       mesh.elements
%                       mesh.vertices
%                       mesh.boundary
%       Output :
%               V : New mesh structure containing the added data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(k == 1)
    V = mesh;
elseif(k == 2)
    ndofs = length(mesh.vertices);
    nel = length(mesh.elements);
    
    total_sides = sum(cellfun('length',mesh.elements));
    sides = zeros(total_sides,2);
    sideid = 0;
    for i=1:nel
        verts = mesh.elements{i};
        elementsides = [verts, circshift(verts,-1)];
        nsides = length(mesh.elements{i});                
        sides(sideid+1:sideid+nsides,:) = ...
            sides(sideid+1:sideid+nsides,:) + elementsides;        
        sideid = sideid + (nsides);
    end  
    numbering = zeros(total_sides,1);    
    count = 1;
    for i=1:length(sides)
        if(numbering(i) == 0)
            numbering(i) = count;
            count = count+1;
        end        
        sidei = sides(i,:);
        for j=1:length(sides)
            sidej = sides(j,:);            
            if(isequal(sidei, circshift(sidej',1)'))
                if(numbering(j) == 0)
                    numbering(j) = numbering(i);                    
                end
            end
        end        
    end    
    midpoint_X = 0.5*(mesh.vertices(sides(:,1),1)+mesh.vertices(sides(:,2),1));
    midpoint_Y = 0.5*(mesh.vertices(sides(:,1),2)+mesh.vertices(sides(:,2),2));
    mp = [midpoint_X, midpoint_Y];
    
    %%%%%%% Creation of new mesh %%%%%%%%% 
    newnumbering = numbering + ndofs;
    
    %%%%%%% Create new-element cell array by adding the new points
    elementdofids = (1:ndofs) + max(newnumbering);
    newelements = mesh.elements;
    newnodes = newnumbering(1:length(mesh.elements{1}));
    newelements{1} = [newelements{1}; newnodes; elementdofids(1)];
    for i=2:nel
        newnodes = newnumbering(length(vertcat(mesh.elements{1:i-1}))+1 ...
        :length(vertcat(mesh.elements{1:i-1}))+length(mesh.elements{i}));
    
       newelements{i} = [newelements{i}; newnodes; elementdofids(i)]; 
    end
    
    %%%%%%% Add the new vertices 
    noofnewvertices = max(numbering);
    newvertices = zeros(size(mesh.vertices,1)+noofnewvertices,2);
    newvertices(1:size(mesh.vertices,1),:) = mesh.vertices;
    for id = 1:noofnewvertices
        index = find(numbering == id, 1);
        newvertices(size(mesh.vertices,1)+id, :) = mp(index, :);
    end
    
    %%%%%%% Add the boundary vertices  
    newboundary = zeros(2*length(mesh.boundary),1);
    newboundary(1:length(mesh.boundary)) = mesh.boundary;
    count = length(mesh.boundary)+1;
    for i=1:total_sides
        if(ismember(sides(i,1),mesh.boundary)...
                && ismember(sides(i,2),mesh.boundary))
            newboundary(count) = newnumbering(i);
            count = count + 1;
        end
    end
    
   V = struct('boundary',{newboundary},...
       'elements',{newelements},'vertices',{newvertices});   
end

end