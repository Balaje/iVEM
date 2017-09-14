function [area, centroid, diameter] = geo(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function that generates the geometric data for the elements
%
%       geo(mesh, el_id):
%           Input:
%               1. (mesh)mesh-data : Structure argument containing the data 
%                               of the mesh. For more information, see doc
%               2. el_id     : The id of the desired element.
%           Output: 
%               The area, centroid(1x2), diameter of the element
%
%       geo(mesh)
%           Input: 
%               1. (mesh)mesh-data : Structure argument as above.
%               2. Cell array containing the data for the area, centroid,
%                  diameter. For more information, type 
%                           help cell
%                  in the terminal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin==2) %% IF THE USER WANTS DATA FOR THE PARTICULAR ELEMENT.
    mesh = varargin{1};
    elid = varargin{2};    
    % Collect the required information
    elem = mesh.elements{elid};
    verts = mesh.vertices(elem,:);
    nsides = length(elem);    
    % 1. Calculate the area of the element.
    area_comp = verts(:,1).*verts([2:end,1],2)...
        - verts([2:end,1],1).*verts(:,2);
            % This will be required to compute the centroid and area.
    area = 0.5*abs(sum(area_comp));    
    % 2. Calculate the centroid of the element
    centroid = sum((verts+verts([2:end,1],:)).*repmat(area_comp,1,2))...
        /(6*area);    
    % 3. Compute the diameter of the element.
    diameter = 0; %
    for i = 1:(nsides-1)
        for j = (i+1):nsides
            diameter = max(diameter, norm(verts(i,:) - verts(j,:)));
        end
    end
    
    
    
else %% IF THE USER WANTS DATA FOR ALL THE ELEMENTS    
    % Initialize the data-structures
    mesh = varargin{1};
    nel = length(mesh.elements);    
    area = cell(length(nel));
    centroid = cell(length(nel));
    diameter = cell(length(nel));    
    % Start computing for all elements
    for el = 1:nel
        % Collect the required information
        elem = mesh.elements{el};
        verts = mesh.vertices(elem,:);
        nsides = length(elem);
        % 1. Calculate the area of the element.
        area_comp = verts(:,1).*verts([2:end,1],2)...
            - verts([2:end,1],1).*verts(:,2);
        % This will be required to compute the centroid and area.
        area{el} = 0.5*abs(sum(area_comp));       
        % 2. Calculate the centroid of the element
        centroid{el} = sum((verts+verts([2:end,1],:)).*repmat(area_comp,1,2))...
            /(6*area{el});        
        % 3. Compute the diameter of the element.
        diameter{el} = 0; %
        for i = 1:(nsides-1)
            for j = (i+1):nsides
                diameter{el} = max(diameter{el}, norm(verts(i,:) - verts(j,:)));
            end
        end
    end    
end
end