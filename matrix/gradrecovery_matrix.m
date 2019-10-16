function [Mx, My, Q] = gradrecovery_matrix(mesh,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that returns the Gradient Recovery Matrices  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nel = size(mesh.elements,1);
nodes = size(mesh.vertices,1);

[~, ~, diameter] = geo(mesh);

if(k==1)
    Mx = zeros(nodes);
    My = zeros(nodes);
    Q = zeros(nodes);
elseif(k==2)
    ndofs = nel+nodes;
    Mx = zeros(ndofs);
    My = zeros(ndofs);
    Q = zeros(ndofs);
end

for el = 1:nel
    nsides = length(mesh.elements{el});
    %usol = uh(mesh.elements{el});
    
    [~,~,Me,G,~,B,~] = elem_matrices(mesh,el,k,@(x,y) 1);
    Qe = diag(sum(Me,1));
    
    if(k==1)
        Dx = repmat([0, 1/diameter{el}, 0], nsides, 1);
        Dy = repmat([0, 0, 1/diameter{el}], nsides, 1);
    elseif(k==2)
        elem = mesh.elements{el};
        eldofs = length(elem);
%         verts = mesh.vertices(elem(1:(eldofs-1)/2),:); % Vertices
%         nsides = length(verts);
        v = mesh.vertices(mesh.elements{el}(1:end-1),:); % Vertices + Midpoints
        npolys = 6; % Number of scaled monomials        
        [~, centroid, diameter] = geo(mesh,el);
        
        Dx = zeros(eldofs,npolys);
        Dy = zeros(eldofs,npolys);
        for i=1:eldofs-1
            vert = v(i,:);
            Dx(i,:) = [0, 1/diameter, 0, 2/diameter^2*(vert(1)-centroid(1))...
                , 1/diameter^2*(vert(2)-centroid(2)), 0];
            Dy(i,:) = [0, 0, 1/diameter, 0, 1/diameter^2*(vert(1)-centroid(1)) ...
                , 2/diameter^2*(vert(2)-centroid(2))];
        end
        % Last DOF (Moments)
        Dx(eldofs,:) = [0, 1/diameter, 0, 0, 0, 0];
        Dy(eldofs,:) = [0, 0, 1/diameter, 0, 0 ,0];
    end
    pi0uh = (G\B);
    
    Mxe = Qe*Dx*pi0uh;
    Mye = Qe*Dy*pi0uh;
    
    Mx(mesh.elements{el},mesh.elements{el}) = ...
        Mx(mesh.elements{el},mesh.elements{el}) + Mxe;
    My(mesh.elements{el},mesh.elements{el}) = ...
        My(mesh.elements{el},mesh.elements{el}) + Mye;
    
    Q(mesh.elements{el},mesh.elements{el}) = ...
        Q(mesh.elements{el},mesh.elements{el}) + Me;
end


end