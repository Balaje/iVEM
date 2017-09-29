function [l2,h1] = l2error(mesh,exactsol,uh,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute the L^2 error and the H^1 error  
%       [l2, h1] = l2error(mesh,exactsol,uh,k)
%    Input :
%           1. mesh structure: (see mesh/ folder).
%           2. exactsol: function handle (Example: 
%                       @(x,y)sin(pi*x)*sin(pi*y).
%           3. uh: Approximate solution vector.
%           4. k: Degree of the polynomial used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l2 = 0;
h1 = 0;
[area, centroid, diameter] = geo(mesh);
nel = size(mesh.elements);

l2consistency = 0;
h1consistency = 0;
stability = 0;

Q = quadrature_rule(6,2);
qw = Q(:,1); qx = Q(:,2); qy = Q(:,3);

modwrap = @(x,a) mod(x-1,a) + 1;
for i=1:nel
    usol = uh(mesh.elements{i});
    
    [~,~,~,G,D,B,~] = elem_matrices(mesh,i,k,@(x,y) 1);
    
    %%% Finding out the projection of the approx and exact solution
    pi0uh = (G\B)*usol;
    pi0u = l2projector(exactsol,i,mesh,k);
    
    if(k==1)
        nsides = size(mesh.elements{i},1);
        verts = mesh.vertices(mesh.elements{i},:);
        Dexact = zeros(nsides,1);
        for v = 1:nsides
            vert = verts(v,:);
            next = verts(modwrap(v+1, nsides), :);
            
            %%% To find the consistency term
            [l2err, h1err] = integrate2dtri_linear(centroid{i}, vert, next,...
                pi0uh, pi0u, diameter{i});
            l2consistency = l2consistency + l2err;
            h1consistency = h1consistency + h1err;
            
            Dexact(v) = exactsol(vert(1),vert(2));
        end
    elseif(k==2)
        eldofs = size(mesh.elements{i},1);
        nsides = (eldofs-1)/2;
        verts = mesh.vertices(mesh.elements{i}(1:end-1),:);
        Dexact = zeros(eldofs,1);
        for v=1:eldofs-1
            vert = verts(v,:);
            next = verts(modwrap(v+1, nsides), :);            
            %%% To find the consistency term
            if(v <= nsides)
                [l2err, h1err] = integrate2dtri_quadratic(centroid{i}, vert, next,...
                    pi0uh, pi0u, diameter{i});
                l2consistency = l2consistency + l2err;
                h1consistency = h1consistency + h1err;
            end
            
            Dexact(v) = exactsol(vert(1),vert(2));
        end
        %%% Find the last DOF.
        Dexact(eldofs) = 0;
        for v=1:nsides            
            vert = verts(v,:);
            next = verts(modwrap(v+1, nsides),:);   
            A = [1, centroid{i}; 1, vert; 1, next];
            triarea = 0.5*abs(det(A));
            for q=1:length(qw)
                xhat = (vert(1)-centroid{i}(1))*qx(q) + (next(1)-centroid{i}(1))*qy(q) + centroid{i}(1);
                yhat = (vert(2)-centroid{i}(2))*qx(q) + (next(2)-centroid{i}(2))*qy(q) + centroid{i}(2);
                Dexact(eldofs) = Dexact(eldofs) + 2*triarea*qw(q)*exactsol(xhat,yhat); 
            end
        end
        Dexact(eldofs) = (1/area{i})*Dexact(eldofs);
    end
    Dapprox = (D*(G\B))*usol;     
    %%% Finding out the stability term
    stability = stability + sum((Dexact-Dapprox).^2);
    
    %%% The equivalent norm corresponding to the L^2 norm.
    l2 = l2consistency + area{i}*stability;
    h1 = h1consistency + stability;
end

l2 = sqrt(l2);
h1 = sqrt(h1);
end

function [l2,h1] = integrate2dtri_linear(centroid, vert, next, uh, u, diameter)
Q = quadrature_rule(6,2);
qw = Q(:,1);  qx = Q(:,2);  qy = Q(:,3);
A = [1, centroid; 1, vert; 1, next];
area = 0.5*abs(det(A));
error = uh - u;

erroru = @(x,y) (error(1) + error(2)*(x-centroid(1))/diameter ...
    + error(3)*(y-centroid(2))/diameter);
errorux = error(2)/diameter;
erroruy = error(3)/diameter;

l2 = 0;
h1 = 0;
for q=1:6
    xhat = (vert(1)-centroid(1))*qx(q) + (next(1)-centroid(1))*qy(q) + vert(1);
    yhat = (vert(2)-centroid(2))*qx(q) + (next(2)-centroid(2))*qy(q) + vert(2);
    
    l2 = l2 + qw(q)*2*area*(erroru(xhat,yhat))^2;
    h1 = h1 + qw(q)*2*area*(errorux^2 + erroruy^2);
end
end

function [l2,h1] = integrate2dtri_quadratic(centroid, vert, next, uh, u, diameter)
Q = quadrature_rule(6,2);
qw = Q(:,1);  qx = Q(:,2);  qy = Q(:,3);
A = [1, centroid; 1, vert; 1, next];
area = 0.5*abs(det(A));
error = uh - u;

l2 = 0;
h1 = 0;
for q=1:6
    xhat = (vert(1)-centroid(1))*qx(q) + (next(1)-centroid(1))*qy(q) + vert(1);
    yhat = (vert(2)-centroid(2))*qx(q) + (next(2)-centroid(2))*qy(q) + vert(2);
    
    m = [1; (xhat-centroid(1))/diameter; (yhat-centroid(2))/diameter;
        ((xhat-centroid(1))/diameter)^2;
        ((xhat-centroid(1))/diameter)*((yhat-centroid(2))/diameter);
        ((yhat-centroid(2))/diameter)^2];
    Mx = [0; 1/diameter; 0; (2/diameter^2)*(xhat-centroid(1)); ...
        (1/diameter^2)*(yhat-centroid(2)); 0];
    My = [0; 0; 1/diameter; 0; (1/diameter^2)*(xhat-centroid(1)); ...
        (2/diameter^2)*(yhat-centroid(2))];
    
    l2 = l2 + qw(q)*2*area*(dot(error,m))^2;
    h1 = h1 + qw(q)*2*area*(dot(error,Mx)^2 + dot(error,My)^2);
end
end