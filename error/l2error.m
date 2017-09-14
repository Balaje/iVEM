function [l2,h1] = l2error(mesh,exactsol,uh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute the L^2 error and the H^1 error  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l2 = 0;
h1 = 0;
[area, centroid, diameter] = geo(mesh);
nel = size(mesh.elements);    

l2consistency = 0;
h1consistency = 0;
stability = 0;

modwrap = @(x,a) mod(x-1,a) + 1;
for i=1:nel
    usol = uh(mesh.elements{i});
    
    [~,~,~,G,D,B,~] = elem_matrices(mesh,i,1,@(x,y) 1);
    
    %%% Finding out the projection of the approx and exact solution
    pi0uh = (G\B)*usol;
    pi0u = l2projector(exactsol,i,mesh);        
        
    nsides = size(mesh.elements{i},1);
    verts = mesh.vertices(mesh.elements{i},:);    
    Dexact = zeros(nsides,1);    
    for v = 1:nsides
        vert = verts(v,:);
        next = verts(modwrap(v+1, nsides), :);        

        %%% To find the consistency term
        [l2err, h1err] = quadrature(centroid{i}, vert, next,...
                        pi0uh, pi0u, diameter{i}); 
        l2consistency = l2consistency + l2err;
        h1consistency = h1consistency + h1err;
        
        Dexact(v) = exactsol(vert(1),vert(2));        
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

function [l2,h1] = quadrature(centroid, vert, next, uh, u, diameter)

qw = [1/6, 1/6, 1/6];
qx = [2/3, 1/6, 1/6];
qy = [1/6, 1/6, 2/3];

A = [1, centroid; 1, vert; 1, next];
area = 0.5*abs(det(A));
error = uh - u;

erroru = @(x,y) (error(1) + error(2)*(x-centroid(1))/diameter ...
    + error(3)*(y-centroid(2))/diameter);
errorux = error(2)/diameter;
erroruy = error(3)/diameter;

l2 = 0;
h1 = 0;
for q=1:3
    xhat = (vert(1)-centroid(1))*qx(q) + (next(1)-centroid(1))*qy(q) + vert(1);
    yhat = (vert(2)-centroid(2))*qx(q) + (next(2)-centroid(2))*qy(q) + vert(2);
    
    l2 = l2 + qw(q)*2*area*(erroru(xhat,yhat))^2;
    h1 = h1 + qw(q)*2*area*(erroru(xhat,yhat)^2 + errorux^2 + erroruy^2);
end
end