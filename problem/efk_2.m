function [u,v] = efk_2(mesh, f, g_Du, g_Dv, u0, tf, t0, delt, plotmesh)

nverts = length(mesh.vertices);
nel = length(mesh.elements);
ndof = nel+nverts;

ntimes = fix((tf-t0)/delt);
modwrap = @(x,a) mod(x-1,a) + 1;

un = zeros(ndof,1);
un(1:nverts) = u0(mesh.vertices(:,1), mesh.vertices(:,2));
vn = zeros(length(un),1);
v = zeros(length(vn),1);

for el=1:nel
   eldofs = length(mesh.elements{el});
   nsides = (eldofs-1)/2;
   verts = mesh.vertices(mesh.elements{el}(1:nsides),:);
   [area, centroid, ~] = geo(mesh,el);
   for v=1:nsides
       vert = verts(v,:);
       next = verts(modwrap(v+1, nsides), :);
       C = {centroid, vert, next};
       un(ndof-nverts+el) = un(ndof-nverts+el) + (1/area)*(integrate2dfunctri(u0,C,6));
   end
end

uu = ones(length(un),1);

t = 0;

for time = 1:ntimes
    %%% Begin Tolerance loop
    error = 100;
    tol = 1e-7;
    
    t = t + delt;
    while error > tol
        K = sparse(ndof,ndof);
        M = sparse(ndof,ndof);
        C = sparse(ndof,ndof);
        F = zeros(ndof,1);
        
        for el = 1:nel
            g = @(x,y) f(x,y,t);
            
            [~, centroid, diameter] = geo(mesh, el);
            [Ke, Fe, Me, G, ~, B, ~] = elem_matrices(mesh,el,2,g);
            
            K(mesh.elements{el},mesh.elements{el}) = ...
                K(mesh.elements{el},mesh.elements{el}) + Ke;
            
            M(mesh.elements{el},mesh.elements{el}) = ...
                M(mesh.elements{el},mesh.elements{el}) + Me;
            
            F(mesh.elements{el},1) = F(mesh.elements{el},1) + Fe;
            %%%%% The nonlinear term
            eldofs = length(mesh.elements{el});
            nsides = (eldofs-1)/2;
            verts = mesh.vertices(mesh.elements{el}(1:nsides),:);
            projection = (G\B)*uu(mesh.elements{el});
            
            projector = G\B;
            Ce = efk_matrix(centroid, diameter, projector, projection, nsides, verts);
            
            C(mesh.elements{el},mesh.elements{el}) = ...
                C(mesh.elements{el},mesh.elements{el}) + Ce;
        end
        
        boundaryvals_p = g_Du(mesh.vertices(mesh.boundary,1),...
            mesh.vertices(mesh.boundary,2), t);
        
        F = delt*F + M*un...
            - (M(:,mesh.boundary) + 0.01*delt*K(:,mesh.boundary)*(M(:,mesh.boundary)\K(:,mesh.boundary))...
                + delt*K(:,mesh.boundary) + delt*C(:,mesh.boundary) - delt*M(:,mesh.boundary))*boundaryvals_p;
        
        freenodes = mesh.boundary;
        solnodes = setdiff(1:ndof, freenodes);
        
        u = zeros(ndof,1);
        u(solnodes) = (M(solnodes,solnodes) + 0.01*delt*K(solnodes,solnodes)*...
            (M(solnodes,solnodes)\K(solnodes,solnodes)) + delt*K(solnodes,solnodes)...
            + delt*C(solnodes,solnodes) - delt*M(solnodes,solnodes))\F(solnodes);
        
        u(freenodes) = g_Du(mesh.vertices(mesh.boundary,1),...
            mesh.vertices(mesh.boundary,2), t);
        
        error = max(abs(uu - u));
        uu = u;
    end
    
    un = u;
    v(mesh.boundary) = g_Dv(mesh.vertices(mesh.boundary,1),...
            mesh.vertices(mesh.boundary,2), t);
    v(solnodes) = (M(solnodes,solnodes)\K(solnodes,solnodes))*u(solnodes);
    
%     cla
%     plot_solution(plotmesh,u(1:length(plotmesh.vertices)));
%     grid on
    str = ['Solution at time t = ',num2str(t)];
    disp(str)
%     title(str);
%     view([-112,13]);
    pause(1);
end
end

%%% The nonlinearity
function V = h(u)
V = u^2;
end

function Ce = efk_matrix(centroid, diameter, projector, projection, nsides, verts)
modwrap = @(x,a) mod(x-1,a) + 1;
Q = quadrature_rule(6,2);
qw = Q(:,1);    qx = Q(:,2);    qy = Q(:,3);
H = zeros(6);
for alpha = 1:6
    for beta = 1:6
        for v=1:nsides
            vert = verts(v,:);
            next = verts(modwrap(v-1,nsides),:);
            area = 0.5*abs(det([1, centroid; 1, vert; 1, next]));
            for q = 1:length(qw)
                xhat = (vert(1)-centroid(1))*qx(q) + (next(1)-centroid(1))*qy(q) + centroid(1);
                yhat = (vert(2)-centroid(2))*qx(q) + (next(2)-centroid(2))*qy(q) + centroid(2);
                M = [1; (xhat-centroid(1))/diameter; (yhat-centroid(2))/diameter;
                    ((xhat-centroid(1))/diameter)^2;
                    ((xhat-centroid(1))/diameter)*((yhat-centroid(2))/diameter);
                    ((yhat-centroid(2))/diameter)^2];                
                proj = dot(M,projection);
                H(alpha, beta) = H(alpha, beta) + 2*area*qw(q)*(M(beta)*M(alpha))*h(proj);                
            end
        end
    end
end

Ce = projector'*(H)*projector;
end