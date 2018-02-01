function [u,v] = rosenau_2(mesh, f, g_D, u0, tf, t0, delt, plotmesh)

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
            Ce = conv_matrix(centroid, diameter, projector, projection, nsides, verts);
            
            C(mesh.elements{el},mesh.elements{el}) = ...
                C(mesh.elements{el},mesh.elements{el}) + Ce;
        end
        
        boundaryvals_p = g_D(mesh.vertices(mesh.boundary,1),...
            mesh.vertices(mesh.boundary,2), t);
        
        freenodes = mesh.boundary;
        solnodes = setdiff(1:ndof, freenodes);
        
        vn(mesh.boundary) = boundaryvals_p;
        vn(solnodes) = (M(solnodes,solnodes)\K(solnodes,solnodes))*un(solnodes);
        
        F = delt*F + (M*un + K*vn)...
            - (M(:,mesh.boundary) + K(:,mesh.boundary)*(M(:,mesh.boundary)\K(:,mesh.boundary))...
                - delt*C(:,mesh.boundary))*boundaryvals_p;               
        
        u = zeros(ndof,1);
        u(solnodes) = (M(solnodes,solnodes) + ...
            K(solnodes,solnodes)*(M(solnodes,solnodes)\K(solnodes,solnodes))...
            - delt*C(solnodes,solnodes))\F(solnodes);
        
        error = max(abs(uu - u));
        uu = u;
    end
    
    un = u;
    boundaryvals_p = g_D(mesh.vertices(mesh.boundary,1),...
            mesh.vertices(mesh.boundary,2), t);
        
    freenodes = mesh.boundary;
    solnodes = setdiff(1:ndof, freenodes);
    
    v(mesh.boundary) = boundaryvals_p;
    v(solnodes) = (M(solnodes,solnodes)\K(solnodes,solnodes))*u(solnodes);
    
%     cla
%     plot_solution(plotmesh,u(1:length(plotmesh.vertices)));
%     grid on
    str = ['Solution at time t = ',num2str(t)];
%     title(str);
%     view([-112,13]);
    disp(str);
    pause(1);
end
end

%%% The nonlinearity
function V = h(u)
V = 2*u;
end

function Ce = conv_matrix(centroid, diameter, projector, projection, nsides, verts)
modwrap = @(x,a) mod(x-1,a) + 1;
Q = quadrature_rule(6,2);
qw = Q(:,1);    qx = Q(:,2);    qy = Q(:,3);
Hx = zeros(6);    Hy = zeros(6);
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
                Mx = [0; 1/diameter; 0; (2/diameter^2)*(xhat-centroid(1)); ...
                    (1/diameter^2)*(yhat-centroid(2)); 0];
                My = [0; 0; 1/diameter; 0; (1/diameter^2)*(xhat-centroid(1)); ...
                    (2/diameter^2)*(yhat-centroid(2))];
                proj = dot(M,projection);
                Hx(alpha, beta) = Hx(alpha, beta) + 2*area*qw(q)*(Mx(beta)*M(alpha))*h(proj);
                Hy(alpha, beta) = Hy(alpha, beta) + 2*area*qw(q)*(My(beta)*M(alpha))*h(proj);
            end
        end
    end
end

Ce = projector'*(Hx+Hy)*projector;
end