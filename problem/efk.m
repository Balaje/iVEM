function [u,v] = efk(mesh, f, g_Du, g_Dv, u0, tf, t0, delt)

ndof = length(mesh.vertices);
nel = length(mesh.elements);

ntimes = fix((tf-t0)/delt);

un = u0(mesh.vertices(:,1), mesh.vertices(:,2));
uu = ones(length(un),1);
v = zeros(length(un),1);

t = 0;

for time = 1:ntimes
    %%% Begin Tolerance loop
    error = 100;
    tol = 1e-7;
    
    t = t + delt;
    while error > tol
        K = zeros(ndof,ndof);
        M = zeros(ndof,ndof);
        C = zeros(ndof,ndof);
        F = zeros(ndof,1);
        
        for el = 1:nel
            g = @(x,y) f(x,y,t);
            
            [~, centroid, diameter] = geo(mesh, el);
            [Ke, Fe, Me, G, ~, B, ~] = elem_matrices(mesh,el,1,g);
            
            K(mesh.elements{el},mesh.elements{el}) = ...
                K(mesh.elements{el},mesh.elements{el}) + Ke;
            
            M(mesh.elements{el},mesh.elements{el}) = ...
                M(mesh.elements{el},mesh.elements{el}) + Me;
            
            F(mesh.elements{el},1) = F(mesh.elements{el},1) + Fe;
            %%%%% The nonlinear term
            nsides = length(mesh.elements{el});
            verts = mesh.vertices(mesh.elements{el},:);
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
    
    cla
    plot_solution(mesh,u)
    grid on
    str = ['Solution at time t = ',num2str(t)];
    title(str);
    view([-112,13]);
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
H = zeros(3,3);
for alpha = 1:3
    for beta = 1:3
        for v=1:nsides
            vert = verts(v,:);
            next = verts(modwrap(v-1,nsides),:);
            area = 0.5*abs(det([1, centroid; 1, vert; 1, next]));
            for q = 1:length(qw)
                xhat = (vert(1)-centroid(1))*qx(q) + (next(1)-centroid(1))*qy(q) + centroid(1);
                yhat = (vert(2)-centroid(2))*qx(q) + (next(2)-centroid(2))*qy(q) + centroid(2);
                M = [1; (xhat-centroid(1))/diameter; (yhat-centroid(2))/diameter];
                proj = dot(M,projection);                
                H(alpha, beta) = H(alpha, beta) + 2*area*qw(q)*(M(beta)*M(alpha))*h(proj);                
            end
        end
    end
end

Ce = projector'*(H)*projector;
end