function [Ke,Fe,Me,G,D,B,H] = elem_matrices(mesh,id,poly_deg,f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function that returns the Elemental Data:
%       Input :
%               mesh, element_id, poly_degree_used,
%               Right hand side function : f (Provide arbitrary argument if
%               not used.
%       Output :
%               Ke : Elemental Stiffness Matrix
%               Fe : Elemental Load Vector
%               Me : Elemental Mass Matrix
%               Matrices G,D,B,H - In the same order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% BEGIN Constructing matrices for k = 1
if(poly_deg==1)
    elem = mesh.elements{id};
    verts = mesh.vertices(elem,:);
    nsides = length(elem);
    npolys = 3;    
    [area, centroid, diameter] = geo(mesh,id);
    modwrap = @(x,a) mod(x-1,a) + 1;    
    D = zeros(nsides, npolys);
    D(:, 1) = 1;
    B = zeros(npolys, nsides);
    B(1, :) = 1/nsides;    
    linear_polynomials = {[0,0], [1,0], [0,1]};
    for v = 1:nsides
        vert = verts(v, :); % This vertex and its neighbours
        prev = verts(modwrap(v-1, nsides), :);
        next = verts(modwrap(v+1, nsides), :);
        normal = [next(2)-prev(2), prev(1)-next(1)]; % Average of normals on edges
        for poly_id = 2:npolys % Only need to loop over non-constant polynomials
            poly_degree = linear_polynomials{poly_id};
            monomial_grad = poly_degree / diameter; % Gradient of a linear polynomial is constant
            D(v, poly_id) = dot((vert - centroid), poly_degree) / diameter;
            B(poly_id, v) = 0.5 * dot(monomial_grad, normal);
        end
    end
    %%%% Routine to compute H
    H = zeros(3,3);
    H(1,1) = area;
    for j=2:3
        for k=2:3
            for v = 1:nsides
                vert = verts(v,:);
                next = verts(modwrap(v-1, nsides), :);                
                H(j,k) = H(j,k) + integrate2dtri({centroid,vert,next},...
                    linear_polynomials{j},linear_polynomials{k},diameter,1);
            end
        end
    end    
    projector = (B*D)\B;
    stabilising_term = (eye(nsides) - D * projector)' * (eye(nsides) - D * projector);
    G = B*D;
    Gt = G;
    Gt(1, :) = 0;
    Ke = projector'*Gt*projector + stabilising_term;
    Fe = ones(nsides,1)*f(centroid(1),centroid(2))*area/nsides;    
    %%%%%% Mass Matrix
    C = H*(G\B);
    ritz_projector = D*projector;
    Me = C'*(H\C) + area*(eye(nsides)-ritz_projector)'*(eye(nsides)-ritz_projector);
end
%%%% END  Constructing matrices for k = 1
%%%% BEGIN Constructing matrices for k = 2
if(poly_deg == 2)
    elem = mesh.elements{id};
    eldofs = length(elem);
    verts = mesh.vertices(elem(1:(eldofs-1)/2),:); % Vertices
    nsides = length(verts);
    v = mesh.vertices(mesh.elements{id}(1:end-1),:); % Vertices + Midpoints    
    npolys = 6; % Number of scaled monomials
    quad_polynomials = {[0,0],[1,0],[0,1],[2,0],[1,1],[0,2]}; % Description
    [area, centroid, diameter] = geo(mesh,id);
    modwrap = @(x,a) mod(x-1,a) + 1;
    %%%% Find out D matrix (DOF Matrix)
    % 1 - (Last-1) DOFs (Values)
    D = zeros(eldofs,npolys);
    for i=1:eldofs-1
        vert = v(i,:);
        for alpha = 1:npolys
            poly_degree = quad_polynomials{alpha};
            D(i,alpha) = ((vert(1)-centroid(1))/diameter)^(poly_degree(1))*...
                ((vert(2)-centroid(2))/diameter)^(poly_degree(2));
        end
    end
    % Last DOF (Moments)
    for alpha = 1:npolys
        for k = 1:length(verts)
            vert = verts(k,:);
            next = verts(modwrap(k+1, nsides), :);
            D(eldofs,alpha) = D(eldofs,alpha) + ...
                (1/area)*integrate2dtri({centroid, vert, next},...
                [0,0], quad_polynomials{alpha}, diameter, 6);
        end
    end    
    %%%% Find out B matrix
    B = zeros(npolys, eldofs);
    for i=1:eldofs-1
        if(i <= nsides)            
            vert = verts(i,:);
            next = verts(modwrap(i+1, nsides), :);
            prev = verts(modwrap(i-1, nsides), :);             
            normal1 = [next(2) - vert(2), vert(1) - next(1)];
            normal2 = [vert(2) - prev(2), prev(1) - vert(1)];
            normal1 = normal1/norm(normal1);
            normal2 = normal2/norm(normal2);
            for alpha = 1:npolys
                B(alpha,i) = ...
                    integrate1dline(prev, vert, alpha, 3, {centroid, diameter, normal2}) + ...
                    integrate1dline(vert, next, alpha, 1, {centroid, diameter, normal1});
            end
        elseif(i > nsides)     
            vert = verts(i-nsides,:);
            next = verts(modwrap(i-nsides+1, nsides), :);            
            normal = [next(2) - vert(2), vert(1) - next(1)];            
            normal = normal/norm(normal);
            for alpha = 1:npolys
                B(alpha,i) = integrate1dline(vert, next, alpha, 2, {centroid, diameter, normal});
            end
        end        
    end
    lapm = [0;0;0;2/diameter^2;0;2/diameter^2];
    B(:,eldofs) = -area*lapm;
    B(1,eldofs) = 1;    
    G = B*D;
    projector = G\B;
    stabilising_term = (eye(eldofs) - D * projector)' * (eye(eldofs) - D * projector);
    Gt = G;
    Gt(1, :) = 0;
    Ke = projector'*Gt*projector + stabilising_term;
    Fe = zeros(size(Ke,1),1);
    for i=1:nsides
        vert = v(i,:);
        next = v(modwrap(i-1, nsides), :);
        Fe = Fe + loadvec(projector, centroid, vert, next, diameter, f);
    end    
    %%%%% Temporary variables - not assigned yet
    Me = Ke;
    H = G;
end
end

function V = integrate2dtri(C,alpha,beta,diameter,qrule)
if(qrule==1)
    P1 = C{1};  P2 = C{2};  P3 = C{3};
    centroid = P1;
    tricentroid = 1/3*(P1+P2+P3);
    A = [1, P1; 1, P2; 1, P3];
    area = 0.5*abs(det(A));
    ma = @(x,y)((x-centroid(1))/diameter)^(alpha(1))*((y-centroid(2))/diameter)^(alpha(2));
    mb = @(x,y)((x-centroid(1))/diameter)^(beta(1))*((y-centroid(2))/diameter)^(beta(2));
    V = area*ma(tricentroid(1),tricentroid(2))*mb(tricentroid(1),tricentroid(2));
elseif(qrule==6)
    P1 = C{1};   P2 = C{2};   P3 = C{3};
    centroid = P1;
    area = 0.5*abs(det([1, P1; 1, P2; 1, P3]));    
    ma = @(x,y)((x-centroid(1))/diameter)^(alpha(1))*((y-centroid(2))/diameter)^(alpha(2));
    mb = @(x,y)((x-centroid(1))/diameter)^(beta(1))*((y-centroid(2))/diameter)^(beta(2));    
    Q = quadrature_rule(6,2); % 6th order triangular quadrature;
    qw = Q(:,1);  qx = Q(:,2); qy = Q(:,3); % Obtain Quadrature points and weights    
    V = 0;
    for q = 1:qrule
        xhat = (P2(1)-P1(1))*qx(q) + (P3(1)-P1(1))*qy(q) + P1(1);
        yhat = (P2(2)-P1(2))*qx(q) + (P3(2)-P1(2))*qy(q) + P1(2);        
        V = V + qw(q)*2*area*ma(xhat,yhat)*mb(xhat,yhat);
    end
end
end

function V = integrate1dline(p1, p2, alpha, basis, linegeo)
x1 = p1(1); x2 = p2(1); y1 = p1(2); y2 = p2(2);
dx = norm(p2-p1);
Q = quadrature_rule(2,1); % 2 Point Quadrature rule on a line  
qx = Q(:,2); qw = Q(:,1); % Get points and weights
qG = length(qx);
centroid = linegeo{1}; diameter = linegeo{2}; normal = linegeo{3};
V = 0;
for q = 1:qG
    xhat = (x2-x1)/2*qx(q) + (x2+x1)/2;
    yhat = (1+qx(q))/2*(y2-y1) + y1;
    Mx = [0; 1/diameter; 0; (2/diameter^2)*(xhat-centroid(1)); ...
        (1/diameter^2)*(yhat-centroid(2)); 0];
    My = [0; 0; 1/diameter; 0; (1/diameter^2)*(xhat-centroid(1)); ...
        (2/diameter^2)*(yhat-centroid(2))];
    dmdn = Mx(alpha)*normal(1) + My(alpha)*normal(2);    
    phi = [-0.5*qx(q)+0.5*qx(q)^2, 1-qx(q)^2, 0.5*qx(q)+0.5*qx(q)^2];
    V = V + qw(q)*(dmdn*phi(basis))*0.5*dx;
end

end

function V = loadvec(projector, centroid, vert, next, diameter, f)
P1 = centroid; P2 = vert; P3 = next;
A = [1, P1; 1, P2; 1, P3];
area = 0.5*abs(det(A));
Q = quadrature_rule(6,2); % 6th order triangular quadrature;
qw = Q(:,1);  qx = Q(:,2); qy = Q(:,3); % Obtain weights and quadrature points
qG = length(qw);
ndofs = size(projector,2);
V = zeros(ndofs,1);
for i=1:ndofs
    for q = 1:qG
        xhat = (P2(1)-P1(1))*qx(q) + (P3(1)-P1(1))*qy(q) + P1(1);
        yhat = (P2(2)-P1(2))*qx(q) + (P3(2)-P1(2))*qy(q) + P1(2);        
        m = [1; (xhat-centroid(1))/diameter; (yhat-centroid(2))/diameter;
            ((xhat-centroid(1))/diameter)^2;
            ((xhat-centroid(1))/diameter)*((yhat-centroid(2))/diameter);
            ((yhat-centroid(2))/diameter)^2];
        projector_i = projector'*m;
        V(i) = V(i) + 2*area*qw(q)*f(xhat,yhat)*projector_i(i,:);
    end
end
end