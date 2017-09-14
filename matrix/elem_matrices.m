function [Ke,Fe,Me,G,D,B,H] = elem_matrices(mesh,id,poly_deg,f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function that returns the Stiffness matrix
%

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
               
               H(j,k) = H(j,k) + integrate_linear({centroid,vert,next},...
                   linear_polynomials{j},linear_polynomials{k},diameter);
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

end

function V = integrate_linear(C,alpha,beta,diameter)
P1 = C{1};
P2 = C{2};
P3 = C{3};
centroid = P1;
tricentroid = 1/3*(P1+P2+P3);
A = [1, P1; 1, P2; 1, P3];
area = 0.5*abs(det(A));
ma = @(x,y)(1/diameter^2)*(x-centroid(1))^(alpha(1))*(y-centroid(2))^(alpha(2));
mb = @(x,y)(1/diameter^2)*(x-centroid(1))^(beta(1))*(y-centroid(2))^(beta(2));
V = area*ma(tricentroid(1),tricentroid(2))*mb(tricentroid(1),tricentroid(2));
end