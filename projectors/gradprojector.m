function V = gradprojector(u,id,mesh,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to construct the gradient L^2-projection
%           pi-delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(k==1)
    element = mesh.elements{id};
    nsides = size(element,1);
    verts = mesh.vertices(element,:);    
    [~,~,~,G,~,B,~] = elem_matrices(mesh,id,k,u);    
    dofu = zeros(nsides,1);
    for i=1:nsides
        vert = verts(i,:);
        dofu(i) = u(vert(1), vert(2));
    end
    b = B*dofu;
    V = G\b;
elseif(k==2)    
    Q = quadrature_rule(6,2);
    qw = Q(:,1); qx = Q(:,2); qy = Q(:,3);
    [area,centroid,~] = geo(mesh,id);
    modwrap = @(x,a) mod(x-1,a) + 1;
    element = mesh.elements{id};
    eldofs = size(element,1);
    nsides = (eldofs-1)/2;
    verts = mesh.vertices(element(1:end-1),:);    
    [~,~,~,G,~,B,~] = elem_matrices(mesh,id,k,u);
    dofu = zeros(eldofs,1);    
    for i=1:eldofs-1
        vert = verts(i,:);
        dofu(i) = u(vert(1), vert(2));
    end
    dofu(eldofs) = 0;
    for v=1:nsides           
        vert = verts(v,:);
        next = verts(modwrap(v+1, nsides),:);
        A = [1, centroid; 1, vert; 1, next];
        triarea = 0.5*abs(det(A));
        for q=1:length(qw)
            xhat = (vert(1)-centroid(1))*qx(q) + (next(1)-centroid(1))*qy(q) + centroid(1);
            yhat = (vert(2)-centroid(2))*qx(q) + (next(2)-centroid(2))*qy(q) + centroid(2);
            dofu(eldofs) = dofu(eldofs) + 2*triarea*qw(q)*u(xhat,yhat);
        end
    end
    dofu(eldofs) = (1/area)*dofu(eldofs);
    b = B*dofu;
    V = G\b;
end

end