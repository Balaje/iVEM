function u = heat(mesh, f, g_D, u0, tf, t0, delt)
ndof = length(mesh.vertices);
nel = length(mesh.elements);

ntimes = fix((tf-t0)/delt);

un = u0(mesh.vertices(:,1), mesh.vertices(:,2));

t = 0;

for time = 1:ntimes
    
    K = sparse(ndof,ndof);
    M = sparse(ndof,ndof);
    F = sparse(ndof,1);
    
    for el = 1:nel
        g = @(x,y) f(x,y,t);        
        [Ke, Fe, Me, ~, ~, ~] = elem_matrices(mesh,el,1,g);
        
        K(mesh.elements{el},mesh.elements{el}) = ...
            K(mesh.elements{el},mesh.elements{el}) + Ke;
        
        M(mesh.elements{el},mesh.elements{el}) = ...
            M(mesh.elements{el},mesh.elements{el}) + Me;
        
        F(mesh.elements{el},1) = F(mesh.elements{el},1) + Fe;
    end
    
    boundaryvals_p = g_D(mesh.vertices(mesh.boundary,1),...
        mesh.vertices(mesh.boundary,2), t);
    
    F = delt*F + M*un...
        - (M(:,mesh.boundary) + delt*K(:,mesh.boundary))*boundaryvals_p;
    
    freenodes = mesh.boundary;
    solnodes = setdiff(1:ndof, freenodes);
    
    u = zeros(ndof,1);
    u(solnodes) = (M(solnodes,solnodes) + delt*K(solnodes,solnodes))\F(solnodes);
    u(freenodes) = g_D(mesh.vertices(mesh.boundary,1),...
        mesh.vertices(mesh.boundary,2), t);
    
    un = u;
    t = t + delt;
end

end