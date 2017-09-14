function u = heat(mesh,f,g_Dp,u0,delt)
ndof = length(mesh.vertices);
nel = length(mesh.elements);

K = sparse(ndof,ndof);
M = sparse(ndof,ndof);
F = sparse(ndof,1);

for el = 1:nel
    [Ke, Fe, Me, ~, ~, ~] = elem_matrices(mesh,el,1,f);
    K(mesh.elements{el},mesh.elements{el}) = ...
        K(mesh.elements{el},mesh.elements{el}) + Ke;
    
    M(mesh.elements{el},mesh.elements{el}) = ...
        M(mesh.elements{el},mesh.elements{el}) + Me;   
        
    F(mesh.elements{el},1) = F(mesh.elements{el},1) + Fe;
end

boundaryvals_p = g_Dp(mesh.vertices(mesh.boundary,1),...
                   mesh.vertices(mesh.boundary,2));
                              
F = delt*F + M*u0...
    - (M(:,mesh.boundary) + delt*K(:,mesh.boundary))*boundaryvals_p;

freenodes = mesh.boundary;
solnodes = setdiff(1:ndof, freenodes);

u = zeros(ndof,1);
u(solnodes) = (M(solnodes,solnodes) + delt*K(solnodes,solnodes))\F(solnodes);
u(freenodes) = g_Dp(mesh.vertices(mesh.boundary,1),...
                   mesh.vertices(mesh.boundary,2));
end