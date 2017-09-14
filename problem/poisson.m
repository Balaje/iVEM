function u = poisson(mesh,f,g_D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to solve the Poisson problem using the Virtual Element Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Get all the geometric data of the mesh
ndof = length(mesh.vertices);
nel = length(mesh.elements);

K = sparse(ndof,ndof);
F = sparse(ndof,1);

for el = 1:nel
    [Ke, Fe, ~, ~, ~, ~, ~] = elem_matrices(mesh,el,1,f);
    K(mesh.elements{el},mesh.elements{el}) = ...
        K(mesh.elements{el},mesh.elements{el}) + Ke;
        
    F(mesh.elements{el},1) = F(mesh.elements{el},1) + Fe;
end

boundaryvals = g_D(mesh.vertices(mesh.boundary,1),...
                   mesh.vertices(mesh.boundary,2));
F = F - K(:,mesh.boundary)*boundaryvals;

freenodes = mesh.boundary;
solnodes = setdiff(1:ndof, freenodes);

u = zeros(ndof,1);
u(solnodes) = K(solnodes,solnodes)\F(solnodes);
u(freenodes) = g_D(mesh.vertices(mesh.boundary,1),...
                   mesh.vertices(mesh.boundary,2));

end