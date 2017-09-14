function V = gradprojector(u,id,mesh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to construct the gradient L^2-projection
%           pi-delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

element = mesh.elements{id};
nsides = size(element,1);
verts = mesh.vertices(element,:);

[~,~,~,G,~,B,~] = elem_matrices(mesh,id,1,u);

dofu = zeros(nsides,1);
for i=1:nsides
    vert = verts(i,:);
    dofu(i) = u(vert(1), vert(2));
end
b = B*dofu;
V = G\b;

end