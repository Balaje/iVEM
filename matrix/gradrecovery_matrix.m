function [Mx, My, Q] = gradrecovery_matrix(mesh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that returns the Gradient Recovery Matrices  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nel = size(mesh.elements,1);
nodes = size(mesh.vertices,1);

[~, ~, diameter] = geo(mesh);

Mx = zeros(nodes);
My = zeros(nodes);
Q = zeros(nodes);

for el = 1:nel
    nsides = length(mesh.elements{el});
    %usol = uh(mesh.elements{el});
       
    [~,~,Me,G,~,B,~] = elem_matrices(mesh,el,1,@(x,y) 1);    
    Qe = diag(sum(Me,1));

    Dx = repmat([0, 1/diameter{el}, 0], nsides, 1);
    Dy = repmat([0, 0, 1/diameter{el}], nsides, 1);
    pi0uh = (G\B);
    
    Mxe = Qe*Dx*pi0uh;
    Mye = Qe*Dy*pi0uh;
        
    Mx(mesh.elements{el},mesh.elements{el}) = ...
        Mx(mesh.elements{el},mesh.elements{el}) + Mxe;
    My(mesh.elements{el},mesh.elements{el}) = ...
        My(mesh.elements{el},mesh.elements{el}) + Mye;
    
    Q(mesh.elements{el},mesh.elements{el}) = ...
        Q(mesh.elements{el},mesh.elements{el}) + Me;
end


end