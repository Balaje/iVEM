function V = l2projector(u,id,mesh,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to construct the standard L^2-projection
%           pi-0
% --- For k=1 
%           pi-delta = pi-0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = gradprojector(u,id,mesh,k);
end