function timedependent
clc
clear
format short
close all

syms x y

mesh = load('../mesh/non-convex.mat');
%%% Discretizing time
t0 = 0;
tf = 1;
delt = 0.1;
u0 = init(mesh.vertices(:,1), mesh.vertices(:,2));

N = fix((tf-t0)/delt);

t = 0;
for m=1:N
    v1 = f(x,y,t+delt);
    v2 = matlabFunction(v1);
    v3 = g_D(x,y,t+delt);
    v4 = matlabFunction(v3);
    
    u = heat(mesh,v2,v4,u0,delt);   
    
    figure(1)    
    subplot(1,2,1);
    cla
    plot_solution(mesh,u);
    grid on
    str = ['Solution of Heat equation at t = ', num2str(t)];
    title(str,'FontSize',14,'interpreter','tex');
    view([-55,14])
    subplot(1,2,2);
    plot_solution(mesh,zeros(length(u),1));
    title('Virtual Element Mesh','FontSize',14);
    
    pause(0.01)
    
    t = t+delt;
    u0 = u;
end

% syms x y
% v = f(x,y,1);
% v1 = matlabFunction(v)
end

function v = init(x,y)
v = sin(pi*x).*sin(pi*y);
end

function v = f(x,y,t)
v = exp(-t)*(2*pi^2-1)*sin(pi*x).*sin(pi*y);
end

function v = g_D(x,y,t)
v = sin(pi*x).*sin(pi*y)*exp(-t);
end