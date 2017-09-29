function timedependent
clc
clear
format short
close all

mesh = load('squares_linear25');
t0 = 0;  tf = 0.5;  delt = 0.1;

u_final = heat(mesh, @f, @g_D, @u0, tf, t0, delt);

plot_solution(mesh, u_final);

end

function v = u0(x,y)
v = sin(pi*x).*sin(pi*y);
end

function v = f(x,y,t)
v = exp(-t)*(2*pi^2-1)*sin(pi*x).*sin(pi*y);
end

function v = g_D(x,y,t)
v = 0.*x.*y*t;
end