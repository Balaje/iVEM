function squares
clear
clc
close all
format short

mesh1 = load('squares_linear10.mat');
mesh2 = load('squares_linear15.mat');
mesh3 = load('squares_linear20.mat');
mesh4 = load('squares_linear25.mat');
mesh5 = load('squares_linear30.mat');
mesh6 = load('squares_linear35.mat');
mesh7 = load('squares_linear40.mat');


mesh = {mesh1, mesh2, mesh3, mesh4, mesh5, mesh6, mesh7};
nref = length(mesh);

l2err = zeros(nref,1);
h1err = zeros(nref,1);
N = zeros(nref,1);
l2order = zeros(nref-1,1);
h1order = zeros(nref-1,1);

for it = 1:nref
    usol = poisson(mesh{it}, @f, @g_D);
    
    [l2err(it),h1err(it)] = l2error(mesh{it},@exact,usol);
    N(it) = size(mesh{it}.elements,1);
    
    figure(1)
    subplot(1,2,1)
    cla
    plot_solution(mesh{it},usol);
    title('Solution of poisson equation','FontSize',14);
    grid on
    view([-47,10])
    subplot(1,2,2)
    cla
    plot_solution(mesh{it},zeros(length(usol),1));
    title('Virtual Element Mesh','FontSize',14);
    
    pause(0.01);
    hold off
end

figure(2)
subplot(1,2,1)
loglog(1./N.^(0.5), l2err, 'b*-', 1./N.^0.5, 1./N, 'k-');
grid on
title('L^2 order of convergence','FontSize',14);
xlabel('Mesh size','FontSize',14);
ylabel('Error','FontSize',14);
legend('||e||_{L^2}','Line of slope = 2')
subplot(1,2,2)
loglog(1./N.^(0.5), h1err, 'ro-', 1./N.^0.5, 1./N.^0.5, 'k-');
grid on
title('H^1 order of convergence','FontSize',14);
xlabel('Mesh size','FontSize',14);
ylabel('Error','FontSize',14);
legend('||e||_{H^1}','Line of slope = 1')


for it = 1:nref-1
    l2order(it) = log(l2err(it)/l2err(it+1))/(0.5*log(N(it+1)/N(it)));
    h1order(it) = log(h1err(it)/h1err(it+1))/(0.5*log(N(it+1)/N(it)));
end

fprintf('L^2 order and H^1 order: \n');
disp([l2order, h1order]);

end

function v = f(x,y)
v = 2*pi^2*sin(pi*x).*sin(pi*y);
end

function v = exact(x,y)
v = sin(pi*x).*sin(pi*y);
end

function v = g_D(x,y)
v = 0*x.*y;
end