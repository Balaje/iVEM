function voronoi_2
clear
clc
close all
format short

mesh1 = load('voronoi_quadratic100.mat');
mesh2 = load('voronoi_quadratic225.mat');
mesh3 = load('voronoi_quadratic400.mat');
mesh4 = load('voronoi_quadratic625.mat');
mesh5 = load('voronoi_quadratic900.mat');

plotmesh1 = load('voronoi100.mat');
plotmesh2 = load('voronoi225.mat');
plotmesh3 = load('voronoi400.mat');
plotmesh4 = load('voronoi625.mat');
plotmesh5 = load('voronoi900.mat');

mesh = {mesh1, mesh2, mesh3, mesh4, mesh5};
plotmesh = {plotmesh1, plotmesh2, plotmesh3, plotmesh4, plotmesh5};
nref = length(mesh);

l2err = zeros(nref,1);
h1err = zeros(nref,1);
N = zeros(nref,1);
l2order = zeros(nref-1,1);
h1order = zeros(nref-1,1);

for it = 1:nref
    usol = poisson_2(mesh{it}, @f, @g_D);
    
    [l2err(it),h1err(it)] = l2error(mesh{it},@exact,usol,2);
    N(it) = size(mesh{it}.elements,1);
    
    noofvertices = length(plotmesh{it}.vertices);
    u_verts = usol(1:noofvertices);
    
    figure(1)
    subplot(1,2,1)
    cla
    plot_solution(plotmesh{it},u_verts);
    title('Solution of poisson equation','FontSize',14);
    grid on
    view([-112,13])
    subplot(1,2,2)
    cla
    plot_solution(plotmesh{it},zeros(length(u_verts),1));
    title('Virtual Element Mesh','FontSize',14);
    
    pause(0.01);
    hold off
end

figure(2)
subplot(1,2,1)
loglog(1./N.^(0.5), l2err, 'b*-', 1./N.^0.5, 1./N.^1.5, 'k-');
grid on
title('L^2 order of convergence','FontSize',14);
xlabel('Mesh size','FontSize',14);
ylabel('Error','FontSize',14);
legend('||e||_{L^2}','Line of slope = 3')
subplot(1,2,2)
loglog(1./N.^(0.5), h1err, 'ro-', 1./N.^0.5, 1./N.^1, 'k-');
grid on
title('H^1 order of convergence','FontSize',14);
xlabel('Mesh size','FontSize',14);
ylabel('Error','FontSize',14);
legend('||e||_{H^1}','Line of slope = 2')


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