function simple_2
clear
close all
format long

mesh = load('voronoi_quadratic');
plotmesh = load('voronoi'); % Must input the corresponding mesh only.

usol = poisson_2(mesh, @f, @g_D);

exact = @(x,y)sin(pi*x)*sin(pi*y);
[l2err,h1err] = l2error(mesh,exact,usol,2);

%%% To plot, we select the solution points on the vertices only
noofvertices = length(plotmesh.vertices);
u_verts = usol(1:noofvertices);

figure(1)
subplot(1,2,1);
plot_solution(plotmesh,u_verts);
grid on
str = 'Solution of Poisson equation';
title(str,'FontSize',14,'interpreter','tex');
view([-112,13])
l2 = num2str(double(l2err));
h1 = num2str(double(h1err));
txt = ['||e||_{L^2} = ',l2,'   |e|_{H^1} = ',h1];
text(1,1,1,txt,'FontSize',11)
subplot(1,2,2);
plot_solution(plotmesh,zeros(length(u_verts),1));
hold on
scatter(mesh.vertices(:,1),mesh.vertices(:,2),5,'rs');
[~,centroid,~] = geo(mesh);
for i=1:length(centroid)
    scatter(centroid{i}(1), centroid{i}(2),4,'b+')
end
title('Virtual Element Mesh','FontSize',12);
hold off

end

function v = f(x,y)
v = 2*pi^2*sin(pi*x).*sin(pi*y);
end

function v = g_D(x,y)
v = 0*x.*y;
end