function gradrecovery2
clc
clear
close all
format long

mesh = load('squares_quadratic30.mat');
plotmesh = load('squares_linear30.mat');

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
text(1,1,max(usol)+0.05,txt,'FontSize',11)
subplot(1,2,2);
plot_solution(plotmesh,zeros(length(u_verts),1));
hold on
scatter(mesh.vertices(:,1),mesh.vertices(:,2),5,'rs');
[~,centroid,~] = geo(mesh);
for i=1:length(centroid)
    scatter(centroid{i}(1), centroid{i}(2),4,'b+')
end
title('Virtual Element Mesh','FontSize',14);
hold off

[Mx, My, M] = gradrecovery_matrix(mesh,2);

D = diag(1./sum(M,1));
uhx = D*Mx*usol;
uhy = D*My*usol;

[l2err_x, ~] = l2error(mesh, @(x,y) pi*cos(pi*x)*sin(pi*y), uhx, 2);
[l2err_y, ~] = l2error(mesh, @(x,y) pi*sin(pi*x)*cos(pi*y), uhy, 2);

%fprintf("L2 Error = (%d,\t%d)\n",l2err_x,l2err_y);

UHX = uhx(1:noofvertices);
UHY = uhy(1:noofvertices);

figure(2)
subplot(1,2,1);
plot_solution(plotmesh,UHX);
grid on
str = 'Recovered x-gradient';
l2 = num2str(l2err_x);
txt = ['||e||_{L^2} = ',l2];
text(1,1,3,txt,'FontSize',11)
title(str,'FontSize',14,'interpreter','tex');
view([-55,14])
subplot(1,2,2);
plot_solution(plotmesh,UHY);
grid on
str = 'Recovered y-gradient';
l2 = num2str(l2err_y);
txt = ['||e||_{L^2} = ',l2];
text(1,1,3,txt,'FontSize',11)
title(str,'FontSize',14,'interpreter','tex');
view([-55,14])


end

function v = f(x,y)
v = 2*pi^2*sin(pi*x).*sin(pi*y);
end

function v = g_D(x,y)
v = 0*x.*y;
end