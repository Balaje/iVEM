function simple
clc
clear
%close all
format short

mesh = load('../mesh/voronoi.mat');
usol = poisson(mesh, @f, @g_D);

exact = @(x,y)sin(pi*x)*sin(pi*y);
[l2err,h1err] = l2error(mesh,exact,usol);

figure(3)
subplot(1,2,1);
plot_solution(mesh,usol);
grid on
str = 'Solution of Poisson equation';
title(str,'FontSize',14,'interpreter','tex');
view([-55,14])
txt = ['||e||_{L^2} = ',num2str(l2err); '||e||_{H^1} = ',num2str(h1err)];
text(1,1,1,txt,'FontSize',11)
subplot(1,2,2);
plot_solution(mesh,zeros(length(usol),1));
title('Virtual Element Mesh','FontSize',14);
end

function v = f(x,y)
v = 2*pi^2*sin(pi*x).*sin(pi*y);
end

function v = g_D(x,y)
v = 0*x.*y;
end