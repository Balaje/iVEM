function gradrecovery
clc
clear
close all
format long

mesh = load('squares_linear10.mat');
usol = poisson(mesh, @f, @g_D);

exact = @(x,y)sin(pi*x)*sin(pi*y);
[l2err,h1err] = l2error(mesh,exact,usol,1,false);

figure(1)
subplot(1,2,1);
plot_solution(mesh,usol);
grid on
str = 'Solution of Poisson equation';
title(str,'FontSize',14,'interpreter','tex');
view([-112,13])
l2 = num2str(l2err);
h1 = num2str(h1err);
txt = ['||e||_{L^2} = ',l2,'   |e|_{H^1} = ',h1];
text(1,1,max(usol)+0.05,txt,'FontSize',11)
subplot(1,2,2);
plot_solution(mesh,zeros(length(usol),1));
title('Virtual Element Mesh','FontSize',14);

[Mx, My, M] = gradrecovery_matrix(mesh,1);

D = diag(1./sum(M,1));
uhx = D*Mx*usol;
uhy = D*My*usol;

[l2err_x, ~] = l2error(mesh, @(x,y) pi*cos(pi*x)*sin(pi*y), uhx, 1);
[l2err_y, ~] = l2error(mesh, @(x,y) pi*sin(pi*x)*cos(pi*y), uhy, 1);

figure(2)
subplot(1,2,1);
plot_solution(mesh,uhx);
grid on
str = 'Recovered x-gradient';
l2 = num2str(l2err_x);
txt = ['||e||_{L^2} = ',l2];
text(1,1,3,txt,'FontSize',11)
title(str,'FontSize',14,'interpreter','tex');
view([-55,14])
subplot(1,2,2);
plot_solution(mesh,uhy);
grid on
str = 'Recovered y-gradient';
l2 = num2str(l2err_y);
txt = ['||e||_{L^2} = ',l2];
text(1,1,3,txt,'FontSize',12)
title(str,'FontSize',14,'interpreter','tex');
view([-55,14])


end

function v = f(x,y)
v = 2*pi^2*sin(pi*x).*sin(pi*y);
end

function v = g_D(x,y)
v = 0*x.*y;
end