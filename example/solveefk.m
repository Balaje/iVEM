function solveefk

clc
clear
format short
close all

mesh = load('squares_linear10');
t0 = 0;  tf = 0.1;  delt = 0.01;

[usol,vsol] = efk(mesh, @f, @g_Du, @g_Dv, @u0, tf, t0, delt);
exact = @(x,y) 10*tf*sin(pi*x)*sin(pi*y);

[l2err,h1err] = l2error(mesh,exact,usol,1);

figure(1)
subplot(1,2,1);
plot_solution(mesh,usol);
grid on
str = 'Solution of the BBM equation';
title(str,'FontSize',14,'interpreter','tex');
view([-112,13])
l2 = num2str(double(l2err));
h1 = num2str(double(h1err));
txt = ['||e||_{L^2} = ',l2,'   |e|_{H^1} = ',h1];
text(1,1,max(usol)+0.05,txt,'FontSize',11)
subplot(1,2,2);
plot_solution(mesh,zeros(length(usol),1));
hold on
scatter(mesh.vertices(:,1),mesh.vertices(:,2),5,'rs');
title('Virtual Element Mesh','FontSize',12);
pause(0.01);
hold off

end

function v = u0(x,y)
    v = 0*x.*y;
end

function v = f(x,y,t)    
    v = (2*sin(pi*x)*sin(pi*y)*(50*t*pi^2 - 25*t + t*pi^4 + 2500*t^3*sin(pi*x)^2*sin(pi*y)^2 + 25))/5;
end

function v = g_Du(x,y,t)
    v = 0.*x.*y*t;
end

function v = g_Dv(x,y,t)
    v = 0.*x.*y*t;
end