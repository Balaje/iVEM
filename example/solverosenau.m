function solverosenau

clc
clear
format short
close all

meshes = {'squares_linear10','squares_linear15',...
    'squares_linear20','squares_linear25'};
%meshes = {'non-convex'};

for i=1:length(meshes)
    mesh = load(meshes{i});
    t0 = 0;  tf = 0.1;  delt = 0.01;
    
    [usol,vsol] = rosenau(mesh, @f, @g_D, @u0, tf, t0, delt);
    exact = @(x,y) 2*pi^2*10*tf*sin(pi*x)*sin(pi*y);
    
    [l2err,h1err] = l2error(mesh,exact,vsol,1,false);
    
    figure(i)
    subplot(1,2,1);
    plot_solution(mesh,vsol);
    grid on
    str = 'Solution of the ROSENAU equation';
    title(str,'FontSize',14,'interpreter','tex');
    view([-112,13])
    l2 = num2str(double(l2err));
    h1 = num2str(double(h1err));
    txt = ['||e||_{L^2} = ',l2,'   |e|_{H^1} = ',h1];
    text(1,1,max(vsol)+0.05,txt,'FontSize',11)
    subplot(1,2,2);
    plot_solution(mesh,zeros(length(usol),1));
    hold on
    scatter(mesh.vertices(:,1),mesh.vertices(:,2),5,'rs');
    title('Virtual Element Mesh','FontSize',12);
    pause(0.01);
    hold off
    
end

end

function v = u0(x,y)
v = 0*x.*y;
end

function v = f(x,y,t)
v = 10.*sin(pi*x).*sin(pi*y).*(- 20*pi*sin(pi*(x + y))*t^2 + 4*pi^4 + 1);
end

function v = g_D(x,y,t)
v = 0.*x.*y*t;
end