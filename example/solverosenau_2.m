function solverosenau_2

clc
clear
format short
close all


% plotmeshes = {'squares_linear10','squares_linear15',...
%     'squares_linear20','squares_linear25'};
% meshes = {'squares_quadratic10','squares_quadratic15',...
%     'squares_quadratic20','squares_quadratic25'};

plotmeshes = {'voronoi100'};
meshes = {'voronoi_quadratic100'};

for i=1:length(meshes)
    mesh = load(meshes{i});
    plotmesh = load(plotmeshes{i});
    t0 = 0;  tf = 0.1;  delt = 0.01;
    
    [usol,vsol] = rosenau_2(mesh, @f, @g_D, @u0, tf, t0, delt, plotmesh);
    exact = @(x,y) 10*tf*sin(pi*x)*sin(pi*y);
    
    [l2err,h1err] = l2error(mesh,exact,usol,2,false);
    
    noofvertices = length(plotmesh.vertices);
    u_verts = usol(1:noofvertices);
    
    figure(i)
    subplot(1,2,1);
    plot_solution(plotmesh,u_verts);
    grid on
    str = 'Solution of BBM equation';
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
    pause(0.01);
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