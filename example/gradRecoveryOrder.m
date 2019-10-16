%% Program to compute the rate of convergence for the Gradient Recovery Method

clear
close all

set(0,'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');
set(0,'defaultaxesfontsize',20);
set(0,'defaultlegendfontsize',15);

%% Cell array for the meshes

vem = 2; % Order of VEM 
isInt = false;

if(vem==1)
    meshes = {'squares_linear10.mat','squares_linear15.mat','squares_linear20.mat'...
        'squares_linear25.mat','squares_linear30.mat'};
else    
    meshes = {'squares_quadratic10.mat','squares_quadratic15.mat',...
        'squares_quadratic20.mat','squares_quadratic25.mat','squares_quadratic30.mat'};
end

%meshParam = zeros(length(meshes),1);
meshParam = [1/10,1/15,1/20,1/25,1/30]';
soll2err = zeros(length(meshes),1);
solh1err = zeros(length(meshes),1);
gxl2err = zeros(length(meshes),1);
gyl2err = zeros(length(meshes),1);

f = @(x,y) 2*pi^2.*sin(pi*x).*sin(pi*y);
g_D = @(x,y) 0.*x.*y;

for m=1:length(meshes)
    % Get the error curves for the solution.
    mesh = load(meshes{m});
    [~,~,dia] = geo(mesh); % Obtain the diameter of the all the polygons in the mesh.
    if(vem==2)        
        usol = poisson_2(mesh, f, g_D);
    else
        usol = poisson(mesh, f, g_D);
    end
    exact = @(x,y)sin(pi*x)*sin(pi*y);
    [soll2err(m),solh1err(m)] = l2error(mesh,exact,usol,vem,false);    
     
    % Gradient Recovery for the solution and the l2-error
    [Mx, My, M] = gradrecovery_matrix(mesh,vem);    
    D = diag(1./sum(M,1));
    
    uhx = D*Mx*usol;
    uhy = D*My*usol;
    [gxl2err(m), ~] = l2error(mesh, @(x,y) pi*cos(pi*x)*sin(pi*y), uhx, vem, isInt);
    [gyl2err(m), ~] = l2error(mesh, @(x,y) pi*sin(pi*x)*cos(pi*y), uhy, vem, isInt);
    
     gL2err = sqrt(gxl2err.^2 + gyl2err.^2);    
end

%% Plot the error curves
fig=figure(1);
loglog(meshParam,gL2err,'r+-','linewidth',1.5);
hold on
loglog(meshParam,solh1err,'mo-','linewidth',1.5);
hold on
loglog(meshParam,2.1*meshParam.^(1.5),'k','linewidth',1);
grid on
legend('$||\nabla u - \nabla u_h ||_{0,h}$',...
    '$|u-u_h|_{1,h}$','$O(h^{1.5})$','Location','best');
legend('boxoff');

%% Print the rates of convergence.
fprintf('The L2 rate of the recovered gradient..\n');
disp(log(gL2err(1:end-1)./gL2err(2:end))./(log(meshParam(1:end-1)./meshParam(2:end))));
fprintf('The L2 rate of the solution..\n');
disp(log(soll2err(1:end-1)./soll2err(2:end))./(log(meshParam(1:end-1)./meshParam(2:end))));
fprintf('The H1 rate of the solution..\n');
disp(log(solh1err(1:end-1)./solh1err(2:end))./(log(meshParam(1:end-1)./meshParam(2:end))));