clear all
close all
clc


setup = GetDefaultInputs()
model = MakeModel(setup)
model = ApplyLoadsAndICs(model, setup);
[model, results] = SolvePDE(model, setup);

Tinterp = InterpolateResults(results, setup);

hfig1 = PlotGeoemtry(model);
hfig2 = PlotNodalSolutions(results);
hfig3 = PlotFieldSolution(model, results, -1)  % -1 for last frame
hfig4 = PlotInterpolatedSolution(Tinterp, setup, -1);


%%

function hfig = PlotInterpolatedSolution(Tinterp, setup, frame)



xinterp = setup.geometry.xinterp;
yinterp = setup.geometry.yinterp;

[X, Y] = meshgrid(xinterp, yinterp);

if frame ==-1
    frame  = size(Tinterp, 2);
end 

Tinterp = Tinterp(:, frame);
Tinterp = reshape(Tinterp, size(X, 1), []);

hfig = figure();
colormap('jet')
contourf(X, Y, Tinterp, 250, 'LineColor','none')
colorbar()
axis equal

end 



function Tinterp = InterpolateResults(results, setup)

[X,Y] = meshgrid(setup.geometry.xinterp, setup.geometry.yinterp);
interpPoints = [X(:),Y(:)]';
Nt = 1:length(results.SolutionTimes);
Tinterp = interpolateSolution(results,interpPoints, Nt);

end 


function setup = GetDefaultInputs()

setup = [];

syms t

% Boundary conditions and loads
setup.loads.v = 1+heaviside(t-1000)*100;      % fluid velocity (m/s)  (must be 2 to 20)
setup.loads.U1 = @(location, state) 1000*sin(.001*state.time) + 1000;
setup.loads.U2 = @(location, state) 1000*sin(.001*state.time) + 1000;
setup.loads.U3 = @(location, state) 1000*sin(.001*state.time) + 1000;
setup.loads.U4 = @(location, state) 1000*sin(.001*state.time) + 1200;
setup.ICs.Tinit = 2000;

% Material properties
setup.material.k = 400;             % thermal conductivity of copper, W/(m-K)
setup.material.rho = 8960;           % density of copper, kg/m^3
setup.material.cp = 386;         % specific heat, J/(kg-K)
setup.material.tz = 0.01;               % plate thickness in meters
setup.material.sig = 5.670373e-8;  % Stefan-Boltzmann constant, W/(m^2-K^4)
setup.material.h = 1+setup.loads.v;               % convection coefficient, W/(m^2-K)
setup.material.Ta = 100;             % the ambient temperature
setup.material.emiss = 0.5;                % emissivity of the plate surface

% Geometry and mesh
setup.geometry.width = 1;                  %  plate width, m
setup.geometry.height = 1;                 %  plate height, m
setup.geometry.elSize = 0.05;               %  element size, m
setup.geometry.xinterp = 0:setup.geometry.elSize:setup.geometry.width;               %  solution interpolation points
setup.geometry.yinterp = 0:setup.geometry.elSize:setup.geometry.height;              %  solution interpolation points

% solver
setup.solver.tend = 10000;                 % end time (s)
setup.solver.dt = 50;                      % time step (s)
setup.solver.RelativeTolerance = 1.0e-3;   % solver relative tolerance
setup.solver.SolverOptions.AbsoluteTolerance = 1.0e-3;  % solver absolute tolerance

end 

function model = MakeModel(setup)

syms T(t,x,y)

emiss = setup.material.emiss;
sig = setup.material.sig;
tz = setup.material.tz;
hc = setup.material.h ;
Ta = setup.material.Ta;
rho = setup.material.rho;
cp = setup.material.cp ;
k = setup.material.k;
width = setup.geometry.width;
height = setup.geometry.height;
elSize = setup.geometry.elSize;

Qc = hc*(T - Ta);  % heat transfer due to convection
Qr = emiss*sig*(T^4 - Ta^4); % heat transfer due to radiation
pdeeq = (rho*cp*tz*diff(T,t) - k*tz*laplacian(T,[x,y]) + 2*Qc + 2*Qr); % pde for temperature of a thin plate

symCoeffs = pdeCoefficients(pdeeq,T,'Symbolic',true);
coeffs = pdeCoefficientsToDouble(symCoeffs);
numberOfPDE = 1;
model = createpde(numberOfPDE);

% make mesh
gdm = [3 4 0 width width 0 0 0 height height]';
g = decsg(gdm,'S1',('S1')');

geometryFromEdges(model,g);
msh = generateMesh(model,'Hmax',elSize);
specifyCoefficients(model,'m',coeffs.m,'d',coeffs.d, ...
    'c',coeffs.c,'a',coeffs.a,'f',coeffs.f);

end 


function hfig = PlotGeoemtry(model)

hfig = figure();
subplot(1,2,1);
pdegplot(model,'EdgeLabels','on'); 
axis([-0.1 1.1 -0.1 1.1]);
title('geometry');
 
subplot(1,2,2);
pdeplot(model); 
axis equal
axis([-0.1 1.1 -0.1 1.1]);
title('mesh');
xlabel('x');
ylabel('y');

end 

function model = ApplyLoadsAndICs(model, setup)


U = {setup.loads.U1, setup.loads.U2, setup.loads.U3, setup.loads.U4};
Tinit = setup.ICs.Tinit ;
% tend = setup.solver.tend;          
% dt = setup.solver.dt;    
% RelativeTolerance = setup.solver.RelativeTolerance ;
% AbsoluteTolerance = setup.solver.SolverOptions.AbsoluteTolerance ;


applyBoundaryCondition(model,'dirichlet','edge',1,'u',U{1});
applyBoundaryCondition(model,'dirichlet','edge',2,'u',U{2});
applyBoundaryCondition(model,'dirichlet','edge',3,'u',U{3});
applyBoundaryCondition(model,'dirichlet','edge',4,'u',U{4});

setInitialConditions(model,Tinit);

end 

function [model, results] = SolvePDE(model, setup)

tend = setup.solver.tend;          
dt = setup.solver.dt;    
RelativeTolerance = setup.solver.RelativeTolerance ;
AbsoluteTolerance = setup.solver.SolverOptions.AbsoluteTolerance ;

tspan = 0:dt:tend;
model.SolverOptions.RelativeTolerance = RelativeTolerance; 
model.SolverOptions.AbsoluteTolerance =  AbsoluteTolerance;

results = solvepde(model,tspan);

end 


function hfig = PlotNodalSolutions(results)
u = results.NodalSolution;
tlist = results.SolutionTimes;
hfig = figure(); 
plot(tlist,u(:,:));
grid on
title 'nodal temperatures';
xlabel 'time (s)';
ylabel 'temperature';
end 

function hfig = PlotFieldSolution(model, results, frame)

tlist = results.SolutionTimes;
u = results.NodalSolution;

if frame == -1
    frame = length(tlist);
end 

%%
hfig = figure;
pdeplot(model,'XYData',u(:,frame),'Contour','on','ColorMap','jet');
title(sprintf('Plate Temperature\nt = %d seconds', ...
  tlist(1,frame)));
xlabel 'x'
ylabel 'y'
axis equal;

end 


