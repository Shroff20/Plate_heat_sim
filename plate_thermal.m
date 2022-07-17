clear all
close all
clc
% 
% syms T(t,x,y)
% syms eps sig tz hc Ta rho Cp k t v
% syms hCoeff(t)
% 
% 
% % Boundary conditions and loads
% v = 1+heaviside(t-1000)*100;      % fluid velocity (m/s)  (must be 2 to 20)
% U1 = @(location, state) 1000*sin(.001*state.time);
% U2 = @(location, state) 1000*sin(.001*state.time);
% U3 = @(location, state) 1000*sin(.001*state.time);
% U4 = @(location, state) 1000*sin(.001*state.time);
% Tinitial = 2000;
% 
% % Material properties
% kThermal = 400;             % thermal conductivity of copper, W/(m-K)
% rhoCopper = 8960;           % density of copper, kg/m^3
% specificHeat = 386;         % specific heat of copper, J/(kg-K)
% thick = 0.01;               % plate thickness in meters
% stefanBoltz = 5.670373e-8;  % Stefan-Boltzmann constant, W/(m^2-K^4)
% hCoeff = 1+v;               % convection coefficient, W/(m^2-K)
% tAmbient = 300;             % the ambient temperature
% emiss = 0.5;                % emissivity of the plate surface
% 
% % Geometry and mesh
% width = 1; % (m)
% height = 1; % (m)
% elSize = 0.1; % element size
% 
% % solver
% tend = 10000;  % end time (s)
% dt = 50; % time step (s)

setup = GetDefaultInputs()

model = makeModel(setup)

hfig = PlotGeoemtry(model);

%%

U = {setup.loads.U1, setup.loads.U2, setup.loads.U3, setup.loads.U4};
Tinit = setup.ICs.Tinit 
tend = setup.solver.tend;          
dt = setup.solver.dt;           
RelativeTolerance = setup.solver.RelativeTolerance ;
AbsoluteTolerance = setup.solver.SolverOptions.AbsoluteTolerance ;


applyBoundaryCondition(model,'dirichlet','edge',1,'u',U{1});
applyBoundaryCondition(model,'dirichlet','edge',2,'u',U{2});
applyBoundaryCondition(model,'dirichlet','edge',3,'u',U{3});
applyBoundaryCondition(model,'dirichlet','edge',4,'u',U{4});

setInitialConditions(model,Tinit);
tlist = 0:dt:tend;

model.SolverOptions.RelativeTolerance = RelativeTolerance; 
model.SolverOptions.AbsoluteTolerance =  AbsoluteTolerance;

R = solvepde(model,tlist);
u = R.NodalSolution;
figure; 
plot(tlist,u(:,:));
grid on
title 'Temperature Along the Top Edge of the Plate as a Function of Time'
xlabel 'Time (s)'
ylabel 'Temperature (K)'


%%
figure;
pdeplot(model,'XYData',u(:,end),'Contour','on','ColorMap','jet');
title(sprintf('Plate Temperature\nt = %d seconds', ...
  tlist(1,end)));
xlabel 'x'
ylabel 'y'
axis equal;


function setup = GetDefaultInputs()

setup = [];

syms t

% Boundary conditions and loads
setup.loads.v = 1+heaviside(t-1000)*100;      % fluid velocity (m/s)  (must be 2 to 20)
setup.loads.U1 = @(location, state) 1000*sin(.001*state.time);
setup.loads.U2 = @(location, state) 1000*sin(.001*state.time);
setup.loads.U3 = @(location, state) 1000*sin(.001*state.time);
setup.loads.U4 = @(location, state) 1000*sin(.001*state.time);
setup.ICs.Tinit = 2000;

% Material properties
setup.material.k = 400;             % thermal conductivity of copper, W/(m-K)
setup.material.rho = 8960;           % density of copper, kg/m^3
setup.material.cp = 386;         % specific heat, J/(kg-K)
setup.material.tz = 0.01;               % plate thickness in meters
setup.material.sig = 5.670373e-8;  % Stefan-Boltzmann constant, W/(m^2-K^4)
setup.material.h = 1+setup.loads.v;               % convection coefficient, W/(m^2-K)
setup.material.Ta = 300;             % the ambient temperature
setup.material.emiss = 0.5;                % emissivity of the plate surface

% Geometry and mesh
setup.geometry.width = 1;                  %  plate width, m
setup.geometry.height = 1;                 %  plate height, m
setup.geometry.elSize = 0.1;               %  element size, m

% solver
setup.solver.tend = 10000;                 % end time (s)
setup.solver.dt = 50;                      % time step (s)
setup.solver.RelativeTolerance = 1.0e-3;   % solver relative tolerance
setup.solver.SolverOptions.AbsoluteTolerance = 1.0e-3;  % solver absolute tolerance

end 

function model = makeModel(setup)

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

