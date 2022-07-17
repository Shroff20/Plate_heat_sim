clear all
close all
clc

syms T(t,x,y)
syms eps sig tz hc Ta rho Cp k t v
syms hCoeff(t)


% Boundary conditions and loads
v = 1+heaviside(t-1000)*100;      % fluid velocity (m/s)  (must be 2 to 20)
U1 = @(location, state) 1000*sin(.001*state.time);
U2 = @(location, state) 1000*sin(.001*state.time);
U3 = @(location, state) 1000*sin(.001*state.time);
U4 = @(location, state) 1000*sin(.001*state.time);
Tinitial = 2000;

% Material properties
kThermal = 400;             % thermal conductivity of copper, W/(m-K)
rhoCopper = 8960;           % density of copper, kg/m^3
specificHeat = 386;         % specific heat of copper, J/(kg-K)
thick = 0.01;               % plate thickness in meters
stefanBoltz = 5.670373e-8;  % Stefan-Boltzmann constant, W/(m^2-K^4)
hCoeff = 1+v;               % convection coefficient, W/(m^2-K)
tAmbient = 300;             % the ambient temperature
emiss = 0.5;                % emissivity of the plate surface

% Geometry and mesh
width = 1; % (m)
height = 1; % (m)
elSize = 0.1; % element size

% solver
tend = 10000;  % end time (s)
dt = 50; % time step (s)


%%  Make 

Qc = hc*(T - Ta);  % heat transfer due to convection
Qr = eps*sig*(T^4 - Ta^4); % heat transfer due to radiation
pdeeq = (rho*Cp*tz*diff(T,t) - k*tz*laplacian(T,[x,y]) + 2*Qc + 2*Qr) % pde for temperature of a thin plate

symCoeffs = pdeCoefficients(pdeeq,T,'Symbolic',true)

symVars = [eps sig tz hc Ta rho Cp k];
symVals = [emiss stefanBoltz thick hCoeff tAmbient rhoCopper specificHeat kThermal];
symCoeffs = subs(symCoeffs,symVars,symVals);
coeffs = pdeCoefficientsToDouble(symCoeffs)
numberOfPDE = 1;
model = createpde(numberOfPDE);


%%
gdm = [3 4 0 width width 0 0 0 height height]';
g = decsg(gdm,'S1',('S1')');

geometryFromEdges(model,g);

figure; 
pdegplot(model,'EdgeLabels','on'); 
axis([-0.1 1.1 -0.1 1.1]);
title('Geometry with Edge Labels Displayed');


msh = generateMesh(model,'Hmax',elSize);
figure; 
pdeplot(model); 
axis equal
title('Plate with Triangular Element Mesh');
xlabel('X-coordinate, meters');
ylabel('Y-coordinate, meters');

%%
specifyCoefficients(model,'m',coeffs.m,'d',coeffs.d, ...
    'c',coeffs.c,'a',coeffs.a,'f',coeffs.f);


U = {U1, U2, U3, U4};
applyBoundaryCondition(model,'dirichlet','edge',1,'u',U{1});
applyBoundaryCondition(model,'dirichlet','edge',2,'u',U{2});
applyBoundaryCondition(model,'dirichlet','edge',3,'u',U{3});
applyBoundaryCondition(model,'dirichlet','edge',4,'u',U{4});

setInitialConditions(model,Tinitial);
%setInitialConditions(model,1000,'edge',1);

tlist = 0:dt:tend;

model.SolverOptions.RelativeTolerance = 1.0e-3; 
model.SolverOptions.AbsoluteTolerance = 1.0e-4;

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

syms T(t,x,y)
syms eps sig tz hc Ta rho Cp k t v
syms hCoeff(t)

% Boundary conditions and loads
setup.loads.v = 1+heaviside(t-1000)*100;      % fluid velocity (m/s)  (must be 2 to 20)
setup.loads.U1 = @(location, state) 1000*sin(.001*state.time);
setup.loads.U2 = @(location, state) 1000*sin(.001*state.time);
setup.loads.U3 = @(location, state) 1000*sin(.001*state.time);
setup.loads.U4 = @(location, state) 1000*sin(.001*state.time);
setup.initialConditions.Tinitial = 2000;

% Material properties
setup.material.kThermal = 400;             % thermal conductivity of copper, W/(m-K)
setup.material.rhoCopper = 8960;           % density of copper, kg/m^3
setup.material.specificHeat = 386;         % specific heat of copper, J/(kg-K)
setup.material.thick = 0.01;               % plate thickness in meters
setup.material.stefanBoltz = 5.670373e-8;  % Stefan-Boltzmann constant, W/(m^2-K^4)
setup.material.hCoeff = 1+v;               % convection coefficient, W/(m^2-K)
setup.material.tAmbient = 300;             % the ambient temperature
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


