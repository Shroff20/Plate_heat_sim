clear all
close all
clc

%% INPUTS

U = cell(5,1);
U{1} = @(location, state) 10*sin(.001*state.time) + 800;                             % temperature of plate edge 1
U{2} = @(location, state) 50*sin(.001*state.time) + 1000;                            % temperature of plate edge 2
U{3} = @(location, state) 100*sin(.001*state.time) + 1200;                           % temperature of plate edge 3
U{4} = @(location, state) 20*sin(.001*state.time) + heaviside(state.time-2000)*1400; % temperature of plate edge 4
U{5} = @(location, state) 1+heaviside(state.time-1000)*20;                           % velocity of fluid for convection calcs
Tinit = 70;                                                                          % initial temperature for plate


%% RUN AND PLOT

[dataTable, detailedSimData] =  RunPlateSim(U, Tinit);

%%

inputs = detailedSimData.dataTable.Inputs;
T = detailedSimData.dataTable.Temperature;
S = detailedSimData.dataTable.Stress;
setup = detailedSimData.setup;
thermalModel = detailedSimData.thermalModel;
structuralModel = detailedSimData.structuralModel;
thermalResults = detailedSimData.thermalResults;
structuralResults = detailedSimData.structuralResults;


hfig = figure();

subplot(4,2,1);
pdegplot(thermalModel,'EdgeLabels','on', 'FaceLabels', 'on', 'FaceAlpha',0.5'); 
axis equal
title('geometry');
 
subplot(4,2,2);
pdeplot3D(thermalModel); 
axis equal
title('mesh');
xlabel('x');
ylabel('y');

subplot(4,2,3)
plot(inputs(:,1:4), 'linewidth', 2)
xlabel('time index')
ylabel('temperature (°F)')
title('Inputs')
yyaxis right
plot(inputs(:,5), 'linewidth', 2)
legend({'T_1', 'T_2', 'T_3', 'T_4', 'velocity'}, 'Location', 'northeastoutside')


subplot(4,2,4)
plot(inputs(:,5), 'linewidth', 2)
legend('fluid velocity')
xlabel('time index')
ylabel('velocity (m/s)')
sgtitle('Results')
title('Inputs')

subplot(4,2,5)
plot(T, 'linewidth', 1)
ylabel('temperature (°F)')
title('Nodal Tempereratures')
xlabel('time index')

subplot(4,2,6)
plot(S, 'linewidth', 1)
ylabel('stress (psi)')
title('Nodal Stresses')
xlabel('time index')


subplot(4,2,7)
xinterp = setup.geometry.xinterp;
yinterp = setup.geometry.yinterp;
[X, Y] = meshgrid(xinterp, yinterp);
Trs = reshape( T(end, :), size(X, 1), []);
colormap('jet')
contourf(X, Y, Trs, 250, 'LineColor','none')
colorbar()
axis equal
title('Final Temperatures')

subplot(4,2,8)
Trs = reshape( S(end, :), size(X, 1), []);
colormap('jet')
contourf(X, Y, Trs, 250, 'LineColor','none')
colorbar()
axis equal
title('Final Von Mises Stresses')




% subplot(4,2,7)
% pdeplot3D(thermalModel,'ColorMapData',thermalResults.Temperature(:, end));
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal;
% title('Final Temperature')
% 
% subplot(4,2,8)
% pdeplot3D(structuralModel,'ColorMapData',structuralResults{end}.VonMisesStress);
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal;
% title('Final Stress')


% hfig1 = PlotGeoemtry(model);
% hfig2 = PlotNodalSolutions(results);
% hfig3 = PlotFieldSolution(model, results, -1);  % -1 for last frame
% hfig4 = PlotInterpolatedSolution(Tinterp, setup, -1);
% hfig5 = PlotInputs(results, inputs);
%%


%%
% figure()
% plot(S')
% 
% 
% figure()
% pdeplot3D(structuralModel, ...
%           "ColorMapData",structuralResults{end}.VonMisesStress, ...
%           "Deformation",structuralResults{end}.Displacement)
%       
% 
% disp('     * done interpolating results')
% PlotInterpolatedSolution(S, setup, -1)      
% PlotInterpolatedSolution(Tinterp, setup, -1)      

%% FUNCTIONS

function [dataTable, detailedSimData] =  RunPlateSim(U, Tinit)
disp('-STARTING--------------------------')
setup = GetDefaultInputs();
setup = ModifyLoadsAndICs(setup, U, Tinit);
thermalModel = MakeModel(setup);
[thermalModel, thermalResults, inputs] = SolveTemperature(thermalModel, setup);
T = InterpolateResults(thermalResults, setup);
[S, structuralResults, structuralModel] = SolveStresses(thermalModel, thermalResults, setup);


dataTable = table();
dataTable.time = thermalResults.SolutionTimes';
dataTable.Temperature = T';
dataTable.Stress = S';
dataTable.Inputs = inputs';
disp(dataTable);

detailedSimData = [];
detailedSimData.setup = setup;
detailedSimData.thermalModel = thermalModel;
detailedSimData.thermalResults = thermalResults;
detailedSimData.structuralModel = structuralModel;
detailedSimData.structuralResults = structuralResults;
detailedSimData.dataTable = dataTable;


disp('-SIMULATION COMPLETE---------------')

end 


function [S, structuralResults, structuralModel] = SolveStresses(thermalModel, thermalResults, setup)

E = setup.material.E;
poisson = setup.material.poisson;
alpha = setup.material.alpha;
density = setup.material.rho;

structuralModel = createpde("structural","static-solid");
structuralModel.Mesh = thermalModel.Mesh;
structuralModel.Geometry = thermalModel.Geometry;

pdegplot(structuralModel,"EdgeLabels","on", "VertexLabels","on", 'FaceLabels', 'on');
structuralProperties(structuralModel,"YoungsModulus",E, "PoissonsRatio",poisson, "MassDensity",density, 'CTE',alpha);
%structuralBC(structuralModel,"Vertex",[1:4],"Constraint","fixed");
structuralBC(structuralModel,"Face",[1:6],'Constraint','roller');    
[X,Y,Z] = meshgrid(setup.geometry.xinterp, setup.geometry.yinterp, setup.geometry.zinterp);
interpPoints = [X(:),Y(:), Z(:)]';

Nt = length(thermalResults.SolutionTimes);
S = nan(size(interpPoints,2), Nt);
structuralResults = cell(1, Nt);
for i = 1:Nt
    structuralBodyLoad(structuralModel,'Temperature',thermalResults,'TimeStep', i);
    structuralResults{i} = solve(structuralModel);  
    S(:,i) = interpolateVonMisesStress(structuralResults{i},interpPoints);
end 

disp('     * done solving structural model')

end 


function hfig = PlotInputs(results, inputs)
hfig = figure();
tspan  = results.SolutionTimes;
subplot(2,1,1)
plot(tspan, inputs(1:4,:), 'linewidth', 2)
legend({'T_1', 'T_2', 'T_3', 'T_4'})
xlabel('t (s)')
ylabel('temperature (°F)')
subplot(2,1,2)
plot(tspan, inputs(5,:), 'linewidth', 2)
legend('fluid velocity')
xlabel('t (s)')
ylabel('velocity (m/s)')
sgtitle('Inputs')
end 



function setup = ModifyLoadsAndICs(setup, U, Tinit)

setup.ICs.Tinit = Tinit;
setup.loads.v = U{5};
setup.loads.BCs = U(1:4);
disp('     * done modifying BCs and loads')

end 

function hfig = PlotInterpolatedSolution(Tinterp, setup, frame)

xinterp = setup.geometry.xinterp;
yinterp = setup.geometry.yinterp;

[X, Y] = meshgrid(xinterp, yinterp);

if frame ==-1
    frame  = size(Tinterp, 2);
end 

Tinterp2 = reshape( Tinterp(:, frame), size(X, 1), []);

hfig = figure();
colormap('jet')
contourf(X, Y, Tinterp2, 250, 'LineColor','none')
colorbar()
axis equal

end 



function Tinterp = InterpolateResults(results, setup)

[X,Y,Z] = meshgrid(setup.geometry.xinterp, setup.geometry.yinterp, setup.geometry.zinterp);
interpPoints = [X(:),Y(:), Z(:)]';
Nt = 1:length(results.SolutionTimes);
Tinterp = interpolateTemperature(results,interpPoints, Nt);
disp('     * done interpolating results')

end 


function setup = GetDefaultInputs()

setup = [];

syms t

% Boundary conditions and loads
setup.loads.v = @(location, state) 1+heaviside(state.time-1000)*200;      % fluid velocity (m/s)  (must be 2 to 20)
U1 = @(location, state) 1000*sin(.001*state.time) + 1000;
U2 = @(location, state) 1000*sin(.001*state.time) + 1000;
U3 = @(location, state) 1000*sin(.001*state.time) + 1000;
U4 = @(location, state) 1000*sin(.001*state.time) + 1200;
setup.loads.BCs = {U1, U2, U3, U4};
setup.ICs.Tinit = 2000;

% Material properties
setup.material.k = 400;             % thermal conductivity, W/(m-K)
setup.material.rho = 8960;           % density , kg/m^3
setup.material.cp = 386;         % specific heat, J/(kg-K)
setup.material.sig = 5.670373e-8;  % Stefan-Boltzmann constant, W/(m^2-K^4)
setup.material.h = @(v) 1+v.^2;               % convection coefficient, W/(m^2-K) as a function of fluid velocity
setup.material.Ta = 70;             % the ambient temperature
setup.material.emiss = 0.5;                % emissivity of the plate surface
setup.material.E = 2E9;          % youngs modulus
setup.material.poisson = .3;       % Poisson ratio
setup.material.alpha = 2.4e-5;       % coeff of thermal expansion
setup.material.Tref = 70;             % reference temperature for thermal expansion

% Geometry and mesh
setup.geometry.width = 1;                  %  plate width, m
setup.geometry.height = 1;                 %  plate height, m
setup.geometry.thickness = .1;                 %  plate thickness, m
setup.geometry.elSize = 0.1;               %  element size, m
setup.geometry.xinterp = 0:setup.geometry.elSize:setup.geometry.width;               %  solution interpolation points
setup.geometry.yinterp = 0:setup.geometry.elSize:setup.geometry.height;              %  solution interpolation points
setup.geometry.zinterp = setup.geometry.thickness/2;                                 %  solution interpolation points

% solver
setup.solver.tend = 3000;                 % end time (s)
setup.solver.dt = 100;                      % time step (s)
setup.solver.RelativeTolerance = 1.0e-2;   % solver relative tolerance
setup.solver.SolverOptions.AbsoluteTolerance = 1.0e-2;  % solver absolute tolerance
disp('     * done getting default inputs')

end 


function thermalModel = MakeModel(setup)

rho = setup.material.rho;
cp = setup.material.cp ;
k = setup.material.k;
width = setup.geometry.width;
height = setup.geometry.height;
thickness = setup.geometry.thickness;
elSize = setup.geometry.elSize;

% Qc = hc*(T - Ta);  % heat transfer due to convection
% Qr = emiss*sig*(T^4 - Ta^4); % heat transfer due to radiation
% pdeeq = (rho*cp*tz*diff(T,t) - k*tz*laplacian(T,[x,y]) + 2*Qc + 2*Qr); % pde for temperature of a thin plate
% symCoeffs = pdeCoefficients(pdeeq,T,'Symbolic',true);

thermalModel = createpde('thermal','transient');

% make mesh
[x,y,z] = meshgrid([0 width], [0 height], [0 thickness]);
x = x(:);
y = y(:);
z = z(:);
K = convhull(x,y,z);
nodes = [x';y';z'];
elements = K';
geometryFromMesh(thermalModel,nodes,elements);
generateMesh(thermalModel, 'Hmax',elSize)

thermalProperties(thermalModel,'ThermalConductivity',k,...
                               'MassDensity',rho,...
                               'SpecificHeat',cp);
end 


function hfig = PlotGeoemtry(model)

hfig = figure();
subplot(1,2,1);
pdegplot(model,'EdgeLabels','on', 'FaceLabels', 'on', 'FaceAlpha',0.5'); 
axis equal
title('geometry');
 
subplot(1,2,2);
pdeplot3D(model); 
axis equal
title('mesh');
xlabel('x');
ylabel('y');

end 


function [model, results, inputs] = SolveTemperature(model, setup)
Ta = setup.material.Ta;
tend = setup.solver.tend;          
dt = setup.solver.dt;    
RelativeTolerance = setup.solver.RelativeTolerance ;
AbsoluteTolerance = setup.solver.SolverOptions.AbsoluteTolerance ;
U = setup.loads.BCs;
Tinit = setup.ICs.Tinit ;
%hc = setup.loads.v;%@(loc, state) setup.material.h(setup.loads.v(loc, state)) ;
hc = @(loc, state) setup.material.h(setup.loads.v(loc, state)) ;

thermalBC(model,'Face',[2] ,'Temperature', U{1});
thermalBC(model,'Face',[3],'Temperature', U{2});
thermalBC(model,'Face',[4],'Temperature', U{3});
thermalBC(model,'Face',[5],'Temperature', U{4});
thermalBC(model,'Face',[1],'ConvectionCoefficient', hc, 'AmbientTemperature',Ta);
thermalBC(model,'Face',[6],'ConvectionCoefficient', hc, 'AmbientTemperature',Ta);

thermalIC(model,Tinit);

tspan = 0:dt:tend;
model.SolverOptions.RelativeTolerance = RelativeTolerance; 
model.SolverOptions.AbsoluteTolerance =  AbsoluteTolerance;

results = solve(model,tspan);

% calculate and record inputs at each timestep
inputs = nan(5, length(tspan));
for i =1:4
    tmp_fuc = setup.loads.BCs{i};
    tmp_state = [];
    tmp_state.time = tspan;
    inputs(i,:) = tmp_fuc([], tmp_state);
end 
inputs(5,:) = setup.loads.v([], tmp_state);


disp('     * done solving pde')

end 


function hfig = PlotNodalSolutions(results)
u = results.Temperature;
tlist = results.SolutionTimes;
hfig = figure(); 
plot(tlist,u(:,:));
grid on
title 'nodal temperatures';
xlabel 'time (s)';
ylabel 'temperature';
end 


function hfig = PlotFieldSolution(model, results, frame)

%tlist = results.SolutionTimes;
 
u = results.Temperature;
if frame == -1
    frame = size(u, 2);
end
u = u(:, frame);

hfig = figure;
pdeplot3D(model,'ColorMapData',u);
xlabel('x')
ylabel('y')
zlabel('z')
axis equal;

end 

