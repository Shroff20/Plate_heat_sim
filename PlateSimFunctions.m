classdef PlateSimFunctions
    %PLATESIMFUNCTIONS Summary of this class goes here
    %   Detailed explanation goes here
    
    
    methods(Static)
        
        function RunSimPipeline(inputs, casename, output_dir)
        
        [dataTable, detailedSimData] =  PlateSimFunctions.RunPlateSim(inputs);
        hfig  = PlateSimFunctions.MakePlotsFromDetailedSimData(detailedSimData);
        PlateSimFunctions.SaveDataAndPlots(char(casename), output_dir, dataTable, hfig)

        end 
        
        
        function SaveDataAndPlots(casename, output_dir, dataTable, hfig)
        mkdir([output_dir, filesep, casename])
        csv_name = [output_dir, filesep,casename , filesep , 'data.csv'];
        pdf_name = [output_dir, filesep, casename , filesep , 'plots.pdf'];
        writetable(dataTable, csv_name)
        set(hfig,'PaperSize',[20 20]); 
        print(hfig,pdf_name,'-dpdf', '-painters', '-fillpage') 
        disp(['     * Saved output to ', [output_dir, filesep, casename]])

        end 


        function inputs = GetExampleInput()
        % cols are T1, T2, T3, T4, fluid velocity
        tspan = 0:1:600;  % time is in units of dt seconds per point in tspan
        inputs = nan(length(tspan), 5);
        inputs(:, 1) = 70 +  heaviside(tspan-100)*800; 
        inputs(:, 2) = 70 +  heaviside(tspan-200)*900;
        inputs(:, 3) = 70 +  heaviside(tspan-300)*1000;
        inputs(:, 4) = 70 +  heaviside(tspan-400)*1100;
        inputs(:, 5) = 1+heaviside(tspan-500)*10;
        end 


        function hfig  = MakePlotsFromDetailedSimData(detailedSimData)


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
        title('Geometry');

        subplot(4,2,2);
        pdeplot3D(thermalModel); 
        axis equal
        title('Mesh');
        xlabel('x');
        ylabel('y');

        subplot(4,2,3)
        plot(inputs(:,1:4), 'linewidth', 2)
        xlabel('time index')
        ylabel('temperature (°F)')
        title('Inputs')
        yyaxis right
        plot(inputs(:,5), 'linewidth', 2)
        ylabel('velocity (m/s)')
        legend({'T_1', 'T_2', 'T_3', 'T_4', 'velocity'}, 'Location', 'northeast')

        subplot(4,2,5)
        plot(T, 'linewidth', 1)
        ylabel('temperature (°F)')
        title('Nodal Tempereratures')
        xlabel('time index')

        subplot(4,2,7)
        plot(S, 'linewidth', 1)
        ylabel('stress (psi)')
        title('Nodal Stresses')
        xlabel('time index')

        subplot(4,2,6)
        xinterp = setup.geometry.xinterp;
        yinterp = setup.geometry.yinterp;
        [X, Y] = meshgrid(xinterp, yinterp);
        Trs = reshape( T(end, :), size(X, 1), []);
        colormap('jet')
        contourf(X, Y, Trs, 100, 'LineColor','none')
        colorbar()
        axis equal
        title('Final Temperatures')

        subplot(4,2,8)
        Srs = reshape( S(end, :), size(X, 1), []);
        colormap('jet')
        contourf(X, Y, Srs, 100, 'LineColor','none')
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

        end 


        function [dataTable, detailedSimData] =  RunPlateSim(inputs)
        disp('-STARTING--------------------------')
        setup = PlateSimFunctions.GetDefaultInputs();
        setup = PlateSimFunctions.ModifyLoadsAndICs(setup, inputs);
        thermalModel = PlateSimFunctions.MakeModel(setup);
        [thermalModel, thermalResults, inputs] = PlateSimFunctions.SolveTemperature(thermalModel, setup);
        T = PlateSimFunctions.InterpolateResults(thermalResults, setup);
        [S, structuralResults, structuralModel] = PlateSimFunctions.SolveStresses(thermalModel, thermalResults, setup);


        dataTable = table();
        dataTable.time = single(thermalResults.SolutionTimes');
        dataTable.Temperature = single(T');
        dataTable.Stress = single(S');
        dataTable.Inputs = single(inputs');
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
        disp('     * starting to solve structural model')

        E = setup.material.E;
        poisson = setup.material.poisson;
        alpha = setup.material.alpha;
        density = setup.material.rho;
        Tref = setup.material.Tref;

        structuralModel = createpde("structural","static-solid");
        structuralModel.Mesh = thermalModel.Mesh;
        structuralModel.Geometry = thermalModel.Geometry;
        structuralModel.ReferenceTemperature = Tref;  % referance temperature for 0 thermal expansion
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



        function setup = ModifyLoadsAndICs(setup, inputs)


        U = cell(5,1);
        tspan = (0:size(inputs,1)-1)*setup.solver.dt;

        for i = 1:size(inputs,2)
            U{i} = @(loc, state) interp1(tspan,inputs(:, i),state.time,'previous', 'extrap');
        end

        setup.loads.v = U{5};
        setup.loads.BCs = U(1:4);
        setup.solver.tend = tspan(end);
        disp('     * done modifying BCs and loads')

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
        setup.material.emiss = 0.25;                % emissivity of the plate surface
        setup.material.E = 2E9;          % youngs modulus
        setup.material.poisson = .3;       % Poisson ratio
        setup.material.alpha = 2.4e-5;       % coeff of thermal expansion
        setup.material.Tref = 70;             % reference temperature for thermal expansion

        % Geometry and mesh
        setup.geometry.width = 1;                  %  plate width, m
        setup.geometry.height = 1;                 %  plate height, m
        setup.geometry.thickness = .1;                 %  plate thickness, m
        setup.geometry.elSize = 0.25;               %  element size, m
        setup.geometry.xinterp = linspace(0, setup.geometry.width, 11);               %  solution interpolation points
        setup.geometry.yinterp = linspace(0, setup.geometry.height, 11);                    %  solution interpolation points
        setup.geometry.zinterp = setup.geometry.thickness/2;                                 %  solution interpolation points

        % solver
        setup.solver.tend = 5000;                 % end time (s)
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
        sig = setup.material.sig;

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
        generateMesh(thermalModel, 'Hmax',elSize);

        thermalProperties(thermalModel,'ThermalConductivity',k,...
                                       'MassDensity',rho,...
                                       'SpecificHeat',cp);

        thermalModel.StefanBoltzmannConstant = sig           ;      
        disp('     * done making thermal model')

        end 



        function [model, results, inputs] = SolveTemperature(model, setup)
        disp('     * starting to solve thermal model')

        Ta = setup.material.Ta;
        tend = setup.solver.tend;          
        dt = setup.solver.dt;    
        RelativeTolerance = setup.solver.RelativeTolerance ;
        AbsoluteTolerance = setup.solver.SolverOptions.AbsoluteTolerance ;
        U = setup.loads.BCs;
        Tinit = setup.ICs.Tinit ;
        hc = @(loc, state) setup.material.h(setup.loads.v(loc, state)) ;

        thermalBC(model,'Face',[2] ,'Temperature', U{1});
        thermalBC(model,'Face',[3],'Temperature', U{2});
        thermalBC(model,'Face',[4],'Temperature', U{3});
        thermalBC(model,'Face',[5],'Temperature', U{4});
        thermalBC(model,'Face',[1],'ConvectionCoefficient', hc, 'AmbientTemperature',Ta);
        thermalBC(model,'Face',[6],'ConvectionCoefficient', hc, 'AmbientTemperature',Ta);

        thermalIC(model,Ta);

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


        disp('     * done solving thermal model')

        end 

        
        
        

    end
end

