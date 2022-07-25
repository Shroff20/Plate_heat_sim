clear all
close all
clc


output_dir = 'C:\Users\ssmee\Google Drive\Sameer''s Files\Jobs\Dynamics_and_machine_learning\Plate_NN\Plate_heat_sim';
casename = 'test3';
inputs = PlateSimFunctions.GetExampleInput();  % timepoints x 5 matrix of inputs
% columns are [T1, T2, T3, T4, velocity]  where T# is the plate edge
% temperature , and velocity is the fluid velocity for convection

% Limits:  T = [100, 1000]
%          v = [1, 10]


PlateSimFunctions.RunSimPipeline(inputs, casename, output_dir)


%% step inputs
N = 5 % steps in each input (N^5 cases)

T_step_config = table();
T_step_config.input = ["T1"; "T2"; "T3"; "T4"; "v"];
T_step_config.min = [100; 100; 100; 100; 0];
T_step_config.max = [1000; 1000; 1000; 1000; 10];
T_step_config.step_mag = [1000; 1000; 1000; 1000; 10]/10;


disp(T_step_config)

elements = {};
for i = 1:height(T_step_config)
    elements{i} = linspace(T_step_config.min(i), T_step_config.max(i), N);
end 
combinations = cell(1, numel(elements)); %set up the varargout result
[combinations{:}] = ndgrid(elements{:});
combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false); %there may be a better way to do this
result = [combinations{:}]; % NumberOfCombinations by N matrix. Each row is unique.

T_initial_conditions = table();
for i = 1:height(T_step_config)
    T_initial_conditions.(T_step_config.input(i)) = result(:,i);
end 
disp(T_initial_conditions)

