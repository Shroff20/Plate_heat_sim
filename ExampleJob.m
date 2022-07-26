clear all
close all
clc


output_dir = 'C:\Users\ssmee\Google Drive\Sameer''s Files\Jobs\Dynamics_and_machine_learning\Plate_NN\Plate_heat_sim';
casename = 'exmaple_output';
inputs = PlateSimFunctions.GetExampleInput();  % timepoints x 5 matrix of inputs
% columns are [T1, T2, T3, T4, velocity]  where T# is the plate edge
% temperature , and velocity is the fluid velocity for convection

% Limits:  T = [100, 1000]
%          v = [0, 10]

PlateSimFunctions.RunSimPipeline(inputs, casename, output_dir)


