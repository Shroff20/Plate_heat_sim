clear all
close all
clc

output_dir = 'C:\Users\ssmee\Google Drive\Sameer''s Files\Jobs\Dynamics_and_machine_learning\Plate_NN\Plate_heat_sim\step_jobs';
N_jobs = 100;  % number of step jobs to create and run
N_timepoints_ss = 100;  % number of points until the system reaches steady state (will hold initial conditions this long, and allow transient to run this long)


%% make step inputs

T_step_config = table();
T_step_config.input = ["T1"; "T2"; "T3"; "T4"; "v"];
T_step_config.min = [100; 100; 100; 100; 0];
T_step_config.max = [1000; 1000; 1000; 1000; 10];
T_step_config.step_mag_max = [1000; 1000; 1000; 1000; 10]/2;
disp(T_step_config)
Ninputs = height(T_step_config);

T_jobs = table();
T_jobs.casename = strings(N_jobs, 1);
T_jobs.U = cell(N_jobs, 1);
T_jobs.step_var = nan(N_jobs,1);
T_jobs.step_mag = nan(N_jobs, 1);


for j = 1:N_jobs
    
    U0 = T_step_config.min + rand(Ninputs, 1).*(T_step_config.max - T_step_config.min);
    step_var = randi(Ninputs);
    step_mag = zeros(Ninputs, 1);
    step_mag(step_var) = T_step_config.step_mag_max(step_var)*rand();

    U = U0*ones(1, N_timepoints_ss*2);
    U(:,N_timepoints_ss:end) = U(:,N_timepoints_ss:end) + step_mag;
    
    T_jobs.U{j} = U';
    T_jobs.step_var(j) = step_var;
    T_jobs.step_mag(j) = step_mag(step_var);
    T_jobs.casename(j) = "step_" + num2str(j);
    
end 

disp(T_jobs)

%% run jobs step inputs

for j = 1:height(T_jobs)
    disp("STARTING JOB: " + T_jobs.casename(j) + " =======================================================")
    PlateSimFunctions.RunSimPipeline( T_jobs.U{j}, T_jobs.casename(j), output_dir)
end 