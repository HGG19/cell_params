% Open sim_batter_system.prj before executing this script!
clearvars;
clc;

time_and_current = csvread('soc_current.csv');
time_and_current = [time_and_current(1:end, 1), time_and_current(1:end, 3)];

disp('Load system parameters.')
run parameters/cell_parameters_SLPBB042126HN.m              % Load cell parameters
run parameters/system_parameters.m                          % Load system parameters
t_init = 'parameters/system_initialization.m';              % Data with initial state for electrical and thermal simulation
t_ambient = 'parameters/ambient_temperature_parameters.m';  % Data specifying the ambient temperature around the cells and the heat transfer coefficients. Needed for all thermal simulations.


% Different options for post-processing of simulation results and settings.

SimPara.LoggingOutput    = true;        % Log all variables specified in 'sim_battery_system.slx/Monitoring and Logging/Logging of Variables' in timeseries format
SimPara.OutputDataType   = 'double';    % Data type of simulation output. Note: Simulation is always done in double precision
SimPara.OutputDecimation = 100;         % Log (and therefore plot and save) only every n-th value.

SimPara.PlotResults      = true;        % Plot all available data (everything that gets logged in the model). If set 'false' there are no plot, or, if plots are saved they are closed after saving.
SimPara.SavePlots        = true;       % Save all plots in a folder with current time in ..\simulation_results (.fig format)

SimPara.SaveSettings     = true;       % Save generated parameter structs in the save file
SimPara.SaveResults      = true;        % Save simulation results in the save file

SimPara.ClearWorkspace   = false;       % Clear workspace after simulation (and saving)


%% Thermal system simulation (constant ambient temperatures)

% Enable or disable thermal simulation by setting "true" or "false". 
% If disabled: T_cell remains constant with the values specified in 'Initial State'.

% Those values are defined in 'system_parameters.m'
% Hint: Enable all thermal system variables with dummy values to avoid errors.

SimPara.thermal_sim_enable = false;  


%% Geometric system information

% Only needed and called if advanced thermal sim or temp sensor submodule is used.

t_geometric_info = 'parameters/geometric_system_parameters.m';


%% Advanced thermal simulation

% Use advanced simulation model to consider thermal interaction between
% cells. If this is set to "false" only the simple thermal model is used.

% If you set this to 'true' you must provide additional thermal
% parameters in thermal_system_parameters.m! Because additional geometric
% information about the system is needed SysPara.p and SysPara.s from
% 'system_parameters.m' get overwritten if this system is activated!

% Note: If activated this has a significant impact on execution time!

SimPara.heat_exchange_enable = false;        % Turn heat exchange between cells on/off

%run parameters/heat_transfer_parameters.m   % Parameters for heat transfer submodel


%% Use temperature sensors submodule

% Simulate temperature sensors in the battery system. The temperature values
% of the sensors are written to 'T_sensors'. Refer to the script called 
% below for further information.

% Because additional geometric information about the system is needed 
% SysPara.p and SysPara.s from 'system_parameters.m' get overwritten if 
% this system is activated!

SimPara.TempSensors_enable = false;     % Turn temp. sensors on/off
run parameters/TempSens_parameters.m    % Parameters for temperature sensors


%% Load Spectrum Analysis of simulation data

% Perform Load Spectrum Analysis (LSA) on the simulation results.

% This is done during runtime and will massively decrease the amount of
% logging data, if you turn the timeseries logging of. You find the LSA
% subsystem here: sim_battery_system.slx/Monitoring and Logging/

% Spectrum Analysis
% Comment out all values including the subsystems you don't need.

% Note: Especially the three parametric LSA significantly add to the
% computational load for large battery systems! Test this impact before you
% turn on everything!
SimPara.LoadSpectra_enable = false;     % Turn LSA on/off
run parameters/LSA_parameters.m         % Parameters for Load Spectrum Analysis


%% System initialization
disp('Load initial system values.')
run(t_init);
run(t_ambient);


%% Calculate the statistical parameter deviation for each cell in the system
disp('Calculate statistical parameter deviations.')
run subroutines/deviations_stat.m

%% Run simulation
disp('Battery system simulation.')


SimPara.OutputDataType   = 'double';
SimPara.ClearWorkspace   = true;
SimPara.LoggingOutput    = true; 

all_load_cycles = [time_and_current(1:end, 1), zeros(size(time_and_current, 1), 1), -time_and_current(1:end, 2)]; % Load cycle format [time [s], volatege [V], current [A]]. Negative values for discharge :)

all_SOC_data = [];
SysPara.BatStateInit.electrical.SOC = 1;
for i = 1:size(all_load_cycles, 1)
    SimPara.t_step  = 0.1;     % Set simulation step-size in s. Note the simulation may get unstable, if you use RC-elements with small time constants.
    SimPara.t_sim   = 10;   % Set end time of simulation in s
    load_cycle = [all_load_cycles(i, 1), all_load_cycles(i, 2), all_load_cycles(i, 3)];
    sim('sim_battery_system');
    run stop_function
    disp("Iteration: " + i + ", Time: " + load_cycle(1) + ", Current: " + load_cycle(3) + ", SOC: " + SOC.Data)
    SysPara.BatStateInit.electrical.SOC = SOC.data;
end
time_data = SOC.Time;
soc_data = SOC.Data;
soc_data = reshape(soc_data, [size(soc_data, 3), 1]);
time_data = reshape(time_data, [size(time_data, 1), 1]);
time_and_soc = [time_data, soc_data];
all_SOC_data = [all_SOC_data; time_and_soc];
disp(all_SOC_data)
%SysPara.BatStateInit.electrical.SOC = all_SOC_data(end, 2);
%disp(all_SOC_data)


disp(SysPara.BatStateInit.electrical)


%% Plot
run subroutines/get_save_folder_and_time.m
run subroutines/saving.m
run subroutines/plotting.m


clearvars t_* % Delete temporary variables
if SimPara.ClearWorkspace == true
    clearvars
end
disp('Simulation run finished.')
