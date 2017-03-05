close all, clear all

mkdir ./data
mkdir ./images
addpath('./data')
addpath('./images')

% Obtain data for steady state, necessary only for the first time, otherwise set to 0.
RUN_STEADY = 0;

% Run FOM to generate snapshot for the first time or when FOM_span/CFL_FOM is changed,
% otherwise set to 0.
RUN_FOM = 0;
CFL_FOM = 1;
FOM_span = 0.1;     % ending time of FOM

% Obtain analytical jacobian needed by Gauss-Newton method, necessary only for the first time, otherwise set to 0.
RUN_JACOBIAN = 0;

if RUN_STEADY==1
    run Solve_Steady_State
end

if RUN_FOM==1
    run FOM
end

if RUN_JACOBIAN==1
    jacobian_van_leer
end

clear all

ROM_METHOD = 2;     % 1:RK4, 2: Gauss-Newton 
CFL_ROM = 1;
ROM_span = 0.2;     % ending time of ROM
ML = 151;           % number of modes kept
RESTART = 1;        % 0 for starting ROM at the end of FOM, 1 for at t=0  

ROM(ROM_METHOD,CFL_ROM,ROM_span,ML,RESTART)