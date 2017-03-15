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

% Obtain analytical jacobian, necessary only for the first time, otherwise set to 0.
RUN_JACOBIAN = 0;
JACOBIAN_FLUX_SCHEME = 1;   % 1 for Van Leer, 2 for Roe

ROM_METHOD = 1;     % 1 for RK4, 2 for GN
CFL_ROM = 5;
ROM_span = 0.2;     % ending time of ROM
ML = 151;           % number of modes kept
RESTART = 0;        % 0 for starting ROM at the end of FOM, 1 for at t=0

if RUN_STEADY==1
    Solve_Steady_State
end

if RUN_FOM==1
    FOM(CFL_FOM,FOM_span)
end

if RUN_JACOBIAN==1
    switch JACOBIAN_FLUX_SCHEME
        case 1
            run jacobian_van_leer
        case 2
            run jacobian_roe
    end
end

switch ROM_METHOD
    case 1
        ROM(ROM_METHOD,CFL_ROM,ROM_span,ML,RESTART)
    case 2
        ROM(ROM_METHOD,CFL_ROM,ROM_span,ML,RESTART,JACOBIAN_FLUX_SCHEME)
end
