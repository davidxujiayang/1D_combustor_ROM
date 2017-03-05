%% constants

% OUTPUT OF CANTERA

Rgas = 8.314e3;                  % Gas constant: the mixture molecular weight is expressed in g/mol hence the 1e3 factor

minf = 0.32;                    % mass flow rate
mf = 0.027;
mout = minf+mf;

h_f = 5e7;                 % Heat of formation of the combustion products 
            
M_ox  = 22.113;                   % Oxidizer molecular weight (42% O, 58% H20)
Cp_ox = 1.797e3;                 % Oxidizer mass specific heat
Cp_f  = 4.668e3;

M_pr  = 21.337;                     % Combustion product molecular weight
Cp_pr = 2.608e3;                   % Combustion product specific heat

f2ox = 0.424*0.25;    % equivalent fuel to oxidizer ratio
Yout = (minf-mf/f2ox)/(minf+mf);
phi = mf/f2ox/minf;     % equivalent ratio

% h_ox_ref = -6.54588e+006;
% h_pr_ref = -1.09064e+007;
% T_ref  = 1030;                   % Reference temperature
% T_ref2 = 2430;

% A_ox  = -Cp_ox*T_ref;            % h_ox = Cp_ox * T + A_ox 
% A_pr  =  h_f - Cp_pr*T_ref;       % h_pr = Cp_pr * T + A_pr

% dC = Cp_ox - Cp_pr;             
% dA = A_ox - A_pr;