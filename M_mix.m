function M = M_mix(M_ox,M_pr,Y)
% Computes the mixture molecular weight

M1 = Y*(1/M_ox) + (1-Y)*(1/M_pr);
M = 1./M1;
