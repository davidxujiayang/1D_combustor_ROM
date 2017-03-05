function gam = gamma_mix(M,Y)

run 'constantsch4.m'

X = Y*M_ox./M;

Cp = X*Cp_ox + (1-X)*Cp_pr;

gam = Cp./(Cp-(Rgas./M));

end