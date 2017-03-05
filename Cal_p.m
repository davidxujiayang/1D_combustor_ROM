function p = Cal_p(U)
run 'constantsch4.m'
U = reshape(U,[],4);
r = U(:,1);
u = U(:,2)./r;
Y = U(:,4)./r;
M = M_mix(M_ox,M_pr,Y);
gamma = gamma_mix(M,Y);
p = (gamma-1).*(U(:,3) - 0.5*U(:,2).*u);
end