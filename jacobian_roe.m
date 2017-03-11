function jacobian_roe
run constantsch4
syms Q1 Q2 Q3 Q4 gam Q1L Q2L Q3L Q4L Q1R Q2R Q3R Q4R gammaL gammaR real;
F=[Q2;
    Q2.^2./Q1+(gam-1).*(Q3-Q2.^2/2./Q1);
    Q2./Q1.*(gam*Q3-(gam-1)*Q2.^2/2./Q1);
    Q4.*Q2./Q1
    ];
dFdQ = jacobian(F, [Q1;Q2;Q3;Q4]);

rL = Q1L;
uL = Q2L/rL;
YL = Q4L/rL;

ML = 1/(YL*(1/M_ox) + (1-YL)*(1/M_pr));
XL = YL*M_ox/ML;
CpL = XL*Cp_ox + (1-XL)*Cp_pr;
gammaL = CpL/(CpL-(Rgas/ML));
pL = (gammaL - 1).*(Q3L - 0.5*Q2L*uL);

rHL = Q3L + pL;
HL = rHL/rL;

rR = Q1R;
uR = Q2R/rR;
YR = Q4R/rR;

MR = 1/(YR*(1/M_ox) + (1-YR)*(1/M_pr));
XR = YR*M_ox/MR;
CpR = XR*Cp_ox + (1-XR)*Cp_pr;
gammaR = CpR/(CpR-(Rgas/MR));
pR = (gammaR - 1).*(Q3R - 0.5*Q2R*uR);

rHR = Q3R + pR;
HR = rHR/rR;

di     = sqrt(Q1R/Q1L);
d1     = 1.0./(1.0+di);

ri     = simplify(sqrt(rL*rR));
ui     = simplify((di*uR + uL)*d1);
Hi     = simplify((di*HR + HL)*d1);
Yi     = simplify((di*YR + YL)*d1);
Mi = simplify(1/(Yi*(1/M_ox) + (1-Yi)*(1/M_pr)));
Xi = simplify(Yi*M_ox/Mi);
Cpi = simplify(Xi*Cp_ox + (1-Xi)*Cp_pr);
gammai = simplify(Cpi/(Cpi-(Rgas/Mi)));
gmi = gammai-1.0;
pi     = simplify(gmi/gammai*ri*Hi);

af     = 0.5*ui*ui;
ucp    = simplify(ui);
c2     = gmi.*(Hi - af);           % speed of sound squared
ci     = simplify(sqrt(c2));
ci1    = 1.0/ci;

% eigenvalues
l1 = ucp+ci;
l2 = ucp-ci;
l3 = ucp;

l1 = abs(l1); l2 = abs(l2); l3 = abs(l3);

% average and half-difference of 1st and 2nd eigs
s1    = 0.5*(l1 + l2);
s2    = 0.5*(l1 - l2);

% left eigenvector product generators
G1    = simplify(gmi*(af*(Q1R-Q1L) - ui*(Q2R-Q2L) + (Q3R-Q3L)));
G2    = simplify(-ucp*(Q1R-Q1L)+(Q2R-Q2L));

% required functions of G1 and G2
C1    = G1*(s1-l3)*ci1*ci1 + G2*s2*ci1;
C2    = G1*s2*ci1          + G2*(s1-l3);
 
% flux assembly
Fi1 = l3*(Q1R-Q1L) + C1;
Fi2 = l3*(Q2R-Q2L) + C1*ui + C2;
Fi3 = l3*(Q3R-Q3L) + C1*Hi + C2*ucp;
Fi4 = l3*(Q4R-Q4L);

Fi = [Fi1;Fi2;Fi3;Fi4];
dFidQL = jacobian(Fi, [Q1L;Q2L;Q3L;Q4L]);
dFidQR = jacobian(Fi, [Q1R;Q2R;Q3R;Q4R]);

dFdQ = func2str(matlabFunction(dFdQ));
dFidQL = func2str(matlabFunction(dFidQL));
dFidQR = func2str(matlabFunction(dFidQR));

save('./data/dFdQ_roe.mat','dFidQL','dFidQR','dFdQ');
end
