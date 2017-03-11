function dRdQ = jacoRQ(a,scheme)
if nargin < 2
    scheme =1;
end

global nelem nIE A V Uinf Propinf w_f Phi Ustd Uptb Uptb2 Uptb3 n INDEX eL eR
run 'constantsch4.m'
U = Phi*a+Ustd;
r = U(INDEX{1});
u = U(INDEX{2})./r;
Y = U(INDEX{4})./r;
M = M_mix(M_ox,M_pr,Y);
gamma = gamma_mix(M,Y);
p = (gamma-1).*(U(INDEX{3}) - 0.5*U(INDEX{2}).*u);
T = p.*M./r/Rgas;
Prop = [Y,M,gamma,p]; % nelem*4 gas property matrix

% residuals
dRdQ  = zeros(4*nelem,4*nelem);

UL = U(eL);
UR = U(eR);
gammaL=gamma(1:nIE,:);
gammaR=gamma(2:nIE+1,:);

switch scheme
    case 1
        [FL,FR] = jacoFlux_vl(UL, UR, gammaL,gammaR, A(2:nIE+1));
    case 2
        [FL,FR] = jacoFlux_roe(UL, UR, gammaL,gammaR, A(2:nIE+1));
end

dRdQ(eL,eL) = dRdQ(eL,eL) + FL;
dRdQ(eL,eR) = dRdQ(eL,eR) + FR;
dRdQ(eR,eL) = dRdQ(eR,eL) - FL;
dRdQ(eR,eR) = dRdQ(eR,eR) - FR;

% boundary conditions: supersonic outlet -> extrapolate from the two
% adjacent interior cells

gammaL = gamma(end);
Propout = 2*Prop(end-1,:) - Prop(end-2,:); % Propout = [Yout,Mout,gammaout,pout]
gammaout = Propout(3);
uout = 2*u(end) - u(end-1);
Tout = 2*T(end) - T(end-1);
rout = Propout(4)/Tout*Propout(2)/Rgas;
Uout(1,1) = rout;
Uout(2,1) = rout*uout;
Uout(3,1) = 0.5*uout*Uout(2)+Propout(4)/(Propout(3)-1);
Uout(4,1) = rout*Propout(1);

switch scheme
    case 1
        [FL,~] = jacoFlux_vl(U([nelem,2*nelem,3*nelem,4*nelem]), Uout, gammaL,gammaout,A(end));
    case 2
        [FL,~] = jacoFlux_roe(U([nelem,2*nelem,3*nelem,4*nelem]), Uout, gammaL,gammaout,A(end));
end

dRdQ([nelem,2*nelem,3*nelem,4*nelem],[nelem,2*nelem,3*nelem,4*nelem]) = dRdQ([nelem,2*nelem,3*nelem,4*nelem],[nelem,2*nelem,3*nelem,4*nelem]) + FL;

% boundary conditions: inlet

gammaR=gamma(1);
Uptb(2) = Uptb2(n);
Uptb(3) = Uptb3(n);

switch scheme
    case 1
        [~,FR] = jacoFlux_vl((Uinf + Uptb)',U([1,1+nelem,1+2*nelem,1+3*nelem]), Propinf(3) ,gammaR,A(1));
    case 2
        [~,FR] = jacoFlux_roe((Uinf + Uptb)',U([1,1+nelem,1+2*nelem,1+3*nelem]), Propinf(3) ,gammaR,A(1));
end

dRdQ([1,1+nelem,1+2*nelem,1+3*nelem],[1,1+nelem,1+2*nelem,1+3*nelem]) = dRdQ([1,1+nelem,1+2*nelem,1+3*nelem],[1,1+nelem,1+2*nelem,1+3*nelem]) - FR;

% source term
for i = 1:nelem
    dRdQ([i,i+nelem,i+2*nelem,i+3*nelem],:) = dRdQ([i,i+nelem,i+2*nelem,i+3*nelem],:)/V(i);
    C = (A(i+1)-A(i))/V(i)*(gamma(i)-1);
    dRdQ(i+nelem,i) = dRdQ(i+nelem,i) - C/2*(U(i+nelem)/U(i))^2 + w_f(i)*U(i+nelem)/U(i)^2;
    dRdQ(i+nelem,i+nelem) = dRdQ(i+nelem,i+nelem) + C*(U(i+nelem)/U(i)) - w_f(i)/U(i);
    dRdQ(i+nelem,i+2*nelem) = dRdQ(i+nelem,i+2*nelem) - C;
end

dRdQ=dRdQ*Phi;
end