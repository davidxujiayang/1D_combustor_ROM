function [R,p] = Cal_Res(U)
global nelem nIE A V Uinf Propinf Sf V4 phist pstd nread n Uptb Uptb2 Uptb3 w_f alpha t tptb
run 'constantsch4.m'

r = U(:,1);
u = U(:,2)./r;
Y = U(:,4)./r;
M = M_mix(M_ox,M_pr,Y);
gamma = gamma_mix(M,Y);
p = (gamma-1).*(U(:,3) - 0.5*U(:,2).*u);
T = p.*M./r/Rgas;
Prop = [Y,M,gamma,p]; % nelem*4 gas property matrix

% residuals
R  = zeros(nelem, 4);

eL = 1:nIE;
eR = 2:nIE+1;
UL = U(eL,:);
UR = U(eR,:);
PropL=Prop(eL,:);
PropR=Prop(eR,:);

[F, ~] = FluxFunction(UL, UR, PropL,PropR);
FA = F.*repmat(A(eR),1,4);

R(eL,:) = R(eL,:) + FA;
R(eR,:) = R(eR,:) - FA;

% boundary conditions: supersonic outlet -> extrapolate from the two
% adjacent interior cells

PropL = Prop(end,:);
Propout = 2*Prop(end-1,:) - Prop(end-2,:); % Propout = [Yout,Mout,gammaout,pout]
uout = 2*u(end) - u(end-1);
Tout = 2*T(end) - T(end-1);
rout = Propout(4)/Tout*Propout(2)/Rgas;
Uout(1) = rout;
Uout(2) = rout*uout;
Uout(3) = 0.5*uout*Uout(2)+Propout(4)/(Propout(3)-1);
Uout(4) = rout*Propout(1);

[F, ~] = FluxFunction(U(end,:), Uout, PropL,Propout);
R(end,:)=R(end,:) + F*A(end);

% boundary conditions: inlet

PropR=Prop(1,:);
Uptb(2) = Uptb2(n);
Uptb(3) = Uptb3(n);
[F, ~] = FluxFunction(Uinf + Uptb,U(1,:), Propinf,PropR);
R(1,:)=R(1,:) - F*A(1);

% source term
S=zeros(nelem, 4);
S(:,2)=p.*(A(2:end)-A(1:end-1))./V;

Sf(:,2) = w_f.*u;

Sq=zeros(nelem, 4);
Sq(:,3) = alpha*(t(n)>tptb)*Sf(:,3).*(phist(:,nread)./pstd-1);
Ssum = S+Sf+Sq;

R = R - Ssum.*V4;
end