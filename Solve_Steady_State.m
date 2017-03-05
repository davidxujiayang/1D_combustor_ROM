function Solve_Steady_State

%% Chamber Geometry
N = 1200;            % Node Number
Li = 0.1397;          % Injector Length

[x,A] = read_shape(N);
xelem = (x(2:end)+x(1:end-1))/2;
Aelem = (A(2:end)+A(1:end-1))/2;

% source term function
ls = 0.1334-Li; lf = 0.2096-Li; % see EUCASS paper by Dr. Frezzotti
s = (1 + sin(-pi/2 + 2*pi*(xelem - ls)/(lf - ls))).*(xelem < lf).*(xelem > ls);
Int = trapz(xelem,s);

% number of elements, interior edges, and boundary edges
nelem = length(x) - 1;
nIE   = nelem - 1;
nBE   = 2;

% Cell volumes
V = Aelem .* (x(2:end)-x(1:end-1));

%% Gas Properties

run 'constantsch4.m'

%% Inlet States

Yinf = 1;
gamma_inf = gamma_mix(M_ox,Yinf);
pinf = 1.44e6;
Tinf = 1030;

rinf = pinf*M_ox/(Rgas*Tinf);
uinf = minf/rinf/A(1);        % inlet velocity
Stateinf = [Yinf,gamma_inf,pinf];
Propinf = [Yinf,M_ox,gamma_inf,pinf];

%% Initialize

Uinf(1) = rinf;
Uinf(4) = rinf*Yinf;
Uinf(2) = minf/Aelem(1);
M = M_mix(M_ox,M_pr,Yinf);
Uinf(3) = 0.5*uinf*Uinf(2)+pinf/(gamma_inf-1);
Uout = zeros(1,4);

U = repmat(Uinf, nelem, 1);

% % Load backup data
% load('U_steady.mat');

if length(U(:,1))~=nelem
    Nold = length(U(:,1))+1;
    xold = read_shape(Nold);
    xelem_old = (xold(2:end)+xold(1:end-1))/2;
    for i = 1:4
        Unew(:,i)=interp1q([xold(1);xelem_old(2:end-1);xold(end)],U(:,i),xelem);
    end
    U = Unew;
end

Std.Li = Li; Std.A = A; Std.xelem = xelem; Std.Aelem = Aelem; Std.V = V; Std.s = s; Std.Int = Int;
Std.Uinf = Uinf; Std.Stateinf=Stateinf; Std.Propinf = Propinf;
save('./data/Std.mat', 'Std');

%% convergence settings

n = 0;
nmax = 50000;           % Max number of iterations
Rtol = -2;              % tolerance
CFL = 0.5;              % CFL number
Res = [];
Rmax = 1;
Rhist = [];

%% iteration

while (Rmax > Rtol)
    
    n = n+1;
    
    r = U(:,1);
    u = U(:,2)./r;
    m = U(:,2).*Aelem;
    Y = U(:,4)./r;
    M = M_mix(M_ox,M_pr,Y);
    gamma = gamma_mix(M,Y);
    p = (gamma-1).*(U(:,3) - 0.5*U(:,2).*u);
    T = p.*M./r/Rgas;
    
    Prop = [Y,M,gamma,p]; % nelem*3 gas property matrix
    
    % residuals
    R  = zeros(nelem, 4);
    sdA = zeros(nelem, 1); % the product of abs(s) and Area
    
    eL = 1:nIE;
    eR = 2:nIE+1;
    
    UL = U(eL,:);
    UR = U(eR,:);
    PropL=Prop(eL,:);
    PropR=Prop(eR,:);
    
    [F, smax] = FluxFunction(UL, UR, PropL,PropR);
    
    FA = F.*repmat(A(eR),1,4);
    sdAchange = smax.*A(eR);
    
    R(eL,:) = R(eL,:) + FA;
    R(eR,:) = R(eR,:) - FA;
    sdA(eL) = sdA(eL) + sdAchange;
    sdA(eR) = sdA(eR) + sdAchange;
    
    % boundary conditions: supersonic outlet -> extrapolate from the two adjacent interior cells

    PropL = Prop(end,:);
    Propout = 2*Prop(end-1,:) - Prop(end-2,:); % Propout = [Yout,Mout,gammaout,pout]
    uout = 2*u(end) - u(end-1);
    Tout = 2*T(end) - T(end-1);
    rout = Propout(4)/Tout*Propout(2)/Rgas;
    Uout(1) = rout;
    Uout(2) = rout*uout;
    Uout(3) = 0.5*uout*Uout(2)+Propout(4)/(Propout(3)-1);
    Uout(4) = rout*Propout(1);
    
    [F, smax] = FluxFunction(U(end,:), Uout, PropL,Propout);
    R(end,:)=R(end,:) + F*A(end);
    sdA(end) = sdA(end) + smax*A(end);
        
    % boundary conditions: inlet
    
    PropR=Prop(1,:);
    [F, smax] = FluxFunction(Uinf,U(1,:), Propinf ,PropR);
    R(1,:)=R(1,:) - F*A(1);
    sdA(1) = sdA(1) + smax*A(1);
    
    % source term
    S=zeros(nelem, 4);
    S(:,2)=p.*(A(2:end)-A(1:end-1))./V;
    
    Sf=zeros(nelem, 4);
    
    w_f = m(1)*phi*f2ox/Int*s./Aelem;
    w_ox = w_f/f2ox;
    
    Sf(:,1) = w_f;
    Sf(:,2) = w_f.*u;
    Sf(:,3) = w_f.*(h_f + Cp_pr * 298) + w_ox.*(Cp_pr - Cp_ox).* 298;
    Sf(:,4) = -w_ox;
    
    Ssum = S+Sf;
    
    R = R - Ssum.*repmat(V,1,4);
    
    % local timestep
    dt = 2*CFL*V./sdA;
    
    % Calculate maximum residual
    Rmax = log10(max(max(abs(R))));
    
    % Forward Euler
    
    U = U - repmat(dt,1,4).*R./repmat(V,1,4);
    
    % skip logging some iterations
    if (mod(n,1000) ~= 1),	continue;    end
        
    % print out residual
    fprintf(1, 'n = %d, log10(Rmax) = %.10f,\n', n, Rmax);
    fprintf(1, 'm1 = %d, m_max = %d \n', m(1), max(m));
    fprintf(1, 'p1 = %d, p_max = %d \n', p(1), max(p));
    fprintf(1, 'T1 = %d, T_max = %d \n', T(1), max(T));
    
    Rhist = [Rhist; Rmax];
    % write backup state
    save('./data/U_steady.mat', 'U');
end
save('./data/U_steady.mat', 'U');
Std.dt_min = min(dt); Std.CFL = CFL; Std.p = p;
save('./data/Std.mat', 'Std');

figure
[hAx,hLine1,hLine2]=plotyy(xelem,T,xelem,p/1e6);
xlabel('Location (m)','fontsize',15)
ylabel(hAx(1),'Temperature (K)','fontsize',15) % left y-axis
ylabel(hAx(2),'Pressure (MPa)','fontsize',15) % right y-axis
grid on, grid minor
hLine1.LineWidth = 1;
hLine2.LineWidth = 1;
set(gcf,'position',[100 100 600 200])
print(gcf,'-djpeg',sprintf('-r%d',300),'./images/steady_state.jpg');

end
