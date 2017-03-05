global nelem nIE A V Uinf Propinf Sf V4 phist pmean pstd nread n Uptb Uptb2 Uptb3 w_f alpha t tptb Phi Ustd INDEX eL eR

%% Load steady results
run Load_steady
Ustd = reshape(U,[],1);

dt_FOM = CFL_FOM/CFLstd * dtstd;

eL=[]; eR=[];
for i=1:4
    INDEX{i} = (i-1)*nelem+1:i*nelem;
    eL = [eL,(i-1)*nelem + (1:nIE)];
    eR = [eR,(i-1)*nelem + (2:nIE+1)];
end

%% Unsteady heat release parameters

alpha = 3.6; Ct = 1;
F1= 1400; % 1st resonant frequency;
T1 = 1/F1; % 1st resonant period;
tau = Ct*T1;

n_FOM = ceil(FOM_span/dt_FOM)+1;
dn_FOM = ceil(tau/dt_FOM);
t = (0:n_FOM-1)*dt_FOM;

%% Perturbations & steady source term

tptb = 0.008;         % perturbation lasting time
Uptb2 = (t<=tptb)*0.001.*sin(2*pi*F1*t)*Uinf(2);   % a small perturbation term is added at the beginning;
Uptb3 = 2*Uptb2*Uinf(2)/Uinf(1);
Uptb = zeros(1,4);

Sf=zeros(nelem, 4);
w_f = mstd*phi*f2ox/Int*s./Aelem;
w_ox = w_f/f2ox;
Sf(:,1) = w_f;
Sf(:,3) = w_f.*(h_f + Cp_pr * 298) + w_ox.*(Cp_pr - Cp_ox).* 298;
Sf(:,4) = -w_ox;

%% FOM
V4 = repmat(V,1,4);
n = 0;
phist = zeros(nelem,dn_FOM);               % for allocation
pMonhist = zeros(1,n_FOM);
Usol = [];
% load Usol_FOM
while (n < n_FOM)
    n = n+1;
    
    nwrite = mod(n-1,dn_FOM)+1;
    nread = mod(n,dn_FOM)+1;
    
    % RK4
    [R,p] = Cal_Res(U);
    U1 = U - dt_FOM/4*(R./V4);
    phist(:,nwrite) = p;
    pMonhist(n) = p(nMon);
    
    R = Cal_Res(U1);
    U2 = U - dt_FOM/3*(R./V4);
    
    R = Cal_Res(U2);
    U3 = U - dt_FOM/2*(R./V4);
    
    R = Cal_Res(U3);
    U = U - dt_FOM*(R./V4);
    
    % skip logging some iterations
    if (mod(n,100) ~= 1), continue; end
    Utemp = reshape(U,4*nelem,1);
    Usol = [Usol Utemp];
    
    plot(t(1:20:n-1),pMonhist(1:20:n-1) - pMonhist(1));xlabel('time (s)');ylabel('p''(Pa)');
    title(['FOM, t = ' num2str(t(n))]);drawnow;
end
print(gcf,'-djpeg',sprintf('-r%d',300),'./images/FOM_result.jpg')
save('./data/Usol_FOM','n','nwrite','U','Usol','phist','pMonhist','CFL_FOM','FOM_span');

