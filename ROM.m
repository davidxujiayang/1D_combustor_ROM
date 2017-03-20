function ROM(ROM_METHOD,CFL_ROM,ROM_span,ML,RESTART,scheme)
% ROM_METHOD: RK4 = 1; GN = 2;
method = {',RK4',',GN'};
if nargin<5
    RESTART = 0;
end
if nargin<6
    scheme = 1;
end

global nelem nIE A V Uinf Propinf Sf V4 phist pmean pstd nread n Uptb Uptb2 Uptb3 w_f alpha t tptb Phi Ustd INDEX eL eR dFdQ dJdQ2 dJdQ3 dJdQ4

%% Load steady results
run Load_steady

Ustd = reshape(U,[],1);
V4 = repmat(V,1,4);

load Usol_FOM
dt_FOM = CFL_FOM/CFLstd * dtstd;
dt_ROM = CFL_ROM/CFLstd * dtstd;
if ML>nelem
    ML = nelem;
end

if ROM_METHOD == 2
    switch scheme
        case 1
            run Load_Jacobian_vl
        case 2
            run Load_Jacobian_roe
    end
end

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


if RESTART == 1
    n_ROM = ceil((ROM_span)/dt_ROM)+1;
    dn_ROM = ceil(tau/dt_ROM);
    t = (0:n_ROM-1)*dt_ROM;
    n_total = n_ROM;
else
    n_FOM = ceil(FOM_span/dt_FOM)+1;
    dn_FOM = ceil(tau/dt_FOM);
    t = (0:n_FOM-1)*dt_FOM;
    n_ROM = ceil((ROM_span-FOM_span)/dt_ROM)+1;
    dn_ROM = ceil(tau/dt_ROM);
    t = [t t(end)+(1:n_ROM)*dt_ROM];
    n_total = n_FOM+n_ROM;
end


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

%% Convert
% basis construction
[~,Nsnap] = size(Usol);
Usol_uns = Usol - repmat(Ustd,1,Nsnap);

a = zeros(4*ML,1);
U = reshape(U,[],1);
U_uns = U-Ustd;
for i=1:4
    I1 = INDEX{i};
    I2 = (i-1)*ML+1:i*ML;
    
    [svdU,~,~] = svd(Usol_uns(I1,:));
    Phi(I1,I2) = svdU(:,1:ML);
    
    if RESTART == 1, continue, end
    
    svd_a = svdU.'*U_uns(I1);
    a(I2) = svd_a(1:ML);
    U_uns(I1) = Phi(I1,I2)*a(I2);
end
PhiT=Phi.';
phist_old = phist;
nwrite_old = nwrite;

%% ROM
phist = zeros(nelem,dn_ROM);               % for allocation

if RESTART == 1
    n = 0;
else
    nwrite_new = mod(n-1,dn_ROM)+1;
    base_old = (linspace(0,1,dn_FOM))';
    base_new = (linspace(0,1,dn_ROM))';
    for i=1:nelem
        phist_temp=circshift(phist_old(i,:)',[dn_FOM-nwrite_old,0]);
        phist(i,:)=(interp1q(base_old,phist_temp,base_new))';
        phist(i,:)=circshift(phist(i,:),[0,nwrite_new-dn_ROM]);
    end
end

pmean = pstd;
namejpg = ['./images/FOM_span=',num2str(FOM_span),',ML=',num2str(ML),method{ROM_METHOD},',CFL=',num2str(CFL_ROM),'.jpg'];
namefig = ['./images/FOM_span=',num2str(FOM_span),',ML=',num2str(ML),method{ROM_METHOD},',CFL=',num2str(CFL_ROM),'.fig'];
namemat = ['./data/FOM_span=',num2str(FOM_span),',ML=',num2str(ML),method{ROM_METHOD},',CFL=',num2str(CFL_ROM),'.mat'];
%load(namemat);
n_log = 0;

Ihist = [1:min(ML,10),ML+1:ML+min(ML,10),2*ML+1:2*ML+min(ML,10),3*ML+1:3*ML+min(ML,10)];
ahist = a(Ihist);
figure

while (n<n_total)
    
    n = n+1;
    nwrite = mod(n-1,dn_ROM)+1;
    nread = mod(n,dn_ROM)+1;
    
    p = Cal_p(Phi*a+Ustd);
    phist(:,nwrite) = p;
    pMonhist(n) = p(nMon);
    
    switch ROM_METHOD
        case 1
            a = RK4(a,PhiT,dt_ROM);
        case 2
            % Iter_max = num of subiterations
            Iter_max = 4;
            a = GN(a,PhiT,dt_ROM,ML,Iter_max,scheme);
    end
    
    ahist = [ahist a(Ihist)];
    % skip logging some iterations
    if (mod(n,20) ~= 0),	continue;    end
    n_log = n_log+1
    plot(t(1:20:n-1),pMonhist(1:20:n-1) - pstd(nMon));xlabel('time (s)');ylabel('p''(Pa)');
    title(['t_{FOM}=',num2str(FOM_span),', ML=',num2str(ML),method{ROM_METHOD},',CFL=',num2str(CFL_ROM)]);drawnow;
    print(gcf,'-djpeg',sprintf('-r%d',300),namejpg);
    %
    %     saveas(gcf,namefig);
    %     save(namemat,'n','a','ahist','phist','pMonhist','dt_ROM','n_FOM');
end
end
function a = RK4(a,PhiT,dt_ROM)
aold = a;
for subiter = 1:4
    R = Cal_Res_a(a);
    a = aold - dt_ROM/(5-subiter)*PhiT*R;
end
end

function a = GN(a,PhiT,dt_ROM,ML,Iter_max,scheme)
aold = a;
for subiter = 1:Iter_max
    R = Cal_Res_a(a);
    J = jacoRQ(a,scheme);
    Jr = eye(4*ML)/(dt_ROM/(Iter_max + 1 - subiter)) + PhiT*J;
    r = (a-aold)/(dt_ROM/(Iter_max + 1 - subiter)) + PhiT*R;
    da = -Jr\r;
    a = a + da;
end
end

