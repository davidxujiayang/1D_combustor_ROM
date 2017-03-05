function [F, smax] = FluxFunction(UL, UR,  PropL,PropR)
%
% OUTPUTS:
%  F   : the flux out of the left cell (into the right cell)
%  smax: the maximum propagation speed of disturbance
%

run 'constantsch4.m'

% process left state
rL = UL(:,1);
uL = UL(:,2)./rL;
pL = PropL(:,4);
YL = UL(:,4)./rL;
rHL = UL(:,3) + pL;
HL = rHL./rL;

% left flux
FL = UL; % for allocation
FL(:,1) = rL.*uL;
FL(:,2) = UL(:,2).*uL + pL;
FL(:,3) = rHL.*uL;
FL(:,4) = UL(:,4).*uL; % left fuel mass flux


% process right state
rR = UR(:,1);
uR = UR(:,2)./rR;
pR = PropR(:,4);
YR = UR(:,4)./rR;
rHR = UR(:,3) + pR;
HR = rHR./rR;

% right flux
FR = UR; % for allocation
FR(:,1) = rR.*uR;
FR(:,2) = UR(:,2).*uR + pR;
FR(:,3) = rHR.*uR;
FR(:,4) = UR(:,4).*uR; % right fuel mass flux


% difference in states
du = UR - UL;

% Roe average
di     = sqrt(rR./rL);
d1     = 1.0./(1.0+di);

ui     = (di.*uR + uL).*d1;
Hi     = (di.*HR + HL).*d1;
Yi     = (di.*YR + YL).*d1;
Mi     = M_mix(M_ox,M_pr,Yi);
gammai = gamma_mix(Mi,Yi);
gmi = gammai-1.0;

af     = 0.5.*ui.*ui;
ucp    = ui;
c2     = gmi.*(Hi - af);           % speed of sound squared
if (min(min([pL,rL,pR,rR,c2]))<=0),	error 'Non-physical state!';	end
ci     = sqrt(c2);
ci1    = 1.0./ci;

% eigenvalues
l = zeros(length(ucp),3);
l(:,1) = ucp+ci;
l(:,2) = ucp-ci;
l(:,3) = ucp;

% entropy fix
epsilon = ci.*.02;
for i=1:3,
    Il=find((l(:,i)<epsilon) & (l(:,i)>-epsilon));
    if ~isempty(Il)
        l(Il,i) = 0.5.*(epsilon(Il) + l(Il,i).*l(Il,i)./epsilon(Il));
    end
end

l = abs(l); l3 = l(3);

% average and half-difference of 1st and 2nd eigs
s1    = 0.5.*(l(:,1) + l(:,2));
s2    = 0.5.*(l(:,1) - l(:,2));

% left eigenvector product generators
G1    = gmi.*(af.*du(:,1) - ui.*du(:,2) + du(:,3));
G2    = -ucp.*du(:,1)+du(:,2);

% required functions of G1 and G2
C1    = G1.*(s1-l3).*ci1.*ci1 + G2.*s2.*ci1;
C2    = G1.*s2.*ci1          + G2.*(s1-l3);

% flux assembly
F = FL; % for allocation

F(:,1)    = 0.5.*(FL(:,1)+FR(:,1))-0.5.*(l3.*du(:,1) + C1   );
F(:,2)    = 0.5.*(FL(:,2)+FR(:,2))-0.5.*(l3.*du(:,2) + C1.*ui + C2);
F(:,3)    = 0.5.*(FL(:,3)+FR(:,3))-0.5.*(l3.*du(:,3) + C1.*Hi + C2.*ucp  );
F(:,4)    = 0.5.*(FL(:,4)+FR(:,4))-0.5.*(l3.*du(:,4) + C1.*Yi   );

% max wave speed
smax = max(l,[],2);
