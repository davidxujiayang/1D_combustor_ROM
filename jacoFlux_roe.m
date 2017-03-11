function [FL,FR] = jacoFlux_roe(UL, UR, gammaL,gammaR,AM)

global  dFdQ dFidQL dFidQR
run 'constantsch4.m'
N = length(UL)/4;

% left flux
FL = zeros(4*N, 4*N); % for allocation
for i=1:N
    FL([i,i+N,i+2*N,i+3*N],[i,i+N,i+2*N,i+3*N]) = (dFdQ(UL(i),UL(i+N),UL(i+2*N),UL(i+3*N),gammaL(i))-dFidQL(UL(i),UL(i+N),UL(i+2*N),UL(i+3*N),UR(i),UR(i+N),UR(i+2*N),UR(i+3*N)))*AM(i);
end

% right flux
FR = zeros(4*N, 4*N); % for allocation
for i=1:N
    FR([i,i+N,i+2*N,i+3*N],[i,i+N,i+2*N,i+3*N]) = (dFdQ(UR(i),UR(i+N),UR(i+2*N),UR(i+3*N),gammaR(i))-dFidQR(UL(i),UL(i+N),UL(i+2*N),UL(i+3*N),UR(i),UR(i+N),UR(i+2*N),UR(i+3*N)))*AM(i);
end

FL = 0.5*FL;
FR = 0.5*FR;
end

