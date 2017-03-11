function [FL,FR] = jacoFlux_vl(UL, UR, gammaL,gammaR,AM)

global  fa_p fa_m fb_p fb_m...
    faQ1_p faQ2_p faQ3_p faQ1_m faQ2_m faQ3_m...
    fbQ1_p fbQ2_p fbQ3_p fbQ1_m fbQ2_m fbQ3_m ...
    fcQ1_p fcQ2_p fcQ3_p fcQ4_p fcQ1_m fcQ2_m fcQ3_m fcQ4_m ...
    dFdQ Mach

N = length(UL)/4;

% left flux
FL = zeros(4*N, 4*N); % for allocation
ML = zeros(1,N);
for i=1:N
    Q1L = UL(i);Q2L = UL(i+N);Q3L = UL(i+2*N);Q4L = UL(i+3*N);gL = gammaL(i);
    
    ML(i) = Mach(Q1L,Q2L,Q3L,gL);
    
    if ML(i)>=1
        FL([i,i+N,i+2*N,i+3*N],[i,i+N,i+2*N,i+3*N])=dFdQ(Q1L,Q2L,Q3L,Q4L,gL)*AM(i);
        continue
    end
    
    fa = fa_p(Q1L,Q2L,Q3L,gL);
    fb = fb_p(Q1L,Q2L,Q3L,gL);
    ja1 = faQ1_p(Q1L,Q2L,Q3L,gL);
    ja2 = faQ2_p(Q1L,Q2L,Q3L,gL);
    ja3 = faQ3_p(Q1L,Q2L,Q3L,gL);
    jb1 = fbQ1_p(Q1L,Q2L,Q3L,gL);
    jb2 = fbQ2_p(Q1L,Q2L,Q3L,gL);
    jb3 = fbQ3_p(Q1L,Q2L,Q3L,gL);
    
    FL(i,i) = ja1*AM(i);
    FL(i,i+N) = ja2*AM(i);
    FL(i,i+2*N) = ja3*AM(i);
    
    FL(i+N,i) = 1/gL*(fb*ja1+fa*jb1)*AM(i);
    FL(i+N,i+N) = 1/gL*(fb*ja2+fa*jb2)*AM(i);
    FL(i+N,i+2*N) = 1/gL*(fb*ja3+fa*jb3)*AM(i);
    
    FL(i+2*N,i) = 1/(gL^2-1)*(0.5*fb^2*ja1+fa*fb*jb1)*AM(i);
    FL(i+2*N,i+N) = 1/(gL^2-1)*(0.5*fb^2*ja2+fa*fb*jb2)*AM(i);
    FL(i+2*N,i+2*N) = 1/(gL^2-1)*(0.5*fb^2*ja3+fa*fb*jb3)*AM(i);
    
    FL(i+3*N,i) = fcQ1_p(Q1L,Q2L,Q3L,Q4L,gL)*AM(i);
    FL(i+3*N,i+N) = fcQ2_p(Q1L,Q2L,Q3L,Q4L,gL)*AM(i);
    FL(i+3*N,i+2*N) = fcQ3_p(Q1L,Q2L,Q3L,Q4L,gL)*AM(i);
    FL(i+3*N,i+3*N) = fcQ4_p(Q1L,Q2L,Q3L,gL)*AM(i);
end

% right flux
FR = zeros(4*N, 4*N); % for allocation
for i=1:N
    if ML(i)>=1
        continue
    end
    Q1R = UR(i);Q2R = UR(i+N);Q3R = UR(i+2*N);Q4R = UR(i+3*N);gR = gammaR(i);
    
    fa = fa_m(Q1R,Q2R,Q3R,gR);
    fb = fb_m(Q1R,Q2R,Q3R,gR);
    ja1 = faQ1_m(Q1R,Q2R,Q3R,gR);
    ja2 = faQ2_m(Q1R,Q2R,Q3R,gR);
    ja3 = faQ3_m(Q1R,Q2R,Q3R,gR);
    jb1 = fbQ1_m(Q1R,Q2R,Q3R,gR);
    jb2 = fbQ2_m(Q1R,Q2R,Q3R,gR);
    jb3 = fbQ3_m(Q1R,Q2R,Q3R,gR);
    
    FR(i,i) = ja1*AM(i);
    FR(i,i+N) = ja2*AM(i);
    FR(i,i+2*N) = ja3*AM(i);
    
    FR(i+N,i) = 1/gR*(fb*ja1+fa*jb1)*AM(i);
    FR(i+N,i+N) = 1/gR*(fb*ja2+fa*jb2)*AM(i);
    FR(i+N,i+2*N) = 1/gR*(fb*ja3+fa*jb3)*AM(i);
    
    FR(i+2*N,i) = 1/(gR^2-1)*(0.5*fb^2*ja1+fa*fb*jb1)*AM(i);
    FR(i+2*N,i+N) = 1/(gR^2-1)*(0.5*fb^2*ja2+fa*fb*jb2)*AM(i);
    FR(i+2*N,i+2*N) = 1/(gR^2-1)*(0.5*fb^2*ja3+fa*fb*jb3)*AM(i);
    
    FR(i+3*N,i) = fcQ1_m(Q1R,Q2R,Q3R,Q4R,gR)*AM(i);
    FR(i+3*N,i+N) = fcQ2_m(Q1R,Q2R,Q3R,Q4R,gR)*AM(i);
    FR(i+3*N,i+2*N) = fcQ3_m(Q1R,Q2R,Q3R,Q4R,gR)*AM(i);
    FR(i+3*N,i+3*N) = fcQ4_m(Q1R,Q2R,Q3R,gR)*AM(i);
end
end

