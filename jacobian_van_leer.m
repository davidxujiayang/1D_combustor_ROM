run constantsch4
syms Q1 Q2 Q3 Q4 gamma real;

r = Q1;
u = Q2/r;
Y = Q4/r;

p = (gamma - 1)*(Q3 - 0.5*Q2*u);
c = sqrt(gamma*p/r);
M = u/c;

fa_p = 0.25*r*c*(M+1)^2;
fa_m = -0.25*r*c*(M-1)^2;
fb_p = (gamma-1)*M*c+2*c;
fb_m = (gamma-1)*M*c-2*c;
fc_p = fa_p*Y;
fc_m = fa_m*Y;

for i = 1:4
    fname = ['faQ',num2str(i),'_p = func2str(matlabFunction(simplify(diff(fa_p,Q',num2str(i),'))));'];
    eval(fname)
    
    fname = ['faQ',num2str(i),'_m = func2str(matlabFunction(simplify(diff(fa_m,Q',num2str(i),'))));'];
    eval(fname)
    
    fname = ['fbQ',num2str(i),'_p = func2str(matlabFunction(simplify(diff(fb_p,Q',num2str(i),'))));'];
    eval(fname)
    
    fname = ['fbQ',num2str(i),'_m = func2str(matlabFunction(simplify(diff(fb_m,Q',num2str(i),'))));'];
    eval(fname)
    
    fname = ['fcQ',num2str(i),'_p = func2str(matlabFunction(simplify(diff(fc_p,Q',num2str(i),'))));'];
    eval(fname)
    
    fname = ['fcQ',num2str(i),'_m = func2str(matlabFunction(simplify(diff(fc_m,Q',num2str(i),'))));'];
    eval(fname)
end
fa_p = func2str(matlabFunction(simplify(fa_p)));
fa_m = func2str(matlabFunction(simplify(fa_m)));
fb_p = func2str(matlabFunction(simplify(fb_p)));
fb_m = func2str(matlabFunction(simplify(fb_m)));

save('./data/dJdQ.mat',...
    'fa_p','fa_m','fb_p','fb_m',...
    'faQ1_p','faQ2_p','faQ3_p',...
    'faQ1_m','faQ2_m','faQ3_m',...
    'fbQ1_p','fbQ2_p','fbQ3_p',...
    'fbQ1_m','fbQ2_m','fbQ3_m',...
    'fcQ1_p','fcQ2_p','fcQ3_p','fcQ4_p',...
    'fcQ1_m','fcQ2_m','fcQ3_m','fcQ4_m'...
    );
