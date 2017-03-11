global fa_p fa_m fb_p fb_m dFdQ Mach
for i = 1:3
    fname = ['global faQ',num2str(i),'_p'];
    eval(fname)
    
    fname = ['global faQ',num2str(i),'_m'];
    eval(fname)
    
    fname = ['global fbQ',num2str(i),'_p'];
    eval(fname)
    
    fname = ['global fbQ',num2str(i),'_m'];
    eval(fname)
end
for i = 1:4
    fname = ['global fcQ',num2str(i),'_p'];
    eval(fname)
    
    fname = ['global fcQ',num2str(i),'_m'];
    eval(fname)
end

load dFdQ_vl;
fa_p = eval(fa_p);
fa_m = eval(fa_m);
fb_p = eval(fb_p);
fb_m = eval(fb_m);
dFdQ = eval(dFdQ);
Mach = eval(Mach);

for i = 1:3
    fname = ['faQ',num2str(i),'_p = eval(faQ',num2str(i),'_p);'];
    eval(fname)
    
    fname = ['faQ',num2str(i),'_m = eval(faQ',num2str(i),'_m);'];
    eval(fname)
    
    fname = ['fbQ',num2str(i),'_p = eval(fbQ',num2str(i),'_p);'];
    eval(fname)
    
    fname = ['fbQ',num2str(i),'_m = eval(fbQ',num2str(i),'_m);'];
    eval(fname)
end

for i = 1:4
    fname = ['fcQ',num2str(i),'_p = eval(fcQ',num2str(i),'_p);'];
    eval(fname)
    
    fname = ['fcQ',num2str(i),'_m = eval(fcQ',num2str(i),'_m);'];
    eval(fname)
end