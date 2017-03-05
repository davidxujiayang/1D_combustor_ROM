load('U_steady.mat');
load('Std.mat');
run constantsch4
Li = Std.Li; A = Std.A; xelem = Std.xelem; Aelem = Std.Aelem; V = Std.V; s = Std.s; Int = Std.Int;
Uinf = Std.Uinf; Propinf = Std.Propinf; dtstd = Std.dt_min; CFLstd = Std.CFL; pstd = Std.p;
mstd = U(1,2)*A(1);
nelem = length(xelem);
nIE   = nelem - 1;
lMon = 0.5080-Li;           % lMon = monitor position
nMon  = find(xelem>lMon,1,'first');