function [x,A] = read_shape(N)

% Creates the mesh

in2m = 0.0254;     % inches to meters factor

Lc = 15 * in2m;        % Chamber length 
Li = 0.1397;           % Injector Length
Ln = 1.34  * in2m;      % Nozzle Length

Lb1 = 0.25  * in2m;     
Lb2 = 0.5  * in2m;     % Lengths of back-steps

Ri = 0.403 * in2m;           % Injector radius
Rc = 0.885 * in2m;           % Chamber radius
Rn = [0.4095 0.769] * in2m;  % Nozzle radii

dx = (Li + Lc + Lb2 + Ln)/(N-1);
Ni = ceil(Li/dx)+1;
Nb1 = ceil(Lb1/dx);
Nc = ceil((Lc-Lb1)/dx);
Nb2 = ceil(Lb2/dx);
Nn = N - (Ni+Nb1+Nc+Nb2);

x = [linspace(-Li,0,Ni), linspace(0+dx,Lb1,Nb1), linspace(Lb1+dx,Lc,Nc), linspace(Lc+dx,Lc+Lb2,Nb2), linspace(Lc+Lb2+dx,Lc+Lb2+Ln,Nn)]';

f = @(y) ...
    (y <=0)                         .*  Ri + ...
    (y > 0).*(y <= Lb1)             .* (Ri + (Rc - Ri)*sin(y/Lb1 * pi/2)) + ...
    (y > Lb1).*(y <= Lc)            .*  Rc + ... 
    (y > Lc).*(y <= Lc + Lb2)       .* (Rc - (Rc - Rn(1))*sin((y - Lc)/Lb2 * pi/2)) + ...
    (y > Lc + Lb2)                  .* (Rn(1) + (Rn(2) - Rn(1))/Ln * (y - (Lc + Lb2)));

R = f(x); 
A = pi*R.^2;         % cross-section area

end

% dx = (Li+Lc+Ln+lb2)/(N-1);
% ratio=1.2;
% dxm = dx/ratio;
% 
% x = ([-Li       :dx    :0,...
%     0 + dxm       :dxm    : lb1,...
%     lb1 + dx      :dx     : Lc,...
%     Lc + dxm  :dxm    : Lc + lb2,...
%     Lc + lb2 + dxm  :dxm    : (Lc + lb2 + Ln)]).';