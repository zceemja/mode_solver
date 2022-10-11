function [n] = sellmeierEquation(w,wl,xGe,xF)
% Calculates the Sellmeir equation.
%
% Input Parameters:
%   fN      - struct containing the fiber physical parameters
%   wArb    - angular frquency for where the group index calculation
%   wRef    - angular frquency of the input refractive index in fN
%
% Reference:
%   Agrawal nonlinear (1.2.6)
%
%   Created: 09-03-2011 (Author: Filipe Ferreira)
%   Updated: 09-03-2011 - Filipe Ferreira

%% Physical constans
e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = 1/sqrt(e0*u0);

%% The Sellmeir parameters for Pure Silica
Bi0   = [ 0.697668  0.407339  0.889883];
li0   = [ 0.070861  0.113600  9.784231];

%% The Sellmeir parameters for Germanium dopants
dBiGe   = [ 0.031510  0.267302  -0.012950];
dliGe   = [ 0.001677  0.032138   0.318034];

%% The Sellmeir parameters for Germanium dopants
dBiF   = [ -3.234366  0.164911  1.369490];
dliF   = [ -1.108703  0.752919  2.906858];

%% Sellmeier Equation Parameters
Bi = Bi0 + xGe*dBiGe + xF*dBiF;
li = li0 + xGe*dliGe + xF*dliF;
wi = 2*pi*c0 ./ (li*1e-6);

%% Equation evaluation
er = ones(1,length(w));
for k1 = 1:3
    er = er + Bi(k1) * wi(k1).^2 ./ (wi(k1).^2 - w.^2);
end

%% Output
n = sqrt(er);

for k1=1:wl
    er(k1) = 1;
    for k2=1:3
        er(k1) = er(k1) + Bi(k2)*wi(k2)*wi(k2)/(wi(k2)*wi(k2) - w(k1)*w(k1));
    end
    n(k1) = sqrt(er(k1));
end
    
