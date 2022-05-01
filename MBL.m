function [alpha1] = MBL(Rc,Rho,Phi,nRho,lambda,w1,Ex2D,pc,v)
% Reference to calculate bend loss
% J. Sakai and T. Kimura, "Bending loss of propagation modes in arbitraryindex
% profile optical fibers," Appl. Opt., vol. 17, no. 10, pp. 14991506,
% 1978.
% Implementation:  Filipe M. Ferreira (filipef@ieee.org)

%% MBL
nCo  = max(nRho(1));
cr   = w1;
k0   = 2*pi/lambda;
nCl  = nRho(end);

delta = (nCo^2 - nCl^2) / (2*nCl^2); % Relative index difference between nMin and nMax
V     = k0*cr*(nCo^2 - nCl^2)^0.5;    % Normalized frequency
W     = cr*k0*(pc^2  - nCl^2)^0.5;    %

nClad_mask_2D = zeros(size(Rho));       % Cladding mask for Pclad calculation
nClad_mask_2D(Rho > 1) = 1;


S     = 0.5*real(Ex2D.*conj(Ex2D));                           % Poyting vector
P     = doubleIntegralMatrix(S.*Rho,Phi,Rho);                 % Total mode power
Pclad = doubleIntegralMatrix(nClad_mask_2D.*S.*Rho,Phi,Rho);  % Cladding power

if v == 0
    s = 2;
else
    s = 1;
end

gamma = sqrt(pi) * (Pclad/P) / (2 * s * cr* ( besselk(v-1,W)*besselk(v+1,W) - besselk(v,W)^2 ) ) .* ...
    exp(-4*delta*W^3*Rc/(3*cr*V^2)) ./ ( W * ( W*Rc/cr + V^2./(2*delta*W) ).^0.5);


alpha1.alpha_dB_p_turn = (2*pi.*Rc)*10/log(10).*gamma;
alpha1.alpha_dB_p_m    = 10/log(10)*gamma;
alpha1.alpha_dB_p_km   = 10/log(10)*gamma*1000;
