% Credits to Xun Mu
%% Attenuation Coefficient

% Physical constans
e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = 1/sqrt(e0*u0);

incMax = 0.01;
% %% SMF Parameters
lambda  = 1.55e-6;
w0      = 2*pi*c0/lambda;
dl      = 0.0229498e-6;
nCl = sellmeierEquation(w0, 1,0, 0);
core = 4.1e-6;
delta = 0.36/100;

k=1;
for lambda = 1.2e-6:50e-9:1.7e-6
%% MODE SOLVER
nCo = nCl / sqrt(1-2*delta);
[V,U,PC,TNOM] = modeSolver4l_LpSl_Verbose_mex([core,20e-6],[delta,0,0],1.55e-6,lambda,incMax);
[~,~,DL]     = dispersion_nl_LpSl_IMP_mex( [core,20e-6],[delta,0,0],1.55e-6,lambda,1e-9,incMax);
Dispersion(k) = DL*1e6;
%% FIELD DISTRIBUTION CALCULATION

    for k1 = 1:TNOM
    if V(k1) == 0&& U(k1) == 1
  
    [ExAux,rhoAux,nRhoAux,rhoLAux] = field4l_LpSl_Verbose_mex([core,15e-6],[delta,0,0],1.55e-6,lambda,incMax,V(k1),U(k1),PC(k1));

         clear rhoL Ex nRho rho
        rhoL      = rhoLAux-1;
        Ex(:,1)   = ExAux(1:rhoL);
        nRho(:,1) = nRhoAux(1:rhoL);
        rho(:,1)  = rhoAux(1:rhoL);

        rho_new   = linspace(0,rho(end), rhoL);
        Ex_new    = interp1(rho,Ex,rho_new);
        nRho_new  = interp1(rho,nRho,rho_new);
        clear Ex nRho rho
        rho(:,1)  = rho_new;
        nRho(:,1) = nRho_new;
        Ex(:,1)   = Ex_new;
    Rho  = ones(length(rho),1)*rho.';
    end
    end

Tf = 1415;

Rpure = 4.1*1e-4*(Tf+273);

deltan  = (nRho.^2-nCl.^2)./(2*nRho.^2);
for j = 1:length(deltan)
if deltan(j) > 0
    deltaGe(j) = deltan(j);
    deltaF(j)  = 0;
else if deltan(j) < 0
       deltaF(j) = -deltan(j);
       deltaGe(j)  = 0; 
    else
        deltaGe(j) = 0;
        deltaF(j)  = 0;
    end
end
end

R = Rpure*(1+0.62.*deltaGe*100+0.6.*(deltaF*100).^2);
clear deltaGe deltaF
A  = sum(R'.*Ex.^2)*(Rho(1,2)-Rho(1,1));
B  = sum(Ex.^2)*(Rho(1,2)-Rho(1,1));

%% Rayleigh Scattering
%% K. Tsujikawa, K. Tajima, K. Shiraki and I. Sankawa, "Method for Predicting Rayleigh Scattering Loss of Silica-Based Optical Fibers," 
% in Journal of Lightwave Technology, vol. 25, no. 8, pp. 2122-2128, Aug. 2007, doi: 10.1109/JLT.2007.899789.

AlphaRS(k)  = A/B/((lambda*1e6)^4);

%% IR absorption + OH
%% Fugihara, M.C. and Pinto, A.N. (2009), Attenuation fitting functions. 
% Microw. Opt. Technol. Lett., 51: 2294-2296. doi:10.1002/mop.24617

AlphaR(k)  = 8.5e-25/((lambda)^4);
AlphaUV(k) = 19.9628*exp(-5.642e6*lambda);

AlphaIR(k)  = 4.6e-14*exp(1.72e7*lambda);

ATT1(k) = 0.195/((1+((1/lambda-1/1.395e-6)/8830)^2));
ATT2(k) = 0.8324/((1+((1/lambda-1/1.374e-6)/7000)^2));
ATT3(k) = 0.7445*exp(-((1/lambda-1/1.364e-6)/4400)^2);

k = k +1;
end

AlphaOH = ATT1+ATT2+ATT3;
alpha = AlphaIR+AlphaRS+AlphaUV+0.024;
lambda = 1200:50:1700;
figure(1)
plot(lambda,alpha,'LineWidth',2);
hold on
plot(lambda,AlphaRS,'--','LineWidth',2);
plot(lambda,AlphaIR,'--','LineWidth',2);
plot(lambda,AlphaUV,'--','LineWidth',2);
xlabel('Wavelength [nm]')
ylabel('Attenuation [dB/km]')
legend('Alpha','Rayleigh scattering','Infrared absorption','UV absorption')

figure(2)
plot(lambda,Dispersion,'LineWidth',2)
xlabel('Wavelength [nm]')
ylabel('Dispersion [ps/(nm*km)]')