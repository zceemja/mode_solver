% Test mode solver against optiwave for:
%   mode effective index
%   chromatic dispersion
%   macrobend loss
%   mode field diameter
% AUTHOR:  Filipe M. Ferreira (filipef@ieee.org)

clear
clc
close all

load coreWiTrenchwithcorrectsellmeircoeffs2022

%% Physical constans
e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = 1/sqrt(e0*u0);

% Third-order nonlinear optical susceptibilitie
chi3_1064 = 2.0e-22;  % m^2 / V^2
chi3_1907 = 1.6e-22;  % m^2 / V^2

%% Fibre Parameters
ch4L.w1  = 4.1;     % width - core
ch4L.w2  = 0;     % width - core to trench
ch4L.w3  = 0;     % width - trench
ch4L.w4  = 0;     % width - fourth layer (if needed)
ch4L.sI1 = 0.36/100;  % relative refractive index - core
ch4L.sI2 = 0;        % relative refractive index - region in between core and trench
ch4L.sI3 = 0/100; % relative refractive index - trench
ch4L.sI4 = 0;        % relative refractive index - fourth layer (if needed)
ch4L.eI1 = 0;
ch4L.a1  = 1e16;     % grading exponent - core 1.9651
ch4L.a2  = 1e16;     % grading exponent - for region in between core and trench
ch4L.a3  = 1e16;     % grading exponent - trench

%% MODE SOLVER
lambda0     = 1.55e-6;
lambdaArb   = [1.25:0.005:1.7]*1e-6;
fN.w0       = 2*pi*c0/lambda0;

Rc = [38]*1e-3; % calculate MLB @ Rc

fN.incMax = 0.01;
tag = 'c';

for k1 = 1:length(lambdaArb)
    tic
    [V,U,PC,TNOM] = modeSolver4l_LpSl_Verbose_mex(  fN.incMax,ch4L.w1*1e-6,ch4L.w2*1e-6,ch4L.w3*1e-6,...
        ch4L.sI1,ch4L.sI2,ch4L.sI3,ch4L.eI1,...
        ch4L.a1,ch4L.a2,ch4L.a3,2*pi*c0/lambda0,ch4L.w1*1e-6,lambdaArb(k1));
    toc
    pc(k1) = PC(1);

    %% Calculate Chromatic Dispersion
    wArb    = 2*pi*c0/lambdaArb(k1); % test frequency
    df      = 10e9;    % df for beta derivative calculation
    [ng,vg,Dtemp,D1,S1] = calculateD(df,fN.incMax,lambda0,lambdaArb(k1),ch4L);
    D(k1)   = Dtemp;

    %% FIELD DISTRIBUTION CALCULATION
    tic
    [ExAux,rhoAux,nRhoAux,rhoLAux] = field4l_LpSl_Verbose_mex_ext(fN.incMax,ch4L.w1*1e-6,ch4L.w2*1e-6,ch4L.w3*1e-6,...
        ch4L.sI1,ch4L.sI2,ch4L.sI3,ch4L.eI1,...
        ch4L.a1,ch4L.a2,ch4L.a3,2*pi*c0/lambda0,ch4L.w1*1e-6,lambdaArb(k1),V(1),U(1),PC(1));
    toc

    %figure(123);plot(rhoAux(1:rhoLAux-1),ExAux(1:rhoLAux-1)); hold on

    clear rhoL Ex nRho rho
    Ex(:,1)     = ExAux(1:rhoLAux-1);
    nRho(:,1)   = nRhoAux(1:rhoLAux-1);
    rho(:,1)    = rhoAux(1:rhoLAux-1);

    Ex2D    = [Ex]*(cos(V(k1)*linspace(0,2*pi,length(rho)))); Ex2D = Ex2D.';
    x2D     = [rho]*cos(linspace(0,2*pi,length(rho))); x2D = x2D.';
    y2D     = [rho]*sin(linspace(0,2*pi,length(rho))); y2D = y2D.';
    Rho     = ones(length(rho),1)*rho.';
    Phi     = linspace(0,2*pi,length(rho)).'*ones(1,length(rho));

    % Calculate MBL
    [alpha1] = MBL(Rc,Rho,Phi,nRho,lambdaArb(k1),ch4L.w1*1e-6,Ex2D,PC(1),V(1));

    disp(['alpha_dB/km   = ',num2str(alpha1.alpha_dB_p_km)])
    mbl_v1(k1)  = alpha1.alpha_dB_p_km;

    %% Mode Field Diameter
    cr = ch4L.w1*1e-6;
    mfd_v(k1) = cr*2*sqrt((doubleIntegralMatrix(abs(Ex2D).^2.*Rho,Phi,Rho)^2)/(pi*doubleIntegralMatrix(abs(Ex2D).^4.*Rho,Phi,Rho)));
    disp(['Mode Field Diameter [um] = ',num2str(mfd_v(k1) *1e6)])

    %% Effective Mode Area
    cr = ch4L.w1*1e-6;
    Aeff(k1)  = cr^2*((doubleIntegralMatrix(abs(Ex2D).^2.*Rho,Phi,Rho)^2)/(doubleIntegralMatrix(abs(Ex2D).^4.*Rho,Phi,Rho)));
    disp(['Effective Mode Area [um^2] = ',num2str(Aeff(k1) *1e12)])

    Tf = 1415;

    Rpure = 4.1*1e-4*(Tf+273);

    nCl = nRho(end);

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

    AlphaRS(k1)  = A/B/((lambdaArb(k1)*1e6)^4);

    %% IR absorption + OH
    %% Fugihara, M.C. and Pinto, A.N. (2009), Attenuation fitting functions.
    % Microw. Opt. Technol. Lett., 51: 2294-2296. doi:10.1002/mop.24617

    AlphaR(k1)  = 8.5e-25/((lambdaArb(k1))^4);
    AlphaUV(k1) = 19.9628*exp(-5.642e6*lambdaArb(k1));

    AlphaIR(k1)  = 4.6e-14*exp(1.72e7*lambdaArb(k1));

    ATT1(k1) = 0.195/((1+((1/lambdaArb(k1)-1/1.395e-6)/8830)^2));
    ATT2(k1) = 0.8324/((1+((1/lambdaArb(k1)-1/1.374e-6)/7000)^2));
    ATT3(k1) = 0.7445*exp(-((1/lambdaArb(k1)-1/1.364e-6)/4400)^2);

    %% Nonlinear Kerr refractive-index coefficient
    %% U. Gubler and C. Bosshard (2000), "Optical third-harmonic generation of fused silica in gas atmosphere: Absolute value of the third-order nonlinear optical susceptibility"
    % in American Physical Society. doi:10.1103/PhysRevB.61.10702
    
    chi3= lambdaArb(k1) * -4.7450e-17 + 2.5049e-22;
    n2(k1) = 3/(4*e0*c0*PC(1)^2)*chi3;

end
Aeff(2)=NaN;  % FIXME: There is a bug with second value

figure()
subplot(2,3,1)
%plot(lambdaArb,coreWiTrenchwithcorrectsellmeircoeffs2022.pc,'o-'); hold on
plot(lambdaArb,pc,'+-');
ylabel('{\itn_{eff}}')
xlabel('lambda [nm]')
subplot(2,3,2)
%plot(lambdaArb,coreWiTrenchwithcorrectsellmeircoeffs2022.mbl,'o-'); hold on
gamma=(2*pi)./lambdaArb .* (n2./Aeff);
plot(lambdaArb,gamma*1e3,'-+');
ylabel('{\it\gamma} [1/W/km]')
xlabel('lambda [nm]')
subplot(2,3,3)
%plot(lambdaArb,coreWiTrenchwithcorrectsellmeircoeffs2022.D,'o-'); hold on
plot(lambdaArb,D*1e6,'+-');
ylabel('{\itD} [ps/nm/km]')
xlabel('lambda [nm]')
subplot(2,3,4)
%plot(lambdaArb,coreWiTrenchwithcorrectsellmeircoeffs2022.mfd,'o-'); hold on
plot(lambdaArb,Aeff*1e12,'+-');
ylabel('{\itA_{eff}} [\mum^2]')
xlabel('lambda [nm]')
subplot(2,3,5)
plot(lambdaArb,n2,'+-');
ylabel('Nonlinear Kerr Index {\itn_2} [m^2/W]')
xlabel('lambda [nm]')

subplot(2,3,6)
plot(lambdaArb,n2,'+-');
ylabel('Nonlinear Kerr Index {\itn_2} [m^2/W]')
xlabel('lambda [nm]')

subplot(2,3,6)
AlphaOH = ATT1+ATT2+ATT3;
alpha = AlphaIR+AlphaRS+AlphaUV+0.024;
plot(lambdaArb,alpha);
hold on
plot(lambdaArb,AlphaRS,'--');
plot(lambdaArb,AlphaIR,'--');
plot(lambdaArb,AlphaUV,'--');
xlabel('Wavelength [nm]')
ylabel('Attenuation [dB/km]')
%legend('Alpha','Rayleigh scattering','Infrared absorption','UV absorption')



