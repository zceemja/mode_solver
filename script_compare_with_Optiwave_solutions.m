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

%% Fibre Parameters
ch4L.w1  = 4;     % width - core 
ch4L.w2  = 4;     % width - core to trench 
ch4L.w3  = 4;     % width - trench 
ch4L.w4  = 5;     % width - fourth layer (if needed) 
ch4L.sI1 = 0.3/100;  % relative refractive index - core
ch4L.sI2 = 0;        % relative refractive index - region in between core and trench
ch4L.sI3 = -0.7/100; % relative refractive index - trench
ch4L.sI4 = 0;        % relative refractive index - fourth layer (if needed)
ch4L.eI1 = 0;
ch4L.a1  = 1e16;     % grading exponent - core 1.9651
ch4L.a2  = 1e16;     % grading exponent - for region in between core and trench
ch4L.a3  = 1e16;     % grading exponent - trench

%% MODE SOLVER
lambda0     = 1.55e-6;
lambdaArb   = [1.5:0.01:1.6]*1e-6;
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
end

figure()
subplot(2,2,1)
plot(coreWiTrenchwithcorrectsellmeircoeffs2022.lambda,coreWiTrenchwithcorrectsellmeircoeffs2022.pc,'o-'); hold on
plot(coreWiTrenchwithcorrectsellmeircoeffs2022.lambda,pc,'+-'); 
ylabel('{\itn_{eff}}')
xlabel('lambda [nm]')
subplot(2,2,2)
plot(coreWiTrenchwithcorrectsellmeircoeffs2022.lambda,coreWiTrenchwithcorrectsellmeircoeffs2022.mbl,'o-'); hold on
plot(coreWiTrenchwithcorrectsellmeircoeffs2022.lambda,mbl_v1,'+-'); 
ylabel('{\itMBL} [dB/km]')
xlabel('lambda [nm]')
subplot(2,2,3)
plot(coreWiTrenchwithcorrectsellmeircoeffs2022.lambda,coreWiTrenchwithcorrectsellmeircoeffs2022.D,'o-'); hold on
plot(coreWiTrenchwithcorrectsellmeircoeffs2022.lambda,D*1e6,'+-'); 
ylabel('{\itD} [ps/nm/km]')
xlabel('lambda [nm]')
subplot(2,2,4)
plot(coreWiTrenchwithcorrectsellmeircoeffs2022.lambda,coreWiTrenchwithcorrectsellmeircoeffs2022.mfd,'o-'); hold on
plot(coreWiTrenchwithcorrectsellmeircoeffs2022.lambda,mfd_v*1e6,'+-'); 
ylabel('{\itMFD} [\mum]')
xlabel('lambda [nm]')




