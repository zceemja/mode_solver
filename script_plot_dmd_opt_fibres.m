% Plot DMD for optimum fibre profiles
% AUTHOR:  Filipe M. Ferreira (filipef@ieee.org)

clear all
close all

%% Physical constans
e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = 1/sqrt(e0*u0);

%% Source parameters
lambda0  = 1550e-9;
fN.w0       = 2*pi*c0/lambda0;
lambda   = 1550e-9;
fN.w        = 2*pi*c0/lambda;

%% Optimum parameters
load bOptNew;   load cOptNew; load aOptNew; load dOptNew
load dn1optNew; load w1optNew

%% Calculate DMD - example
for KK1 = 1:size(aOpt,1)
    for KK0 = 1:size(aOpt,2)
        r1     = w1opt(KK1,KK0);
        delta1 = dn1opt(KK1,KK0);
        rr1(KK1,KK0) = r1;
        
        %% Characteristics of each of the 4 layers of the refractive index profile
        fN.radius1     = r1;     fN.radius2     = r1;
        fN.deltaN      = delta1;
        
        ch4L.sI1 = delta1;         ch4L.eI1 = 0;             ch4L.a1  = aOpt(KK1,KK0);  ch4L.w1  = w1opt(KK1,KK0);
        ch4L.sI2 = 0;              ch4L.eI2 = 0;             ch4L.a2  = 1e16;           ch4L.w2  = bOpt(KK1,KK0);
        ch4L.sI3 = dOpt(KK1,KK0);  ch4L.eI3 = dOpt(KK1,KK0); ch4L.a3  = 1e16;           ch4L.w3  = cOpt(KK1,KK0);
        ch4L.sI4 = 0;              ch4L.eI4 = 0;             ch4L.a4  = 1e16;           ch4L.w4  = Inf;
        ch4R = ch4L;
        
        %% Calculate fibre modes
        fN.incMax = 0.01; % maximum step, normalised to core radius
        fNmodes = FiberAnalyzerLp(fN,ch4L);
        
%         %% Calculate Chromatic Dispersion
%         lambdaArb = 2*pi*c0/fN.w0; % test frequency
%         df = 10e9;    % df for beta derivative calculation
%         [ng,vg,Dtemp,D1,S1] = calculateD(df,fN.incMax,lambda0,lambdaArb,ch4L);
%         D(KK1,KK0) = Dtemp;
        
        %% Calculate DMD
        [dmdOpt(KK1,KK0)] = calculateDMDmax(fN,ch4L,fNmodes);
    end
end

%% Plot results
figure()
plot(unique(dn1opt)*1000,dmdOpt*1000,'+-',unique(dn1opt),12+0*dmdOpt*1000,'k--'); grid; axis([1 5 0 15])
xlabel('\Delta{\itn_{core}} [\cdot10^3]')
ylabel('{\itDMD} [ps/km]')
leg = legend(num2str([3 6 10 15 21].'));
title(leg,'#modes')