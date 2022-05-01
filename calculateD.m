function [ng,vg,D,D1,S1] = calculateD(df,incMax,lambda0,lambdaT,ch4L)
% Calculates chrometic dispersion
% AUTHOR:  Filipe M. Ferreira (filipef@ieee.org)

%% Physical constans
e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = 1/sqrt(e0*u0);

%% Frequency Range
w0      = 2*pi*c0/lambda0;
% lambdaT = 2*pi*c0./wArb;

w       = 2*pi*c0./lambdaT;
dw      = 2*pi*df;
wArb    = w + [-4:4]*dw;
lArb    = 2*pi*c0./(wArb);

%% Modes
for k1 = 1:length(wArb)
    [V,U,PC,TNOM] = modeSolver4l_LpSl_Verbose_mex(incMax,ch4L.w1*1e-6,ch4L.w2*1e-6,ch4L.w3*1e-6,...
                                                    ch4L.sI1,ch4L.sI2,ch4L.sI3,ch4L.eI1,...
                                                    ch4L.a1,ch4L.a2,ch4L.a3,w0,ch4L.w1*1e-6,lArb(k1));
    modePc(1,k1)         = PC(1);
end

%% Dispersion Calculation
%% dneff_dw calculation
for l1 = 1:1%size(V.',1)
    for l2 = 3:length(modePc)-2
        der1neff_w(l1,l2-2) = ( 8*modePc(l1,l2+1) - 8*modePc(l1,l2-1) - 1*modePc(l1,l2+2) + modePc(l1,l2-2) ) / ( 12 * ( dw ) );
        der2neff_w(l1,l2-2) = (-1*modePc(l1,l2+2) +16*modePc(l1,l2+1) -30*modePc(l1,l2+0) +16*modePc(l1,l2-1) -modePc(l1,l2-2) ) / ( 12 * ( dw )^2 );
    end
end
der1neff0_w = der1neff_w(3);
der2neff0_w = der2neff_w(3);
der1neff0_w = (-modePc(1,9)+4*280/105*modePc(1,8)-280/5*modePc(1,7)+4*280/5*modePc(1,6)-4*280/5*modePc(1,4)+280/5*modePc(1,3)-4*280/105*modePc(1,2)+modePc(1,1))/(280*dw);
der2neff0_w = (-modePc(1,9)+8*560/315*modePc(1,8)-560/5*modePc(1,7)+8*560/5*modePc(1,6)-205*560/72*modePc(1,5)+8*560/5*modePc(1,4)-560/5*modePc(1,3)+8*560/315*modePc(1,2)-modePc(1,1))/(560*(dw)^2);
der3neff0_w = (7*modePc(1,9)-72*modePc(1,8)+169*2*modePc(1,7) -61*8*modePc(1,6)+61*8*modePc(1,4)-169*2*modePc(1,3)+72*modePc(1,2)- 7*modePc(1,1))/(240*(dw )^3);

D1 = -(2*pi/lArb(5)^2)*(2*der1neff0_w+wArb(5)*der2neff0_w)*1e6;
S1 = ((2*pi/lArb(5)^2)^2*c0*(3*der2neff0_w + wArb(5)*der3neff0_w) + (4*pi/lArb(5)^3)*(2*der1neff0_w+wArb(5)*der2neff0_w))/1000;
    
D = -2*pi/lArb(5)^2 * ( 2*der1neff0_w + wArb(5)*der2neff0_w );

%% Group Velocity
ng = modePc(5) + wArb(5)*der1neff0_w;
vg = c0 / ng;
