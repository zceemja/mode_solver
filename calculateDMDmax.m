function [ofValue,DGDlp,test] = calculateDMDmax(fN,ch4L,fNmodes)
% Calculates max DMD across the C-band
% AUTHOR:  Filipe M. Ferreira (filipef@ieee.org)

%% Physical constans
e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = 1/sqrt(e0*u0);

%% Frequency Range
% wArb    = 2*pi*c0/1.565e-6:2*dw:2*pi*c0/1.530e-6+0*dw;
% fN.w0  	= mean(wArb);
lArb    = [1565 1547.5 1530]*1e-9;%2*pi*c0./(wArb);

% dl      = 35e-9/3.9869;
% dw      = 2*pi*c0/(1.55e-6)^2*dl;
% wArb    = 2*pi*c0/1.565e-6:2*dw:2*pi*c0/1.530e-6+0*dw;
% % fN.w    = w;
% fN.w0  	= mean(wArb);
% lArb    = 2*pi*c0./(wArb);

%% Computation of DGDflat
w       = fN.w;
w0      = fN.w0;
lambda0 = 2*pi*c0/w0;

w1  = ch4L.w1*1e-6;
w2  = ch4L.w2*1e-6;
w3  = ch4L.w3*1e-6;
sI1 = ch4L.sI1;
sI2 = ch4L.sI2;
sI3 = ch4L.sI3;
eI1 = ch4L.eI1;
a1  = ch4L.a1;
a2  = 1e16;
a3  = 1e16;
radius = ch4L.w1*1e-6;

tic
[FO,V,U,DMD,TNOM] = fo4l_LpSl_mex(fN.incMax,ch4L.w1*1e-6,ch4L.w2*1e-6,ch4L.w3*1e-6,...
                                    ch4L.sI1,ch4L.sI2,ch4L.sI3,ch4L.eI1,...
                                    ch4L.a1,ch4L.a2,ch4L.a3,2*pi*c0/lambda0,ch4L.w1*1e-6,lambda0,lArb(1),lArb(2),lArb(3));
toc

if TNOM ~= length(fNmodes)
    DGDlp = 1e6*ones(3,length(fNmodes));
    DGDflat =  max(max(abs(abs(DGDlp)*1e12-targetDMD)));
    ofValue =  DGDflat;
    %% Save
    saveComputedData(ofValue, tag1,['ofValue_',fileName1]);
    saveComputedData(DGDlp, tag1,['DGDlp_',fileName1]);
    
    test = 1;
    return
else
    k1 = 1;
    k2 = 1;
    DMD0=DMD(1:3,1:TNOM);
    newDMD = DMD(1:3,1:TNOM);
    TNOM0 = TNOM;
    while k1 <= TNOM0
        if V(k1) > 0
            %newcol
            newDMD0 = newDMD(1:3,1:k2);
            newDMD = [newDMD0,newDMD(1:3,k2:end)];
            DMD = newDMD;
            TNOM = TNOM + 1;
            k2 = k2 + 1;
        end
        k1 = k1 + 1;
        k2 = k2 + 1;
    end
    TNOM = length(DMD);
    DGDlp = DMD;
end

ofValue =  max(max(abs(DGDlp*1e12)));
