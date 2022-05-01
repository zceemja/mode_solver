function modesLP = FiberAnalyzerLp(fN,ch4L)
% Calculates the normalized propagation constants of the LP mode of
% step-index fibers
%
% Input Parameters:
%   fN           - struct containing the fiber physical parameters
%   askRec       - ask for recalculation if necessary
%
% Reference:
%   D. Marcuse, "Theory of Dielectric Optical Waveguides"
%
% Author: Filipe Ferreira

%% Physical constans
e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = 1/sqrt(e0*u0);

%% LP Mode Solver Sliced
w           = fN.w;
w0          = fN.w0;
lambda      = 2*pi*c0/w;
lambda0      = 2*pi*c0/w0;

tic
[V,U,PC,TNOM] = modeSolver4l_LpSl_Verbose_mex(fN.incMax,ch4L.w1*1e-6,ch4L.w2*1e-6,ch4L.w3*1e-6,...
                                              ch4L.sI1,ch4L.sI2,ch4L.sI3,ch4L.eI1,...
                                              ch4L.a1,ch4L.a2,ch4L.a3,2*pi*c0/lambda0,ch4L.w1*1e-6,lambda);
toc
% fprintf(['TNOM: ',num2str(TNOM),'\n'])
for k1 = 1:TNOM
    modesLP(k1).vIndex = V(k1);
    modesLP(k1).uIndex = U(k1);
    modesLP(k1).pc = PC(k1);
    disp(['LP',num2str(V(k1)),num2str(U(k1)),': ',num2str(PC(k1),12)]);
end
if TNOM == 0
    modesLP = -1;
end
