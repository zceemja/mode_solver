function out = doubleIntegralMatrix(F,X,Y)
% 2D integral of function F in a region defined by X and Y.
%
% Input Parameters:
%   F   - function         
%   X	- radial distance coordinate (varies horizontally)
%   Y	- azimuth coordinate (varies vertically)


dx = abs(X(2,1)-X(1,1));
dy = abs(Y(1,2)-Y(1,1));
if sum(X(end,1:end) == 2*pi) == length(X(end,1:end))
    out = sum( sum(F(1:end-1,:)  , 2) )* dx* dy;
else
    out = sum( sum(F  , 2) )* dx* dy;
end