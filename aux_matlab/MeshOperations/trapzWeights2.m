% 3D Integration weight matrix for trapazoidal rule
%
%
%
% 	Written by Oliver T. Schmidt, 2012
%	University of Stuttgart
% 	E-mail: iagoschm@iag.uni-stuttgart.de
%
% 	2012/06 Oliver T. Schmidt           initial upload
%

function [weightXY] = trapzWeights2(x,y)

noX = length(x);
weightX = zeros(noX,1);
weightX(1) = (x(2)-x(1))/2;
for i=2:noX-1
    weightX(i) = (x(i)-x(i-1))/2 + (x(i+1)-x(i))/2;
end
weightX(noX) = (x(noX)-x(noX-1))/2;

noY = length(y);
weightY = zeros(noY,1);
weightY(1) = (y(2)-y(1))/2;
for i=2:noY-1
    weightY(i) = (y(i)-y(i-1))/2 + (y(i+1)-y(i))/2;
end
weightY(noY) = (y(noY)-y(noY-1))/2;

weightXY    = weightX*weightY';