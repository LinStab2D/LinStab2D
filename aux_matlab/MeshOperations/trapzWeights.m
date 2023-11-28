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

function [weightX] = trapzWeights(x)

noX = length(x);
weightX = zeros(noX,1);
weightX(1) = (x(2)-x(1))/2;
for i=2:noX-1
    weightX(i) = (x(i)-x(i-1))/2 + (x(i+1)-x(i))/2;
end
weightX(noX) = (x(noX)-x(noX-1))/2;