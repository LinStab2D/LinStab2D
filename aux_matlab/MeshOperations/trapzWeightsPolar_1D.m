% 2D Integration weight matrix for trapazoidal rule
%
%
%
% OTS 2017

function [weight_r] = trapzWeightsPolar_1D(r)

% rings
nothetar = length(r);
weight_r = zeros(nothetar,1);
% weight_r(1) = pi*((r(2))/2)^2;
weight_r(1) = pi*( r(1) + (r(2)-r(1))/2)^2; % polar singularity treatment removes node at r=0
for i=2:nothetar-1
    weight_r(i) = pi*( r(i) + (r(i+1)-r(i))/2 )^2 - pi*( r(i) - (r(i)-r(i-1))/2 )^2;
end
weight_r(nothetar) = pi*r(end)^2 - pi*( r(end) - (r(end)-r(end-1))/2 )^2;
% return row or column vector consistent with r
weight_r = reshape(weight_r,size(r));