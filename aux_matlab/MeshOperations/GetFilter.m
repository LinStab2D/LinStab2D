function [FILTER,FILTER_ct] = GetFilter(mesh,alpha_filter)
%GETFILTER Summary of this function goes here
%   Detailed explanation goes here

if strcmp(alpha_filter,'none')
    FILTER    = @(x) x ;
    FILTER_ct = @(x) x ;
else
    [M_x,F_x]=CreateFilter(mesh.X(mesh.usedInd),mesh.Y(mesh.usedInd),alpha_filter);
    [M_z,F_z]=CreateFilter(mesh.Y(mesh.usedInd),mesh.X(mesh.usedInd),alpha_filter);


    FILTER    = @(x) M_z \ (M_x \(F_x *(F_z* x)));
    FILTER_ct = @(x) F_z'* (F_x'*(M_x'\(M_z'\x)));
end
end

