function [M,F] = CreateFilter(X,Y,alpha_f)
%CREATEFILTER Summary of this function goes here
%   Detailed explanation goes here

yList = unique(Y);

n=numel(X);
F=sparse(n,n);
M=sparse(n,n);

for iy=1:length(yList)
    y           = yList(iy);
    pos         = find(Y==y);
    
    x           = X(pos);
    [x,order]   = sort(x);
    pos         = pos(order);
    
    %%check that all points are equaly spaced
    nx          = length(x);
    if nx>10
        [M_1D,F_1D] = filter_O10_diffMat(nx,alpha_f);
        F(pos,pos)  = F_1D;
        M(pos,pos)  = M_1D;
    else
        warning('CreateFilter: There are regions with less than 10 points. No filter will be applied in those.')
        F(pos,pos) = eye(nx);
        M(pos,pos) = eye(nx);
    end
end

