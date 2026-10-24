function [M,F] = CreateFilter(X,Y,alpha_f)
%CREATEFILTER Build independent filters on contiguous computational grid lines.
%   X and Y are the retained coordinates of a uniform mesh, before deformation.
%   Masked gaps must terminate both the explicit and implicit filter stencils.

X = X(:);
Y = Y(:);
yList = unique(Y);
xList = unique(X);
dx = Inf;
if numel(xList)>1
    dx = min(diff(xList));
end

n=numel(X);
F=sparse(n,n);
M=sparse(n,n);

for iy=1:length(yList)
    y           = yList(iy);
    pos         = find(Y==y);
    
    x           = X(pos);
    [x,order]   = sort(x);
    pos         = pos(order);
    
    % Use the mesh spacing, not the number of surviving points, to find
    % holes. The tolerance allows for roundoff in the uniform coordinates.
    domList = [0; find(diff(x)>1.5*dx); numel(x)];
    for idom=1:numel(domList)-1
        dom_pos = pos(domList(idom)+1:domList(idom+1));
        nx = numel(dom_pos);
        if nx>10
            [M_1D,F_1D] = filter_O10_diffMat(nx,alpha_f);
            F(dom_pos,dom_pos) = F_1D;
            M(dom_pos,dom_pos) = M_1D;
        else
            warning('CreateFilter: There are regions with fewer than 11 points. No filter will be applied in those.')
            F(dom_pos,dom_pos) = speye(nx);
            M(dom_pos,dom_pos) = speye(nx);
        end
    end
end
