function [Dx,D2x]=gridtrans1D(eta,Deta,D2eta,x)
% generalized 1D grid transformation for coordinate x=x(eta)
% returns transformed differentiation matrices
% NOTE: eta is not needed, just for completeness

x           = x(:);

dxdeta      = Deta*x;
d2xdeta2    = D2eta*x;
detadx      = 1./dxdeta;
d2etadx2    = 1./d2xdeta2;
Dx          = sparse(diag(detadx)  *Deta);
D2x         = sparse(-diag(detadx.^3)*diag(1./d2etadx2)*Deta + diag(detadx.^2)*D2eta);
