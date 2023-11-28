function w = fdweights(xi,x,m)
%FDWEIGHTS Finite difference weights for arbitrary nodes.
% Input:
%   xi  evaluation point for the derivative (scalar)
%   x   nodes (vector)
%   m   order of derivative sought (integer scalar)
% Output:
%   w   weights for the approximation to the mth derivative (vector)
%
% Examples:
%   h = 0.4;  x = h*(0:5)';    % equally spaced nodes
%   w = fdweights(x(2),x,1)    % first derivative at second node
%   [ w'*exp(x)  exp(x(2)) ]   % applied to find d(exp(x))/dx|_{x=h}
%
%   x = [0 0.25 1 1.5 2.5]';   % unequally spaced nodes
%   w = fdweights(0.5,x,0)     % interpolation weights
%   [ w'*sin(x)  sin(0.5) ]    % applied to find sin(0.5)

% This is a compact implementation, NOT an efficient one!
% See Fornberg's "A Practical Guide to Pseudospectral Methods", Cambridge
%                   Univ. Press

% Copyright 2007, Toby Driscoll

p = length(x)-1;
w = zeros(size(x));
x = x-xi;              % translate evaluation point to 0
for k = 0:p
  w(k+1) = weight(x,m,p,k);
end


function c = weight(x,m,j,k)
% Implement recursion for the weights, assuming evaluation point is zero.
% Input:
%   x   nodes (vector, length n)
%   m   order of derivative sought (integer scalar)
%   j   stencil width (i.e., use first j+1 nodes only)
%   k   index of node for this weight (in 0:j)
% Output:
%   c   FD weight 

if (m<0) | (m>j)        % undefined coeffs 
  c = 0;    
elseif (m==0) & (j==0)  % base case of one-point interpolation
  c = 1;   
else                    % generic recursion 
  if k<j
    c = (x(j+1)*weight(x,m,j-1,k) - m*weight(x,m-1,j-1,k))/(x(j+1)-x(k+1));
  else
    beta = prod(x(j)-x(1:j-1)) / prod(x(j+1)-x(1:j));
    c = beta*(m*weight(x,m-1,j-1,j-1) - x(j)*weight(x,m,j-1,j-1));
  end
end