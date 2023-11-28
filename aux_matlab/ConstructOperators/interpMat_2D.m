function [W] = interpMat_2D(x,y,X,Y,method)
% returns interpolation matrix W*f(x,y) = f(X,Y)

%%
% method  = 'cubic';
% f       = peaks;
% x       = 1:49;
% [x,y]   = ndgrid(x,x);
% 
% X       = linspace(0,49,100);
% Y       = linspace(1,50,60);
% [X,Y]   = ndgrid(X,Y);
% F       = interpn(x,y,f,X,Y,method);

%%
% method  = 'linear';
% f       = peaks;
% x       = 1:49;
% [x,y]   = meshgrid(x,x);
% 
% X       = linspace(0,49,100);
% Y       = linspace(1,50,60);
% [X,Y]   = meshgrid(X,Y);
% F       = interp2(x,y,f,X,Y,method);

nx  = size(x,1);
ny  = size(y,2);
nxy = nx*ny;

nX  = size(X,1);
nY  = size(Y,2);
nXY = nX*nY;

W   = zeros(nXY,nxy);
for i=1:ny
    for j=1:nx
        e_ij            = zeros(nx,ny);
        e_ij(j,i)       = 1;
        W_ij            = interp2(x,y,e_ij,X,Y,method);
%         W_ij            = interpn(x,y,e_ij,X,Y,method,0); % for ndgrid
        W(:,(i-1)*nx+j) = W_ij(:);
    end
end
W   = sparse(W);

% figure
% subplot(1,3,1), spy(W)
% subplot(1,3,2), pcolor(X,Y,reshape(W*f(:),nX,nY)), axis equal tight
% subplot(1,3,3), pcolor(x,y,f), axis equal tight