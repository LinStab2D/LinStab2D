function [x] = segmspace(x_segm,Nx_segm,Nx_smooth)
% Smoothed piecewise linspace
%
% Inputs
% ------
% x_segm:   n nodes of th piecewise grid
% Nx_segm:  n-1 number of grid points for each segment
% Nx_smooth: number of times mvg is applied to smooth distribution
%
% OTS, 2017
Nx              = sum(Nx_segm);
N_segm          = length(x_segm)-1;
dz_segm         = diff(x_segm)./Nx_segm;
Nz_cumsum       = cumsum(Nx_segm);
x               = zeros(Nx,1);

x(1:Nx_segm(1))  = linspace(x_segm(1),x_segm(2)-dz_segm(1)/2,Nx_segm(1));
for si = 2:N_segm-1
    x(Nz_cumsum(si-1)+1:Nz_cumsum(si-1)+Nx_segm(si))  = linspace(x_segm(si)+dz_segm(si)/2,x_segm(si+1)-dz_segm(si)/2,Nx_segm(si));
end
si = N_segm;
x(Nz_cumsum(si-1)+1:Nz_cumsum(si-1)+Nx_segm(si))  = linspace(x_segm(si)+dz_segm(si)/2,x_segm(si+1),Nx_segm(si));

% smooth grid in z (only in the interior domain)
for i=1:Nx_smooth, x(2:end-1) = mvg(x(2:end-1)); end

x = x.';

end

