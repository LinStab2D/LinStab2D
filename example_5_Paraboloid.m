set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex')
clear variables
close all
clc
addpath(genpath('aux_matlab'));
addpath(genpath('BaseFlows'));

% Physical parameters
baseFlow.Re     = 5e3;          % Reynolds number
baseFlow.Ma     = 0.01;          % Mach number
baseFlow.Pr     = 0.7;          % Prandtl number
baseFlow.kappa 	= 1.4;          % heat capacity ratio
baseFlow.T_0 	= 293.15;       % temperature

% Perturbation & EVP parameters
freq            = 0;            % frequency
omega           = freq*2*pi;    % angular frequency
m               = 5;           % azimuthal wave number
nEig            = 3;            % Arnoldi method number of eigenvalues

% Domain & grid
Nr              = 100;           % # of grid points (radial)
Nz              = 100;           % # of grid points (streamwise)
FDorder         = 4;            % finite difference order of accuracy

% Flags
verbose         = true;         % visualize grid, base flow and results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create mesh and obtain differentiation matrices                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cartesian mesh in computational domain
y_symmetry      = true;     % use symmetry on y coordinate around y=0 
                            % (for axysymmetric problems)
x_periodicity   = false;    % use periodic b.c. on x
alphaFilter     = 0    ;    % spatial filter coefficient
xrange          = [-1 0 ];  % domain range in x
yrange          = [ 0 1 ];  % domain range in y

cmesh           = CreateMesh(xrange,yrange,Nz,Nr,FDorder, ...     
                          y_symmetry,x_periodicity,alphaFilter); %construct mesh
                     
x               = cmesh.X;           % x,y: Cartesian grid coordinates
% x   = (x.^4).*sign(x)./((x.^4)+.3);
x               = 1./( 1./(x.^4).*sign(x) + 1./x); x(isnan(x))=0;
x               = x/max(abs(x(:)));
y               = cmesh.Y;
y               = 1-(1-y).^1.4;      % spread points more uniformly on the surface
y               = 1-(1-y).^1.4;      % spread points more uniformly on the surface

% Grid transformation to parabolic sement in physical domain
% (   Here an analytical example is used, but X and Y can also be obtained
%        numerically )
Lx              = 2.25;              % wall-normal thickness parameter
Ly              = 1.5;              % Lenght parameter

x               = x*Lx-0.5;
y               = y*Ly;

z               = (x+1i*y);
X               = -real(z.^2);      % X,Y: Parabolic grid coordinates
Y               = -imag(z.^2);

%
mesh            = DeformMesh(cmesh,X,Y);

% Plot meshes in comutational and physical domains
if verbose
    figure
    
    subplot(1,2,1)
    title('Computatinal domain')
    hold on
    plot(cmesh.X, cmesh.Y, '-k', cmesh.X', cmesh.Y', '-k','HandleVisibility','off');
    for bondaries= fields(mesh.idx)'
        ids = mesh.idx.(bondaries{1});
        plot(cmesh.X(ids), cmesh.Y(ids),'linewidth',2);
    end
    xlabel('$\xi$');
    ylabel('$\eta$');
    legend(fields(mesh.idx))
    axis equal tight
    
    subplot(1,2,2)
    title('Physical domain')
    hold on
    plot(mesh.X, mesh.Y, '-k', mesh.X', mesh.Y', '-k','HandleVisibility','off');
    for bondaries= fields(mesh.idx)'
        ids = mesh.idx.(bondaries{1});
        plot(mesh.X(ids), mesh.Y(ids),'linewidth',2);
    end
    legend(fields(mesh.idx),'Location', 'Best')
    xlabel('$x$');
    ylabel('$y$');
    axis equal tight
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Base flow import                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import base flow

% Custumize the next function to read your baseflow
baseFlow    = example_5_readbaseflow(mesh,baseFlow); 

% Viscosity (via Sutherland) and heat conductivity 
% (via constant Prandtl number)
[baseFlow]  = sutherland_air(baseFlow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sponge                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sponge defined in the origial, non-deformed, mesh.
xs          = -0.5-Lx+.1;
xs_trans    = Lx/10;
ys          = Ly-.1;
ys_trans    = 0.1;
spongeAmp   = 5;
mesh.sponge = spongeAmp/2*max(tanh( (y-ys)/ys_trans)+1,  ...
    tanh(-(x-xs)/xs_trans)+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize base flow                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    figure('name','Base Flow')
    vars = {baseFlow.RHO ,'$\rho$';     baseFlow.U,'$U$';
            baseFlow.V,'$V$';           baseFlow.W,'$W$';
            baseFlow.T,'$T$';           mesh.sponge,'sponge';
            baseFlow.MU,'$\mu$';        baseFlow.dmudT,'$\frac{d\mu}{dT}$';
            baseFlow.d2mudT2                                ,'$\frac{d^2\mu}{dT^2}$';
            reshape(mesh.Dy*(baseFlow.W(:)),Nr,Nz)          ,'$\frac{dW}{dy}$';
            reshape(mesh.Dx*(baseFlow.W(:)),Nr,Nz)         ,'$\frac{dW}{dx}$';
            reshape(mesh.Dx*(mesh.Dy*baseFlow.W(:)),Nr,Nz) ,'$\frac{d^2W}{dydx}$'};
    plotFlow(mesh.X,mesh.Y,vars,4,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up linear operator and boundary conditions                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dqdt = L0 q
%%
[L0,idx] 	= GetLinProblem(mesh,baseFlow,'axy',m);

% Enforce Dirichlet b,c on the top, right and left boundaries, for u,v,w
% and T
borders     = 'ltr';  vars = 'uvwT';
[L0,idx_dirchlet] ...
            = BC_Dirichlet(L0,idx,borders,vars);

[W,invW]    = GetCompEnergyNorm(mesh,baseFlow,'axy');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resolvent analysis                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dq/dt = L*q + B*v
% u = C*q
% v: input/forcing, u:output/response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input/forcing matrix
B                       = ones(mesh.ngp*5,1);
B([idx.T_j;idx.rho_j])  = 0; % no temperature and density forcing
B( idx_dirchlet  )      = 0; % no forcing on dirichlet b.c. 
B                       = spdiags(B,0,mesh.ngp*5,mesh.ngp*5);

% Output/observable matrix
C                       = B;    

% Compute optimal forcings and responses 
[gains,responses,forces]                 = resolvent(L0,omega,nEig,W,invW,B,C,mesh.filters);

%% Plot modes and gains
if verbose
    figure('name','Mode gains')
    bar(gains);
    xlabel('mode');
    ylabel('gain');
    title('Resolvent gains')
    
    figure('name','Resolvent forcing and response modes')
    vars = {real(forces(idx.u_j,1)) ,'$f_u^{(1)}$'; real(responses(idx.u_j,1)) ,'$u^{(1)}$';
            real(forces(idx.u_j,2)) ,'$f_u^{(2)}$'; real(responses(idx.u_j,2)) ,'$u^{(2)}$';
            real(forces(idx.u_j,3)) ,'$f_u^{(3)}$'; real(responses(idx.u_j,3)) ,'$u^{(3)}$';};
    plotFlow(mesh.X,mesh.Y,vars,3,2);

end
