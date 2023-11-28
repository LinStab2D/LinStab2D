set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex')
clear variables
close all
clc
addpath(genpath('aux_matlab'));
verbose         = true;

% Physical parameters
baseFlow.Re     = 1e4;          % Reynolds number
baseFlow.Ma     = 0.9;          % Mach number
baseFlow.Pr     = 0.7;          % Prandtl number
baseFlow.kappa 	= 1.4;          % heat capacity ratio
baseFlow.T_0    = 288.15;       % temperature

m               = 0;
St              = 0.39;
omega           = 2*pi*St;
nEig            = 3;            % Arnoldi method number of eigenvalues

% Load interpolated LES mean flow
load('BaseFlows/M09jet.mat')

% Domain & grid
Nr              = numel(r_1D);  % # of grid points (radial)
Nz              = numel(z_1D);  % # of grid points (streamwise)
FDorder         = 4;            % finite difference order of accuracy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create mesh and obtain differentiation matrices                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cartesian mesh in computational domain
y_symmetry      = true;     % use symmetry on y coordinate around y=0 
                            % (for axisymmetric problems)
x_periodicity   = false;    % use periodic b.c. on x
alphaFilter     = 0.0;      % spatial filter coefficient
xrange          = [0 1];    % domain range in x
yrange          = [0 1];    % domain range in y

cmesh           = CreateMesh(xrange,yrange,Nz,Nr,FDorder, ...     
                          y_symmetry,x_periodicity,alphaFilter); %construct mesh

% Cartesian mesh in physical domain                      
[X,Y]           = meshgrid(z_1D,r_1D);
mesh            = DeformMesh(cmesh,X,Y);

% Plot meshes in computational and physical domains
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
U = reshape(U,[numel(r_1D) numel(z_1D)]);
V = reshape(V,[numel(r_1D) numel(z_1D)]);
W = reshape(W,[numel(r_1D) numel(z_1D)]);
T = reshape(T,[numel(r_1D) numel(z_1D)]);
RHO = reshape(RHO,[numel(r_1D) numel(z_1D)]);
spongeFun = reshape(spongeFun,[numel(r_1D) numel(z_1D)]);


baseFlow.RHO    = RHO;    % baseflow density
baseFlow.U      = W;      % baseflow U velocity
baseFlow.V      = U;      % baseflow V velocity
baseFlow.W      = V;      % baseflow W velocity
baseFlow.T      = T;      % baseflow temperature

% Viscosity (via Sutherland) and heat conductivity
% (via constant Prandtl number)
[baseFlow]  = sutherland_air(baseFlow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sponge                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mesh.sponge = spongeFun;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize base flow                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    figure('name','Base Flow')
    vars = {baseFlow.RHO ,'$\rho$';     baseFlow.U,'$U$';
            baseFlow.V,'$V$';           baseFlow.W,'$W$';
            baseFlow.T,'$T$';           mesh.sponge,'sponge';
            baseFlow.MU,'$\mu$';        baseFlow.dmudT,'$\frac{d\mu}{dT}$';
            baseFlow.d2mudT2                               ,'$\frac{d^2\mu}{dT^2}$';
            reshape(mesh.Dy*(baseFlow.W(:)),Nr,Nz)         ,'$\frac{dW}{dy}$';
            reshape(mesh.Dx*(baseFlow.W(:)),Nr,Nz)         ,'$\frac{dW}{dx}$';
            reshape(mesh.Dx*(mesh.Dy*baseFlow.W(:)),Nr,Nz) ,'$\frac{d^2W}{dydx}$'};
    plotFlow(mesh.X,mesh.Y,vars,4,3,[],'linecolor','none');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up linear operator and boundary conditions                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dqdt = L0 q
%%
[L0,idx]            = GetLinProblem(mesh,baseFlow,'axy',m);

% Enforce Dirichlet b,c on the top, right and left boundaries, for u,v,w
% and T
borders             = 'ltr';  
vars                = 'uvwT';
[L0,idx_dirchlet]   = BC_Dirichlet(L0,idx,borders,vars);

[W,invW]            = GetCompEnergyNorm(mesh,baseFlow,'axy');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resolvent analysis                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dq/dt = L*q + B*v
% u = C*q
% v: input/forcing, u:output/response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input/forcing matrix
B                   = mesh.sponge==0;
B                   = kron(B(:),ones(5,1));

% Output/observable matrix
C                   = B;    

% Compute optimal forcings and responses
[S,U,V]             = resolvent(L0,omega,nEig,W,invW,B,C,mesh.filters);

%% Plot modes and gains
if verbose
    figure('name','Mode gains')
    bar(S);
    xlabel('mode');
    ylabel('gain');
    title('Resolvent gains')
    
    figure('name','Resolvent forcing and response modes with filter')
    vars = {real(V(idx.u_j,1)) ,'$f_u^{(1)}$'; real(U(idx.u_j,1)) ,'$u^{(1)}$';
            real(V(idx.u_j,2)) ,'$f_u^{(2)}$'; real(U(idx.u_j,2)) ,'$u^{(2)}$';
            real(V(idx.u_j,3)) ,'$f_u^{(3)}$'; real(U(idx.u_j,3)) ,'$u^{(3)}$';};
        
    plotFlow(mesh.X,mesh.Y,vars,3,2,[],'linecolor','none');
end
