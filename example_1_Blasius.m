set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex')
clear variables
close all
clc
addpath(genpath('aux_matlab'));
verbose         = true;

% Physical parameters
baseFlow.Re     = 6e5;          % Reynolds number
baseFlow.Ma     = 0.3;          % Mach number
baseFlow.Pr     = 0.7;          % Prandtl number
baseFlow.kappa 	= 1.4;          % heat capacity ratio
baseFlow.T_0 	= 273.15 ;      % temperature

delta           = 1.72/sqrt(baseFlow.Re);
omega           = 100*(baseFlow.Re/1e6);
nEig            = 3;            % Arnoldi method number of eigenvalues

%% Coordinate axes
nx              = 300;
ny              = 100;
dx              = 0.01;
x_min           = 5*dx;
x_max           = 1.25+dx;
y_min           = 0;
y_i             = 0.01;     
y_max           = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create mesh and obtain differentiation matrices                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTES: -why have CreateMesh and DeformMesh separately?
%        -put y_symmetry,x_periodicity,alpha_filter in options and default
%         to false?
%        -missing default for sponge function
%        -opts for eigs hardwired in resolvent function

% Cartesian mesh in computational domain
y_symmetry      = false;        % use symmetry on y coordinate around y=0 (for axysymmetric problems)
x_periodicity   = false;        % use periodic b.c. on x
alpha_filter    = 0;            % spatial filter coefficient
xrange          = [x_min x_max];% domain range in x
yrange          = [-1 1];       % domain range in y
FDorder         = 4;
cmesh           = CreateMesh(xrange,yrange,nx,ny,FDorder, ...     
                             y_symmetry,x_periodicity,alpha_filter); %construct mesh                     

% Mapping function as in Hanifi 1996 and mesh deformation
X               = cmesh.X;              % x,y: Cartesian grid coordinates
Y               = cmesh.Y;
a               = y_i*y_max/(y_max-2*y_i);
b               = 1+2*a/y_max;
Y               = a*(1+Y)./(b-Y);  
mesh            = DeformMesh(cmesh,[],Y,true);

%% Blasius Solution
baseFlow        = example_1_readbaseflow(X,Y,baseFlow);
baseFlow        = sutherland_air(baseFlow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize base flow                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    figure('name','Base Flow')
    vars = {baseFlow.RHO ,'$\rho$';
            baseFlow.U,'$U$';
            baseFlow.V,'$V$';
            baseFlow.W,'$W$';
            baseFlow.T,'$T$';           
            baseFlow.MU,'$\mu$';};
    plotFlow(mesh.X,mesh.Y,vars,3,2);
    drawnow
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up linear operator and boundary conditions                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dqdt = L0 q
%%
[L0,idx] 	= GetLinProblem(mesh,baseFlow,'2D',0);

% Enforce Dirichlet b.c.s for u,v,w and T on all boundaries
borders='lrbt';  vars = 'uvwT';
[L0,idx_dirchlet] = BC_Dirichlet(L0,idx,borders,vars);

[W,invW] = GetCompEnergyNorm(mesh,baseFlow,'2D');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resolvent analysis                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dq/dt = L*q + B*v
% u = C*q
% v: input/forcing, u:output/response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input/forcing matrix
B                       = ones(mesh.ngp*5,1);
B([idx.T_j;idx.rho_j;idx.w_j])  = 0; % no temperature and density forcing
B( idx_dirchlet  )              = 0; % no forcing on Dirichlet b.c. 
B                               = spdiags(B,0,mesh.ngp*5,mesh.ngp*5);

% Output/observable matrix
C                       = B;    

% Compute optimal forcings and responses with and without filters (for
% comparison)
[gains,responses,forces]   = resolvent(L0,omega,nEig,W,invW,B,C,mesh.filters);

%% Plot modes and gains
if verbose
    figure('name','Mode gains')
    bar(gains);
    xlabel('mode');
    ylabel('gain');
    title('Resolvent gains')
    
    figure('name','Resolvent forcing and response modes with filter')
    vars = {real(forces(idx.u_j,1)) ,'$f_u^{(1)}$'; real(responses(idx.u_j,1)) ,'$u^{(1)}$';
            real(forces(idx.u_j,2)) ,'$f_u^{(2)}$'; real(responses(idx.u_j,2)) ,'$u^{(2)}$';
            real(forces(idx.u_j,3)) ,'$f_u^{(3)}$'; real(responses(idx.u_j,3)) ,'$u^{(3)}$';};
    plotFlow(mesh.X,mesh.Y,vars,3,2);
end


