clear all
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex')
close all
clc
addpath(genpath('aux_matlab'));
addpath(genpath('BaseFlows'));

% Physical parameters
baseFlow.Re     = 1500;     % Reynolds number
baseFlow.Ma     = 0.5;      % Mach number
baseFlow.Pr     = 0.7;      % Prandtl number
baseFlow.kappa 	= 1.4;      % heat capacity ratio
baseFlow.T_0 	= 300.;     % temperature

% Perturbation & EVP parameters
omega_target    = 1.5+1.i;      % value around which the eigenvalues will be seek
beta            = 0;            % azimuthal wave number
nEig            = 50;           % Arnoldi method number of eigenvalues

% Domain & grid
Nx_cavity       = 80;           % # of grid points in the cavity (stream wise)
Ny_cavity       = 50;           % # of grid points in the cavity (wall normal)
FDorder         = 4;            % finite difference order of accuracy

% Flags
verbose         = true;         % visualize grid, base flow and results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create mesh and obtain differentiation matrices                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cartesian mesh in computational domain
y_symmetry      = true;      % use symmetry on y coordinate around y=0 
                             % (for axysymmetric problems)
x_periodicity   = false;     % use periodic b.c. on x
alpha           = .0;        % spatial filter coefficient

x_cavity_start  = 0.34029;
L_cavity        =  1.9995;   
xrange          = x_cavity_start + [ -1 2]*L_cavity;  % domain range in x
yrange          = [ -1 1 ] ;  % domain range in y

Ny=2*Ny_cavity+1;
Nx=3*Nx_cavity+1;

cmesh           = CreateMesh(xrange,yrange,Nx,Ny,FDorder, ...     
                            y_symmetry,x_periodicity,alpha); %construct mesh
                     
% Mash points outside the domain
x               = cmesh.X;           
y               = cmesh.Y;
mask            = (y<0) & ((x<.34) | (x>2.34) ) ;  
mesh            = MeshMask(cmesh,mask);

%% Mesh Deformation
strech_n        = 3;                     % polinomial strech order
% ---- Deformation in Y                         ----
y               = mesh.Y;
% ---- Concentrate points on the shear layer    ----
a               = 0.5;                   % mesh derivative at the center
y               = (y.^3+a*y)/(1+a);      % map function that concentrates points near y=0
% ---- Stretch mesh on the top border           ----
strech_y        = 0;                     % move points with y>strech_y
strech_L        = 2;                     % Distance on which point will be streched

p               = y>strech_y;            % indices for the points which will be moved
y(p)=y(p)+ strech_L*(  (y(p)-strech_y)/(max(y(p))-strech_y) ).^strech_n;

mesh            = DeformMesh(mesh,[],y); % apply mesh deformation in y

    
%  --- Deformation in X ---
strech_Lx       =   1;                   % Distance on which point will be streched
x               = mesh.X;
% ---- Stretch mesh on the left border          ----
% strech_xl =-0.5;                       % move points with y>strech_y
% p = x<strech_xl;                       % indices for the points which will be moved
% x(p)=x(p) -  strech_Lx*(  (x(p)-strech_xl)/(min(x(p))-strech_xl) ).^strech_n;
% ---- Stretch mesh on the right border          ----
strech_xr       = 2.8;                   % move points with y>strech_y
p               = x>strech_xr;                 % indices for the points which will be moved
x(p)            = x(p) +  strech_Lx*((x(p)-strech_xr)/(max(x(p))-strech_xr)).^strech_n;

mesh            = DeformMesh(mesh,x,[]); % apply mesh deformation in x

% If, instead of the above, we used 
% mesh = DeformMesh(mesh,x,y);    
% the solver is ca be slower and require more memory, as a general the 
% deformation reduces the sparticity of the differentiation matrices.

clear x y 
%%
% Plot meshes in comutational and physical domains
if verbose
    figure
   
    title('Physical domain')
    hold on
    ind = mesh.usedInd;
    plot(mesh.X(ind), mesh.Y(ind), '.k','HandleVisibility','off');
    
    p=mesh.usedInd;
    for bondaries= fields(mesh.idx)'
        ids = mesh.idx.(bondaries{1});
        plot(mesh.X(p(ids)), mesh.Y(p(ids)),'o');
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
baseFlow    = example_4_readbaseflow(mesh,baseFlow); % Custumize the function to read your baseflow

% Viscosity (via Sutherland) and heat conductivity (via constant Prandtl number)
[baseFlow]  = sutherland_air(baseFlow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sponge                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sponge defined in the origial, non-deformed, mesh.
x           = mesh.X;
y           = mesh.Y;

xs_trans    = .5;  xs_offset =  .5;
ys_trans    = .5 ; ys_offset =  1;
spongeAmp   = 5;
ind = mesh.usedInd;

mesh.sponge = nan(size(mesh.X));
mesh.sponge(ind) = ...
        spongeAmp/2*max(tanh( ( y(ind)-max(y(:)) + ys_offset )/ys_trans)+1, ...
                        tanh(-( x(ind)-min(x(:)) - xs_offset)/xs_trans)+1 + ...
                        tanh( ( x(ind)-max(x(:)) + xs_offset)/xs_trans)+1  );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize base flow                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    figure('name','Base Flow')
    vars = {baseFlow.RHO ,'$\rho$';     baseFlow.U,'$U$';
            baseFlow.V,'$V$';           baseFlow.W,'$W$';
            baseFlow.T,'$T$';           mesh.sponge,'sponge';
            baseFlow.MU,'$\mu$';        mesh.W,'$W$'; 
         };
    plotFlow(mesh.X,mesh.Y,vars,3,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up linear operator incl. boundary conditions                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dqdt = L0 q
%%
[L0,idx] 	= GetLinProblem(mesh,baseFlow,'2D',beta);

% Enforce Dirichlet b,c on the top, right and left boundaries, for u,v,w
% and T
borders='lbtrm';  vars = 'uvwT';
[L0,idx_dirchlet] = BC_Dirichlet(L0,idx,borders,vars);

[W,invW] = GetCompEnergyNorm(mesh,baseFlow,'2D');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Global modes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Divide by -1i to convert eigenvalues from s to omega (=-1i s)
L0 = L0/-1i;

%setup eigs options
opts.tol    = 1e-8;
opts.disp   = 2;
opts.p      = 300;


% compute direct and adjoint modes
tic;
[V    ,omega    ] = eigs(L0  ,nEig,omega_target,opts);
omega=diag(omega);

[V_adj,omega_adj] = eigs(invW*L0'*W,nEig,conj(omega_target),opts);
omega_adj=diag(omega_adj);

time_eigsMatlab = toc;

%% Plot and compare eigenvalues. Show some eigen vectors
if verbose    
    [~,iplot    ] = sort(imag(omega),"descend");
    [~,iplot_adj] = sort(imag(conj(omega_adj)),"descend");

    nplot = 3;
    iplot     = iplot(1:nplot);
    iplot_adj = iplot_adj(1:nplot);

%     [~,iplot]=min(abs(lambda-(.9*alpha+.10i)));
    figure('name','Temporal Spectra');
        plot(real(omega),imag(omega),'ob',real(omega_adj),-imag(omega_adj),'xb')
        xlabel('$\omega_r$')
        ylabel('$\omega_i$')

%             clear leg
        leg = {'Direct Spectra','Conjugate of Adjoint Spectra'};
        hold on
        for iiplot=iplot'
            plot(real(omega(iiplot)),imag(omega(iiplot)),'^','MarkerSize',10);
            leg{end+1} = num2str(omega(iiplot));
        end
        grid on
        xlabel('$\omega_r$');
        ylabel('$\omega_i$');
        legend(leg)

    ndofs = size(L0,1);
    used_dofs = 1:ndofs;  used_dofs(idx_dirchlet)=[];
    U = zeros(ndofs,1);
    for iiplot=iplot'
        U=V(:,iiplot);
        figure('name',['Direct eigenmode omega = ' num2str(omega(iiplot),'%.3f')])
        vars = {real(U(idx.rho_j)) ,'$\rho$'; 
                real(U(idx.u_j  )) ,'$u$'; 
                real(U(idx.v_j  )) ,'$v$'; 
                real(U(idx.w_j  )) ,'$w$'; 
                real(U(idx.T_j  )) ,'$T$' };
        title(['$\omega = ' num2str(omega(iiplot)) '$, Direct']);
        axs = plotFlow(mesh.X,mesh.Y,vars,3,2,mesh.usedInd);
    end
    for iiplot=iplot_adj'
        U=V_adj(:,iiplot);
        figure('name',['Adjoint eigenmode omega^* = ' num2str(conj(omega_adj(iiplot)),'%.3f')])
        vars = {real(U(idx.rho_j)) ,'$\rho$'; 
                real(U(idx.u_j  )) ,'$u$'; 
                real(U(idx.v_j  )) ,'$v$'; 
                real(U(idx.w_j  )) ,'$w$'; 
                real(U(idx.T_j  )) ,'$T$' };
        title(['$\omega^* =' num2str(conj(omega_adj(iiplot))) '$, Adjoint']);
        axs = plotFlow(mesh.X,mesh.Y,vars,3,2,mesh.usedInd);  
    end
end