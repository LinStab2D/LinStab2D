set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex')
clear 
close all
clc
addpath(genpath('aux_matlab'));
addpath(genpath('BaseFlows'));

% Physical parameters
baseFlow.Re     = 28000;          % Reynolds number
baseFlow.Ma     = 6;          % Mach number
baseFlow.Pr     = 0.7;          % Prandtl number
baseFlow.kappa  = 1.4;          % heat capacity ratio
baseFlow.T_0    = 273.15 ;       % temperature


% Domain lenghts
Lx = pi ;
Ly = 5  ;

% Perturbation & EVP parameters
freq        = 0;            % frequency
omega       = freq*2*pi;    % angular frequency for spatial stability
alpha       = .6283;        % stream-wize wavenumber for temporal stability
nEig        = 50;           % Arnoldi method number of eigenvalues

nPeriods    = 2  ;           % 0 for periodic solution, 
                            % 1/n for a solution periodic in n domains.
if nPeriods==0 ; floquetExp=0;
else ; floquetExp  = 2*pi/(Lx*nPeriods) ;          % Floquet exponent      
end

% Domain & grid
Nx          = 125*(Lx/pi);           % # of grid points (radial)
Ny          = 250;           % # of grid points (streamwise)
FDorder     = 4;            % finite difference order of accuracy


% Flags
verbose     = true;         % visualize grid, base flow and results
temporalAna = true;
spatialAna  = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create mesh and obtain differentiation matrices                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cartesian mesh in computational domain
y_symmetry      = false;         % use symmetry on y coordinate around y=0 (for axysymmetric problems)
x_periodicity   = true;          % use periodic b.c. on x
alpha_filter    = .0;            % spatial filter coefficient
xrange          = [0  Lx ];  % domain range in x
yrange          = [ -1   1 ];  % domain range in y

cmesh           = CreateMesh(xrange,yrange,Nx,Ny,FDorder, ...     
                             y_symmetry,x_periodicity,alpha_filter); %construct mesh
                     
X   = cmesh.X;           % x,y: Cartesian grid coordinates
Y   = cmesh.Y;

% Mapping function as in Hanifi 1996 and mesh deformation
yi = 2.0                ; ymax = Ly ; 
a  = yi*ymax/(ymax-2*yi); b    = 1+2*a/ymax;
Y=a*(1+Y)./(b-Y);  

mesh    = DeformMesh(cmesh,[],Y,true);

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

baseFlow    = example_3_readbaseflow(mesh,baseFlow); % Custumize the function to read your baseflow

% Viscosity (via Sutherland) and heat conductivity (via constant Prandtl number)
% [baseFlow]  = sutherland_air(baseFlow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sponge                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sponge defined in the origial, non-deformed, mesh.
ys_trans    = .5;
ys          = ymax-1;
spongeAmp   = .1;
mesh.sponge = spongeAmp/2*(tanh( (Y-ys)/ys_trans)+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize base flow                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    figure('name','Base Flow')
    vars = {baseFlow.RHO ,'$\rho$';     baseFlow.U,'$U$';
            baseFlow.V,'$V$';           baseFlow.W,'$W$';
            baseFlow.T,'$T$';           mesh.sponge,'sponge';
            baseFlow.MU,'$\mu$';        baseFlow.dMUdT,'$\frac{d\mu}{dT}$';
            baseFlow.d2MUdT2                                ,'$\frac{d^2\mu}{dT^2}$';
            reshape(mesh.Dy*(baseFlow.W(:)),Nx,Ny)          ,'$\frac{dW}{dy}$';
            reshape(mesh.Dx*(baseFlow.W(:)),Nx,Ny)         ,'$\frac{dW}{dx}$';
            reshape(mesh.Dx*(mesh.Dy*baseFlow.W(:)),Nx,Ny) ,'$\frac{d^2W}{dydx}$'};
    plotFlow(mesh.X,mesh.Y,vars,4,3);
    drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Temporal Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if temporalAna

    % Get linear opreator
    [L,idx]     = GetLinProblem(mesh,baseFlow,'2D',alpha,floquetExp);
    % Enforce Dirichlet b,c on the top, right and left boundaries, for u,v,w
    % and T
    borders='tb';  vars = 'uvwT';
    [L,idx_dirchlet] = BC_Dirichlet(L,idx,borders,vars);

    %%     Compute eigenvalues       %%
    omega_tar   = 0.5498+1.1i;   % value around which the eigenvalues will be seek

    %-------- Classical approach ------

    % Remove constrains lines/rows and get L1 : \omega \hat q = L1 \hat q
    L1 = L/-1i;
%     L1(idx_dirchlet,:)=0;
%     L1(:,idx_dirchlet)=0;
%     L1(idx_dirchlet,idx_dirchlet) = -1e4;

    %setup eigs options
    opts.tol    = 1e-8;
    opts.disp   = 2;
    opts.p      = 1500;


    % compute eigs
    tic;
    [V,omega] = eigs(L1,nEig,omega_tar,opts);
    omega=diag(omega);
    time_eigsMatlab = toc;

    %% Plot and compare eigenvalues. Show some eigen vectors
    if verbose    
        [~,iplot] = sort(imag(omega),"descend");
        iplot=iplot(1:4);

    %     [~,iplot]=min(abs(lambda-(.9*alpha+.10i)));
        figure('name','Temporal Spectra');
            plot(real(omega),imag(omega),'o')
            xlabel('$\omega_r$')
            ylabel('$\omega_i$')
            
            leg{1} = 'Spectra';
            hold on
            for iiplot=iplot'
                plot(real(omega(iiplot)),imag(omega(iiplot)),'o');
                leg{end+1} = num2str(omega(iiplot));
            end
            grid on
            xlabel('$\omega_r$');
            ylabel('$\omega_i$');
            legend(leg)

        ndofs = size(L,1);
        
        [nDoms,extX,extY,extV,ext_idx,ext_usedInd] = ...
            ReconstructFloquetSol(X,Y,V,idx,mesh.usedInd,floquetExp);
        
        U = zeros(ndofs*nDoms,1);
        for iiplot=iplot'
            U=extV(:,iiplot);
            figure('name',['Eigenmode alpha = ' num2str(omega(iiplot),'%.3f')])
            vars = {real(U(ext_idx.rho_j)) ,'$\rho$'; 
                    real(U(ext_idx.u_j  )) ,'$u$'; 
                    real(U(ext_idx.v_j  )) ,'$v$'; 
                    real(U(ext_idx.w_j  )) ,'$w$'; 
                    real(U(ext_idx.T_j  )) ,'$T$' };
            title(['$\alpha' num2str(omega(iiplot)) '$']);
            axs = plotFlow(extX,extY,vars,3,2,ext_usedInd,101,'linecolor','none');
            for i=1:length(axs)
                axs(i).YLim=[0,4];
            end
        end
    end
end    
    
%% Spatial Stability Analysis
if spatialAna
    for wi = linspace(0,0.03,10)
    [Lc,Lw,R,idx]   = GetSpatialLinProblem(mesh,baseFlow,'2D',floquetExp);

    % Enforce Dirichlet bc on the top, right and left boundaries, for u,v,w
    % and T
    borders='tb';  vars = 'uvwT';
%     [~,idx_dirchlet] = BC_Dirichlet(L,idx,borders,vars);
    % use most unstable mode from temporal analysis for validation. Spatial
    % analysis should yieald a mode with alpha corresponding to the input
    % of the temporal analysis
    w = real(omega(iplot(1)))+wi*1i; 
    
    alpha_target = alpha;   % target alpha (around which we look for modes)
    L = Lc+Lw*w ;
    [L,R,~] = BC_Dirichlet_SpatialStability(L,R,idx,borders,vars);

    alphas = eigs(L,R,nEig,alpha_target,opts);


    if verbose
        figure('name','Spatial spectra')
            plot(real(alphas),imag(alphas),'ob');
            hold on;
            plot(real(alpha),imag(alpha),'r^','markerSize',15);
            xlabel('$\alpha_r$');
            ylabel('$\alpha_i$');
            legend('Spatial spectra','target');
            grid on;
            title(sprintf('wi=%0.4f',wi))
    end
end
end