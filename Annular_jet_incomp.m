set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex')
clear variables
close all
clc
addpath(genpath('aux_matlab'));
addpath(genpath('BaseFlows/vfinal'));

verbose = false;
save_result = true;
compute_adj = false;

Re_list = [2000]; 
S_list = ["0"];

for k=1:length(S_list)
    if S_list(k)=="0"
        m_list = [0,1,2];
    else
        m_list = [0,1,2,-1,-2];
    end
m_list = [-1];
for j=1:length(Re_list)
for i=1:length(m_list)

save_file_name= "Ar01_S"+S_list(k)+"_250x250_m"+num2str(m_list(i))+"_Re"+num2str(Re_list(j))+"_incomp_v3_m-1_HF";

% Physical parameters
baseFlow.Re     = Re_list(j);          % Reynolds number
m               = m_list(i);
St              = 0.2;
omega_target    = St*2*pi + 0.05i; % 0.257644079083815 - 0.000347156751073108i
nEig            = 25;       % Arnoldi method number of eigenvalues
nEig_adj        = 1;

% Load interpolated LES mean flow
path_to_h5 = "Ar01_S"+S_list(k)+"_PINN_1_250x250ext_Rmax5_Zmax5_baseflow_v3.h5";

Vz = h5read(path_to_h5,'//Vz');
Vr = h5read(path_to_h5,'//Vr');
if S_list(k)=="0"
    Vtheta = Vr.*0;
else
    Vtheta = h5read(path_to_h5,'//Vtheta');
end
Vr = Vr';Vtheta = Vtheta';Vz = Vz';
r_1D = h5read(path_to_h5,'//r_1D');
z_1D = h5read(path_to_h5,'//z_1D');

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

% Mask points inside nozzle
x               = cmesh.X;           
y               = cmesh.Y;
mask = zeros(size(x),'logical');
Nz_nozzle = numel(z_1D(z_1D<0));
Nr_nozzle = numel(r_1D(r_1D<0.45));
mask(1:Nr_nozzle,1:Nz_nozzle) = true;
cmesh = MeshMask(cmesh,mask);

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
    p=mesh.usedInd;
    for bondaries= fields(mesh.idx)'
        ids = mesh.idx.(bondaries{1});
        scatter(cmesh.X(p(ids)), cmesh.Y(p(ids)),'filled');
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
        scatter(mesh.X(p(ids)), mesh.Y(p(ids)),'filled');
    end
    legend(fields(mesh.idx),'Location', 'Best')
    xlabel('$x$');
    ylabel('$y$');
    axis equal tight
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Base flow import                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = Vz; %reshape(U,[numel(r_1D) numel(z_1D)]);
V = Vr; %reshape(V,[numel(r_1D) numel(z_1D)]);
W = Vtheta; %reshape(W,[numel(r_1D) numel(z_1D)]);

baseFlow.U      = U;      % baseflow U velocity
baseFlow.V      = V;      % baseflow V velocity
baseFlow.W      = W;      % baseflow W velocity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sponge                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spongeFun = zeros(numel(r_1D),numel(z_1D));

for i=1:1:numel(r_1D)
    inlet_sponge = 0.;
    outlet_sponge = 0;
%     spongeFun(i,end-20:end) = linspace(0,outlet_sponge,21); %linspace(sqrt(inlet_sponge),0,10).^2;
%     spongeFun(i,:) = outlet_sponge/2*(1+tanh(2*(z_1D-8.5)));
end

mesh.sponge = spongeFun;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize base flow                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    figure('name','Base Flow')
    vars = {baseFlow.U(mesh.usedInd),'$V_z$';
            baseFlow.V(mesh.usedInd),'$V_r$';           
            baseFlow.W(mesh.usedInd),'$V_\theta$';
            mesh.sponge(mesh.usedInd),'sponge';       
            mesh.Dy*(baseFlow.U(mesh.usedInd))          ,'$\frac{dV_z}{dr}$';
            mesh.Dx*(baseFlow.U(mesh.usedInd))          ,'$\frac{dV_z}{dz}$';
            mesh.D2y*(baseFlow.V(mesh.usedInd))         ,'$\frac{d^2V_r}{dr^2}$';
            mesh.Dx*(baseFlow.V(mesh.usedInd))          ,'$\frac{dV_r}{dz}$';
            mesh.Dy*(baseFlow.W(mesh.usedInd))          ,'$\frac{dV_\theta}{dr}$';
            mesh.Dx*(baseFlow.W(mesh.usedInd))          ,'$\frac{dV_\theta}{dz}$';
            mesh.D2y*(baseFlow.W(mesh.usedInd))         ,'$\frac{d^2W}{dr^2}$'};
%             reshape(mesh.Dy*(baseFlow.U(:)),Nr,Nz)         ,'$\frac{dU}{dy}$';
%             reshape(mesh.Dx*(baseFlow.U(:)),Nr,Nz)         ,'$\frac{dU}{dx}$';
%             reshape(mesh.Dx*(mesh.Dy*baseFlow.U(:)),Nr,Nz) ,'$\frac{d^2U}{dydx}$';
%             reshape(mesh.D2y*(baseFlow.V(:)),Nr,Nz)         ,'$\frac{dV}{dy}$';
%             reshape(mesh.Dx*(baseFlow.V(:)),Nr,Nz)         ,'$\frac{dV}{dx}$';
%             reshape(mesh.Dy*(baseFlow.W(:)),Nr,Nz)         ,'$\frac{dW}{dy}$';
%             reshape(mesh.Dx*(baseFlow.W(:)),Nr,Nz)         ,'$\frac{dW}{dx}$';
%             reshape(mesh.Dx*(mesh.Dy*baseFlow.W(:)),Nr,Nz) ,'$\frac{d^2W}{dydx}$'};
    plotFlow(mesh.X,mesh.Y,vars,4,3,mesh.usedInd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up linear operator and boundary conditions                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    B dqdt = A q
%%
[A,B,idx]            = GetLinProblemIncomp(mesh,baseFlow,'axy',m);

% Enforce Dirichlet b,c on the top, right and left boundaries, for u,v,w
% and T

% borders             = 'r';  
% vars                = 'uvwp';
% [A,B,idx_neumann]   = BC_Neumann_incomp(A,B,idx,borders,vars,mesh.Dx);

borders             = 't';  
vars                = 'p';
[A,B,idx_neumannp2]   = BC_Neumann_incomp(A,B,idx,borders,vars,mesh.Dy);

borders             = 'lr';  
vars                = 'p';
[A,B,idx_neumannp]   = BC_Neumann_incomp(A,B,idx,borders,vars,mesh.Dx);

% borders             = 'r';  
% vars                = 'w';
% [A,B,idx_neumann]   = BC_Neumann_incomp(A,B,idx,borders,vars,mesh.Dx);

% 
% borders             = 'lt';  
% vars                = 'w';
% [A,B,idx_dirchlet1]   = BC_Dirichlet_incomp(A,B,idx,borders,vars);

borders             = 'ltr';  
vars                = 'uvw';
[A,B,idx_dirchlet]   = BC_Dirichlet_incomp(A,B,idx,borders,vars);

[W,invW]            = GetTurbulentEnergyNorm(mesh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Global modes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setup eigs options
opts.tol    = 1e-8;
opts.disp   = 1;
opts.p      = 200; % It was 300


% compute direct and adjoint modes
[V    ,omega    ] = eigs(A,B,nEig,omega_target,opts);
omega=diag(omega);

if compute_adj
    [V_adj,omega_adj] = eigs(invW'*A'*W',invW'*B'*W',nEig,conj(omega_target),opts);
    omega_adj=diag(omega_adj);
end

%% Save data %%
if save_result
    save(save_file_name)
end

end
end
end
%% PLot data %%
if verbose 
    nplot = 2;
    [~,iplot    ] = sort(real(omega),"descend");     
    iplot     = iplot(1:nplot);

    figure('name','Temporal Spectra');
    plot(real(omega)/(2*pi),imag(omega),'ob')
    if compute_adj
        [~,iplot_adj] = sort(imag(conj(omega_adj)),"descend");
        iplot_adj = iplot_adj(1:nplot);
        hold on
        plot(real(omega_adj)/(2*pi),-imag(omega_adj),'xb')
        hold off
    end

    leg = {'Direct Spectra'};
    hold on
    for iiplot=iplot'
        plot(real(omega(iiplot))/(2*pi),imag(omega(iiplot)),'^','MarkerSize',10);
        leg{end+1} = num2str(omega(iiplot));
    end
    grid on
    xlabel('$St$');
    ylabel('$\omega_i$');
    legend(leg)

    ndofs = size(A,1);
    used_dofs = 1:ndofs;  used_dofs(idx_dirchlet)=[];
    U = zeros(ndofs,1);
    for iiplot=iplot'
        U=V(:,iiplot);
        figure('name',['Direct eigenmode omega = ' num2str(omega(iiplot),'%.3f')])
        vars = {real(U(idx.vr_j  )) ,'$v_r$'; 
                real(U(idx.vtheta_j  )) ,'$v_\theta$'; 
                real(U(idx.vz_j  )) ,'$v_z$';
                real(U(idx.p_j  )) ,'$p$'};
        title(['$\omega = ' num2str(omega(iiplot)) '$, Direct']);
        axs = plotFlow(mesh.X,mesh.Y,vars,4,1,mesh.usedInd);

        figure('name','Sponge effect')
        Vr = reshape(real(U(idx.vr_j  )),Nr,Nz);
        P = reshape(real(U(idx.p_j  )),Nr,Nz);
        plot(z_1D,spongeFun(1,:)/spongeFun(1,end),z_1D,Vr(1,:)/max(Vr(1,:)),z_1D,P(1,:)/max(P(1,:)))
        legend('Sponge','Vr','P')

        figure('name','Continuity of response mode')
        var_r = 1./mesh.Y(:).*U(idx.vr_j,1)+mesh.Dy_symm*(U(idx.vr_j,1))+1j*m*1./mesh.Y(:).*U(idx.vtheta_j,1)+mesh.Dx*(U(idx.vz_j,1));
        var_r = reshape(var_r,Nr,Nz);
        pcolor(mesh.X,mesh.Y,real(var_r)); shading interp
        colorbar;
    end
    if compute_adj
        for iiplot=iplot_adj'
            F = V_adj(:,iiplot);
            U = V(:,iiplot);
            figure('name',['Adjoint eigenmode omega^* = ' num2str(conj(omega_adj(iiplot)),'%.3f')])
            u = abs(sqrt(U(idx.vr_j).^2+U(idx.vtheta_j).^2+U(idx.vz_j).^2));
            f = abs(sqrt(F(idx.vr_j).^2+F(idx.vtheta_j).^2+F(idx.vz_j).^2));
            vars = {u.*f ,'$\Lambda$';
                    real(F(idx.vz_j)),'$F_r$';
                    real(F(idx.vr_j)),'$F_z$';
                    real(F(idx.p_j)),'$Mass$'};
            title(['$\omega^* =' num2str(conj(omega_adj(iiplot))) '$, Adjoint']);
            axs = plotFlow(mesh.X,mesh.Y,vars,4,1,mesh.usedInd);  
        end
    end
end

if verbose
    [~,iplot    ] = sort(real(omega),"descend");
    
    Ntheta = 100;
    theta = linspace(0,2*pi,Ntheta);
    [RR, Theta] = meshgrid(r_1D,theta);
    Xplot = RR.*cos(Theta); Yplot = RR.*sin(Theta);

    z_target = 31/62;
    [~,idx_z] = min(abs(z_1D-z_target));
    
    qstep = 3;
    nplot = 3;
    iplot     = iplot(1:nplot);
    for iiplot=iplot'
        Uz = nan(size(mesh.X)); Ur = nan(size(mesh.X)); Utheta = nan(size(mesh.X));
        Uz(mesh.usedInd) = V(idx.vz_j,iiplot);
        Ur(mesh.usedInd) = V(idx.vr_j,iiplot);
        Utheta(mesh.usedInd)= V(idx.vtheta_j,iiplot);
        Uz_ztarget = repmat(Uz(:,idx_z)',Ntheta,1).*exp(1i*m.*Theta); 
        Ur_ztarget = repmat(Ur(:,idx_z)',Ntheta,1).*exp(1i*m.*Theta);
        Utheta_ztarget = repmat(Utheta(:,idx_z)',Ntheta,1).*exp(1i*m.*Theta);
        Ux = Ur_ztarget.*cos(Theta) - Utheta_ztarget.*sin(Theta);
        Uy = Ur_ztarget.*sin(Theta) + Utheta_ztarget.*cos(Theta);
        
        figure('name',['R-Theta direct eigenmode omega = ' num2str(omega(iiplot),'%.3f')])
        pcolor(Xplot,Yplot,real(Uz_ztarget)); shading interp
%         hold on
%         quiver(Xplot(1:qstep:end,1:qstep:end),Yplot(1:qstep:end,1:qstep:end),real(Ux(1:qstep:end,1:qstep:end)),real(Uy(1:qstep:end,1:qstep:end)),0.35,'k')
%         hold off
        title(['St = ' num2str(real(omega(iiplot))/(2*pi)) ', $\sigma= ' num2str(imag(omega(iiplot))) '$ ,Direct']);
        ylabel({'y'});
        xlabel({'x'});
        colormap(bluewhitered(Nr*Ntheta))
        axis([-0.6 0.6 -0.6 0.6])
    end
end