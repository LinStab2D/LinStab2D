set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex')
clear variables
close all
clc
addpath(genpath('aux_matlab'));
addpath(genpath('BaseFlows'));
verbose         = true;
save_result = false;
Dextend = false;
save_file_name= "Resolvent_Ar01_S02_250x250_m-1_Re1500_incomp_v3";

% Physical parameters
baseFlow.Re     = 2000; % Reynolds number
m               = -1;
St              = 0.04;
omega    = 2*pi*St; % 0.257644079083815 - 0.000347156751073108i
nEig            = 4;       % Arnoldi method number of eigenvalues

% Load interpolated LES mean flow
path_to_h5 = "BaseFlows/Ar01_S02_n20_4500Pa_1stmesh_PINN_2_250x250ext_Rmax5_Zmax5_baseflow_v3.h5";
% Vr = h5read(path_to_h5,'/BaseFlow/Vr');
% Vtheta = h5read(path_to_h5,'/BaseFlow/Vt');
% Vz = h5read(path_to_h5,'/BaseFlow/Vz');
% r_1D = h5read(path_to_h5,'/BaseFlow/R');
% z_1D = h5read(path_to_h5,'/BaseFlow/Z');
Vr = h5read(path_to_h5,'//Vr');
Vtheta = h5read(path_to_h5,'//Vtheta');
Vz = h5read(path_to_h5,'//Vz');
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
alphaFilter     = 0.;      % spatial filter coefficient
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

inlet_sponge = 0.;
outlet_sponge = 0.;
spongeFun(end,:) = inlet_sponge;
for i=1:1:numel(r_1D)
    spongeFun(i,1) = inlet_sponge;
    %spongeFun(i,numel(z_1D)-40:numel(z_1D)) = linspace(0.,outlet_sponge,40+1);
    spongeFun(i,:) = (1+tanh(1*(z_1D-10)))*outlet_sponge/2;
end

figure()
plot(z_1D,spongeFun(1,:))

mesh.sponge = spongeFun;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize base flow                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    figure('name','Base Flow')
    vars = {baseFlow.U(mesh.usedInd),'$U$';
            baseFlow.V(mesh.usedInd),'$V$';           
            baseFlow.W(mesh.usedInd),'$W$';
            mesh.sponge(mesh.usedInd),'sponge';       
            mesh.Dy*(baseFlow.U(mesh.usedInd))          ,'$\frac{dU}{dy}$';
            mesh.Dx*(baseFlow.U(mesh.usedInd))          ,'$\frac{dU}{dx}$';
            mesh.Dx*(mesh.Dy*baseFlow.U(mesh.usedInd))  ,'$\frac{d^2U}{dydx}$';
            mesh.D2y*(baseFlow.V(mesh.usedInd))         ,'$\frac{dV}{dy}$';
            mesh.Dx*(baseFlow.V(mesh.usedInd))          ,'$\frac{dV}{dx}$';
            mesh.Dy*(baseFlow.W(mesh.usedInd))          ,'$\frac{dW}{dy}$';
            mesh.Dx*(baseFlow.W(mesh.usedInd))          ,'$\frac{dW}{dx}$';
            mesh.Dx*(mesh.Dy*baseFlow.W(mesh.usedInd))  ,'$\frac{d^2W}{dydx}$'};
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
%    B omega q = A q
%%

[A,B,idx]            = GetLinProblemIncomp(mesh,baseFlow,'axy',m);

% Enforce Dirichlet b,c on the top, right and left boundaries, for u,v,w
% and T

borders             = 't';  
vars                = 'p';
[A,B,idx_neumannp2]   = BC_Neumann_incomp(A,B,idx,borders,vars,mesh.Dy);
% 
borders             = 'lr';  
vars                = 'p';
[A,B,idx_neumannp]   = BC_Neumann_incomp(A,B,idx,borders,vars,mesh.Dx);

borders             = 'ltr';  
vars                = 'uvw';
[A,B,idx_dirchlet]   = BC_Dirichlet_incomp(A,B,idx,borders,vars);

[W,invW]            = GetTurbulentEnergyNorm(mesh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resolvent analysis                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B dq/dt = A*q + C*v
% u = D*q
% v: input/forcing, u:output/response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Input/forcing matrix
C = mesh.sponge(mesh.usedInd)<=2e-3;
C = [C; C; C; 0*C];
% C = spdiags(C,0,length(C),length(C));
% C(idx.bi_vr,idx.bi_vtheta) = -1i*m*eye(length(idx.bi_vr));
% C(idx.bi_vtheta,idx.bi_vr) = -1/(1i*m)*eye(length(idx.bi_vr));

% Output/observable matrix
D = mesh.sponge(mesh.usedInd)>=0;
D = kron(D(:),ones(4,1));   

% Compute optimal forcings and responses
[S,U,V]             = resolvent_incompressible(A,B,omega,nEig,W,invW,C,D);


%% Plot modes and gains
if verbose
    figure('name','Mode gains')
    bar(S);
    xlabel('mode');
    ylabel('gain');
    title('Resolvent gains')
    
    figure('name','Leading resolvent forcing and response modes')
    vars = {real(V(idx.vr_j,1)) ,'$f_r^{(1)}$'; real(U(idx.vr_j,1)) ,'$u_r^{(1)}$';
            real(V(idx.vtheta_j,1)) ,'$f_\theta^{(1)}$'; real(U(idx.vtheta_j,1)) ,'$u_\theta^{(1)}$';
            real(V(idx.vz_j,1)) ,'$f_z^{(1)}$'; real(U(idx.vz_j,1)) ,'$u_z^{(1)}$';
            real(V(idx.p_j,1)) ,'$f_p^{(1)}$'; real(U(idx.p_j,1)) ,'$p^{(1)}$'};
    plotFlow(mesh.X,mesh.Y,vars,4,2,mesh.usedInd);

    figure('name','Continuity of response and forcing modes')
    var_f = 1./mesh.Y(:).*V(idx.vr_j,1)+mesh.Dy_symm*(V(idx.vr_j,1))+1j*m*1./mesh.Y(:).*V(idx.vtheta_j,1)+mesh.Dx*(V(idx.vz_j,1));
    var_f = reshape(var_f,Nr,Nz);
    var_r = 1./mesh.Y(:).*U(idx.vr_j,1)+mesh.Dy_symm*(U(idx.vr_j,1))+1j*m*1./mesh.Y(:).*U(idx.vtheta_j,1)+mesh.Dx*(U(idx.vz_j,1));
    var_r = reshape(var_r,Nr,Nz);
    axs(1) = subplot(2,1,1);
    pcolor(mesh.X,mesh.Y,real(var_r)); shading interp
    colorbar;
    axs(2) = subplot(2,1,2);
    pcolor(mesh.X,mesh.Y,real(var_f)); shading interp
    colorbar;

    Ur = reshape(U(idx.vr_j,1),Nr,Nz);
    Ut = reshape(U(idx.vtheta_j,1),Nr,Nz);
    Uz = reshape(U(idx.vz_j,1),Nr,Nz);
    Up = reshape(U(idx.p_j,1),Nr,Nz);

    figure()
    plot(z_1D,real(var_r(1,:)),z_1D,imag(var_r(1,:)),z_1D,abs(var_r(1,:)))
    legend('Real','Imag','Abs')
    xlabel('r')

    figure()
    plot(r_1D,real(Ur(:,30)),r_1D,real(Ut(:,30)),r_1D,real(Uz(:,30)),r_1D,real(Up(:,30)))
    legend('$u_r$','$u_\theta$','$u_z$','$p$')
    xlabel('r')

    figure()
    plot(z_1D,real(Ur(1,:)),z_1D,real(Ut(1,:)),z_1D,real(Uz(1,:)),z_1D,real(Up(1,:)))
    legend('$u_r$','$u_\theta$','$u_z$','$p$')
    xlabel('z')

    Ntheta = 100;
    theta = linspace(0,2*pi,Ntheta);
    [RR, Theta] = meshgrid(r_1D,theta);
    Xplot = RR.*cos(Theta); Yplot = RR.*sin(Theta);

    z_target = 61/62;
    [~,idx_z] = min(abs(z_1D-z_target));    
    qstep = 3;
    Uz = nan(size(mesh.X)); Ur = nan(size(mesh.X)); Utheta = nan(size(mesh.X));
    Uz(mesh.usedInd) = U(idx.vz_j,1);
    Ur(mesh.usedInd) = U(idx.vr_j,1);
    Utheta(mesh.usedInd)= U(idx.vtheta_j,1);
    Uz_ztarget = repmat(Uz(:,idx_z)',Ntheta,1).*exp(1i*m.*Theta); 
    Ur_ztarget = repmat(Ur(:,idx_z)',Ntheta,1).*exp(1i*m.*Theta);
    Utheta_ztarget = repmat(Utheta(:,idx_z)',Ntheta,1).*exp(1i*m.*Theta);
    Ux = Ur_ztarget.*cos(Theta) - Utheta_ztarget.*sin(Theta);
    Uy = Ur_ztarget.*sin(Theta) + Utheta_ztarget.*cos(Theta);
    
    figure()
    pcolor(Xplot,Yplot,real(Uz_ztarget)); shading interp
%         hold on
%         quiver(Xplot(1:qstep:end,1:qstep:end),Yplot(1:qstep:end,1:qstep:end),real(Ux(1:qstep:end,1:qstep:end)),real(Uy(1:qstep:end,1:qstep:end)),0.35,'k')
%         hold off
    ylabel({'y'});
    xlabel({'x'});
    colormap(bluewhitered(Nr*Ntheta))
    axis([-0.6 0.6 -0.6 0.6])
   
%     figure('name','Second resolvent forcing and response modes')
%     vars = {real(V(idx.vr_j,2)) ,'$f_r^{(2)}$'; real(U(idx.vr_j,2)) ,'$u_r^{(2)}$';
%             real(V(idx.vtheta_j,2)) ,'$f_\theta^{(2)}$'; real(U(idx.vtheta_j,2)) ,'$u_\theta^{(2)}$';
%             real(V(idx.vz_j,2)) ,'$f_z^{(2)}$'; real(U(idx.vz_j,2)) ,'$u_z^{(2)}$';
%             real(V(idx.p_j,2)) ,'$f_p^{(2)}$'; real(U(idx.p_j,2)) ,'$p^{(2)}$'};
%         
%     plotFlow(mesh.X,mesh.Y,vars,4,2,mesh.usedInd);
% 
%     figure('name','Third resolvent forcing and response modes')
%     vars = {real(V(idx.vr_j,3)) ,'$f_r^{(3)}$'; real(U(idx.vr_j,3)) ,'$u_r^{(3)}$';
%             real(V(idx.vtheta_j,3)) ,'$f_\theta^{(3)}$'; real(U(idx.vtheta_j,3)) ,'$u_\theta^{(3)}$';
%             real(V(idx.vz_j,3)) ,'$f_z^{(3)}$'; real(U(idx.vz_j,3)) ,'$u_z^{(3)}$';
%             real(V(idx.p_j,3)) ,'$f_p^{(3)}$'; real(U(idx.p_j,3)) ,'$p^{(3)}$'};
%         
%     plotFlow(mesh.X,mesh.Y,vars,4,2,mesh.usedInd);
    
end
