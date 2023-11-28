set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex')
clear variables
close all
clc
addpath(genpath('aux_matlab'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create mesh and obtain differentiation matrices                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain & grid
% Cartesian mesh in computational domain
y_symmetry      = true;     % use symmetry on y coordinate around y=0 (for axysymmetric problems)
x_periodicity   = false;    % use periodic b.c. on x
alpha           = .0;       % spatial filter coefficient
xrange          = [-1 0 ];  % domain range in x
yrange          = [ 0 1 ];  % domain range in y

Nr              = 50;           % # of grid points (radial)
Nz              = 50;           % # of grid points (streamwise)
FDorder         = 4;            % finite difference order of accuracy

% Cartesian mesh in computational domain

cmesh   = CreateMesh(xrange,yrange,Nz,Nr,FDorder, ...     
                     y_symmetry,x_periodicity,alpha); %construct mesh
x       = cmesh.X;           % x,y: Cartesian grid coordinates
y       = cmesh.Y;

% Grid transformation to parabolic sement in physical domain
d       = 0.2;              % wall-normal thickness parameter
x       = x*d-0.5;
z       = (x+1i*y);
X       = -real(z.^2);      % X,Y: Parabolic grid coordinates
Y       = -imag(z.^2);

mesh    =cmesh;
mesh    = DeformMesh(cmesh,X,Y);
x       = mesh.X;           % x,y: Cartesian grid coordinates
y       = mesh.Y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test Deritives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D       = mesh.Dx;
f       =  cos(x) ; % test function
df_an   = -sin(x) ; % test function analytical derivative

df      = x*0;
df(:)   = D*f(:);

figure('name','Derivatives test')
    vars = {f        ,'$f$'                 ;
            df_an    ,'$df_{analytical}$'   ; 
            df       ,'$df_{numerical}$'    ;
            df-df_an ,'$df_{error}$'            }    ; 
    plotFlow(mesh.X,mesh.Y,vars,2,2)
