function [Wen,invWen] = GetTurbulentEnergyNorm(mesh)
    % [Wen,invWen] = GetCompEnergyNorm(mesh,BF,model)
    % Computes the compressible energy perturbation norm and its inverse. 
    % Inputs : 
    %   mesh : mesh structure
    %   BF   : Baseflow structure
    %   model: model used, '2D' or 'axysymmetric'
    % Outputs: 
    %   Wen,invEn  : Energy norm matrix, already containing the integration
    %                weights and its inverse.


intWeight = mesh.W .*mesh.Y; %int () r dr dz

p=mesh.usedInd;
Wen  = ...
    0.5*[
     intWeight(p)
     intWeight(p)
     intWeight(p)
     intWeight(p)*1e-4
    ];

nDOFs = mesh.ngp*4;
invWen  = spdiags(1./Wen,0,nDOFs,nDOFs);
Wen     = spdiags(Wen   ,0,nDOFs,nDOFs);
% invF = 1./F;