function [Wen,invWen] = GetCompEnergyNorm(mesh,BF,model)
    % [Wen,invWen] = GetCompEnergyNorm(mesh,BF,model)
    % Computes the compressible energy perturbation norm and its inverse. 
    % Inputs : 
    %   mesh : mesh structure
    %   BF   : Baseflow structure
    %   model: model used, '2D' or 'axysymmetric'
    % Outputs: 
    %   Wen,invEn  : Energy norm matrix, already containing the integration
    %                weights and its inverse.

if strcmp(model,'axy')
    intWeight = mesh.W .*mesh.Y; %int () r dr dz
elseif  strcmp(model,'2D')
    intWeight = mesh.W ; %int () dy dz
else
    error('Invalid model type. "axy" and "2D" available.');
end    

p=mesh.usedInd;
Wen  = ...
    0.5*[
    (BF.T(p) ./(BF.RHO(p)*BF.kappa*BF.Ma^2)).* intWeight(p)
     BF.RHO(p).*intWeight(p)
     BF.RHO(p).*intWeight(p)
     BF.RHO(p).*intWeight(p)
    (BF.RHO(p)./(BF.kappa*(BF.kappa-1)*BF.T(p)*BF.Ma.^2)).*intWeight(p)
    ];

nDOFs = mesh.ngp*5;
invWen  = spdiags(1./Wen,0,nDOFs,nDOFs);
Wen     = spdiags(Wen   ,0,nDOFs,nDOFs);
% invF = 1./F;