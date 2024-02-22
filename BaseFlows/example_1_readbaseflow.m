function baseFlow = example_1_readbaseflow(X,Y,baseFlow)
    
    [U,V,T,RHO] ...
                = blasius(X,Y,baseFlow.Ma,baseFlow.Pr,1/baseFlow.Re,baseFlow.kappa);
    baseFlow.RHO    = RHO; 
    baseFlow.U      = U;
    baseFlow.V      = V; 
    baseFlow.W      = zeros(size(X));  
    baseFlow.T      = T;
    end

%%
function [uue, vue, TTe, rhorhoe] = blasius(X, Y, Me, Pr, nu, gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Blasius numerical solution from White's 'Viscous Fluid Flow'
% page 233 table 4-1. But see section 7-3.1 in White for
% *compressible* laminar BL. Momentum equation same as incompressible
% case, just with different similarity variables. Also cf. Anderson's
% 'Hypersonic and high-temperature gas dynamics', Cantwell ch 8, and
% Schlichting page 254--255.
%
% Written by: Brandon Yeung
% Created: 2021 Sep 08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta_blasius = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 ...
    1.6 1.7 1.8 1.9 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6 3.8 4.0 ...
    4.2 4.4 4.6 4.8 5.0 5.2 5.4 5.6 5.8 6.0]';

f_blasius = [0 0.00235 0.00939 0.02113 0.03755 0.05864 0.08439 0.11474 0.14967 0.18911 0.23299 ...
    0.28121 0.33366 0.39021 0.45072 0.51503 0.58296 0.65430 0.72887 0.80644 0.88680 ...
    1.05495 1.23153 1.41482 1.60328 1.79557 1.99058 2.18747 2.38559 2.58450 2.78388 ...
    2.98355 3.18338 3.38329 3.58325 3.78323 3.98322 4.18322 4.38322 4.58322 4.78322]';

fprime_blasius = [0 0.04696 0.09391 0.14081 0.18761 0.23423 0.28058 0.32653 0.37196 0.41672 0.46063 ...
    0.50354 0.54525 0.58559 0.62439 0.66147 0.69670 0.72993 0.76106 0.79000 0.81669 ...
    0.86330 0.90107 0.93060 0.95288 0.96905 0.98037 0.98797 0.99289 0.99594 0.99777 ...
    0.99882 0.99940 0.99970 0.99986 0.99994 0.999971 0.999981 0.999995 0.999998 0.999999]';

fprimeprime_blasius = [0.46960 0.46956 0.46931 0.46861 0.46725 0.46503 0.46173 0.45718 0.45119 0.44363 0.43438 ...
    0.42337 0.41057 0.39598 0.37969 0.36180 0.34249 0.32195 0.30045 0.27825 0.25567 ...
    0.21058 0.16756 0.12861 0.09511 0.06771 0.04637 0.03054 0.01933 0.01176 0.00687 ...
    0.00386 0.00208 0.00108 0.00054 0.00026 0.000119 0.000052 0.000022 0.000009 0.000003]';

uue         = fprime_blasius;
TTe         = 1 + sqrt(Pr) * (gamma-1)/2 * Me^2 * (1 - uue.^2);
rhorhoe     = 1 ./ TTe;

% Map compressible eta to incompressible eta
eta_incomp  = cumtrapz(eta_blasius, TTe);

% Points to sample
x           = X(:);
y           = Y(:);
eta_interp  = y ./ sqrt(2 * nu * x);

% Interpolate onto sampling grid
uue_interp      = interp1(eta_incomp, uue, eta_interp, 'spline', 1);
TTe_interp      = interp1(eta_incomp, TTe, eta_interp, 'spline', 1);
rhorhoe_interp  = interp1(eta_incomp, rhorhoe, eta_interp, 'spline', 1);

% Return profiles in original shapes
uue     = reshape(uue_interp, size(X));
TTe     = reshape(TTe_interp, size(X));
rhorhoe = reshape(rhorhoe_interp, size(X));


% Compute wall-normal velocity from compressible continuity
rhorhoeuue = rhorhoe .* uue;
ddxrhorhoeuue = zeros(size(X));
for j = 1:size(X,1)
    ddxrhorhoeuue(j,:) = gradient(rhorhoeuue(j,:), X(1,:));
end
rhorhoevue = zeros(size(X));
for i = 1:size(X,2)
    rhorhoevue(:,i) = -1*cumtrapz(Y(:,1), ddxrhorhoeuue(:,i));
end
vue = rhorhoevue./rhorhoe;

end
