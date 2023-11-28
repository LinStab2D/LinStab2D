function BF = example_3_readbaseflow(mesh,BF)
    
    X = mesh.X;
    Y = mesh.Y;
    
    BFdata = load('BaseFlows/BL_BF')  ; %reads 'data' from mat file
    
    BF.RHO      = interp2(BFdata.x,BFdata.y,BFdata.rho,X,Y,'spline');       % Baseflow Density 
    BF.W        = interp2(BFdata.x,BFdata.y,BFdata.U,X,Y,'spline');         % Baseflow Streamwise velocity (normal to the x-y computational plane)
    BF.V        = interp2(BFdata.x,BFdata.y,BFdata.V,X,Y,'spline');         % Baseflow Wall normal velocity
    BF.U        = interp2(BFdata.x,BFdata.y,BFdata.W,X,Y,'spline');         % Baseflow spanwise velocity
    BF.T        = interp2(BFdata.x,BFdata.y,BFdata.T,X,Y,'spline');         % Baseflow Temperature
    BF.MU       = interp2(BFdata.x,BFdata.y,BFdata.mu,X,Y,'spline');        % Baseflow Viscosity
    BF.dMUdT    = interp2(BFdata.x,BFdata.y,BFdata.dmudT,X,Y,'spline');     % Baseflow Viscosity first derivative
    BF.d2MUdT2  = interp2(BFdata.x,BFdata.y,BFdata.d2mudT2,X,Y,'spline');   % Baseflow Viscosity second derivative
    
    BF.cv       = 1/(BF.kappa*(BF.kappa-1)*BF.Ma^2); % ideal gas 
    BF.c1       = (BF.kappa-1)*BF.Re*BF.Pr*BF.Ma^2;  % ideal gas 
    BF.c2       = BF.kappa*BF.Ma^2;
    
    % derivatives of the viscosity can be provided here, or overwriten
    % latter using, e.g., the sutherland_air function.
    BF.kappa    = 1.4;      % ideal gas adiabatic coefficient
