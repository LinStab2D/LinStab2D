function [baseFlow] = sutherland_air(baseFlow)
    % function sutherland_air
    % Updates a baseflow viscosity fields based on the temperature field
    % and sutherland's Law.
    % 
    % input  : baseflow
    % output : baseflow with updated viscosity
    fprintf( "--- Setting viscosity using Sutherland's model"); tic();

    kappa   = baseFlow.kappa;
    T_0     = baseFlow.T_0;
    Re      = baseFlow.Re;
    Pr      = baseFlow.Pr;
    Ma      = baseFlow.Ma;

    % Dependant coefficients
    baseFlow.cv 	= 1/(kappa*(kappa-1)*Ma^2);
    baseFlow.c1     = (kappa-1)*Re*Pr*Ma^2;
    baseFlow.c2     = kappa*Ma^2;

    % Sutherland's Law
    S1     = 110.4;                 % Sutherland temperature
    T_ref  = 280.0;                 % reference temperature
    mu_ref = 1.735e-5;              % reference dynamic viscosity

    Tstar  = baseFlow.T*T_0;

    mu_0   = mu_ref*(T_ref+S1)./(T_0  +S1).*(T_0  /T_ref).^(3/2);
    mu     = mu_ref*(T_ref+S1)./(Tstar+S1).*(Tstar/T_ref).^(3/2);
    dmudT  = -mu_ref.*(T_ref+S1)./(Tstar+S1).^2.*(Tstar./T_ref).^(3./2)+3./2.*mu_ref.*(T_ref+S1)./(Tstar+S1).*sqrt(Tstar./T_ref)./T_ref;
    d2mudT2= (2.*mu_ref.*(T_ref+S1)./(Tstar+S1).^3.*(Tstar./T_ref).^(3./2))-3.*mu_ref.*(T_ref+S1)./((Tstar+S1).^2).*sqrt((Tstar./T_ref))./T_ref+3./4.*mu_ref.*(T_ref+S1)./(Tstar+S1).*((Tstar./T_ref).^(-1./2))./(T_ref.^2);

    % non-dimensionalize
    baseFlow.MU     = mu          /mu_0;
    baseFlow.dmudT  = dmudT       /mu_0*T_0;
    baseFlow.d2mudT2= d2mudT2     /mu_0*T_0^2;

    fprintf( ' Done in %.0f seconds.\n',toc); 
