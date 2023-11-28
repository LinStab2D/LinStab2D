function [nDomains,extX,extY,extV,ext_idx,ext_usedInd] = ...
        ReconstructFloquetSol(X,Y,Vper,idx,usedInd,floquetExp)
    % [nDomains,extX,extY,extV,ext_idx,ext_usedInd] = ...   
    %        ReconstructFloquetSol(X,Y,Vper,idx,usedInd,floquetExp)
    % Reconstruct the solution field from its periodic component and
    % floquet exponent.
    % Inputs 
    %   X           : gridpoints x coordiantes
    %   Y           : gridpoints y coordiantes
    %   Vper        : Matrix containing solution vectors on its columns
    %   nr,nc       : number of rowns and columns of subplots. nr*nc >=
    %                   size(vars,1)
    %   idx         : structure containing the indexes corresponding to the 
    %                   each flow variable.
    %   usedInd     : Used index vector (needed when mesh masks are used)
    %   floquetExp  : Floquet exponent
    % Outputs
    %   nDomains    : Number of domains needed to reconstruct the solution
    %   extX,extY,
    %   extV,ext_idx,
    %   ext_usedInd : Extended versions of the input to be passed to
    %                   PlotFlow to plot the solution field.
    
    fprintf( '--- Reconstructing fields from periodic solution and floquet Exponent...  '); tic();
    %get grid mesh sizes
    [ny,nx] = size(X);
    N  = nx*ny;
    %get number of dofs
    ndofs = size(Vper,1);
    
    %get domain lenght in x
    LX = (max(X(:)) - min(X(:))) * (nx+1)/nx;
    
    %get number of domains needed to plot the floquer solution
    if floquetExp==0
        nDomains = 1;
    else 
        nDomains = round(abs(2*pi/(LX*floquetExp)));
    end
        

    extX    = [];
    extY    = [];

    ext_idx = struct('rho_j',[],'u_j',[],'v_j',[],'w_j',[],'T_j',[]);
    ext_usedInd = [];
    for i=1:nDomains
        %extend the mesh
        extX = [ extX , X + (i-1)*LX];
        extY = [ extY , Y           ];
        %extend the solution indexes
        ext_idx.rho_j   = [ext_idx.rho_j , idx.rho_j + ndofs*(i-1)];
        ext_idx.u_j     = [ext_idx.u_j   , idx.u_j   + ndofs*(i-1)];
        ext_idx.v_j     = [ext_idx.v_j   , idx.v_j   + ndofs*(i-1)];
        ext_idx.w_j     = [ext_idx.w_j   , idx.w_j   + ndofs*(i-1)];
        ext_idx.T_j     = [ext_idx.T_j   , idx.T_j   + ndofs*(i-1)];
        
        %extend used indexes vector
        ext_usedInd = [ext_usedInd ; usedInd + N*(i-1)];
        
    end
    
    extV = repmat(Vper,nDomains);
    
    % create floquet exponents
    expFloquet = exp(1i*floquetExp*extX(:)); 
    % expand the exponents to all 5 variables, for all periods
    expFloquet = repmat(expFloquet,5*nDomains,1);
    % make it a diagonal matrix
    expFloquet = spdiags(expFloquet,0,N*5*nDomains,N*5*nDomains);
    %Apply the exponent to all vectors
    extV = expFloquet*extV;
    
        
    fprintf( ' Done in %.0f seconds.\n',toc);

