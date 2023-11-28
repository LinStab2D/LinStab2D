function [Lc,Lw,R,idx,L0,R0,R1,R2] = GetSpatialLinProblem(mesh,BF,model,floquetExp)
    % Constructs spatial linear operator
    % (Lc + omega Lw] q' = alpha R q',
    % where q' = [q; alpha q]
    % Inputs : 
    %       mesh        : mesh object, e.g. as obtained with CreateMesh
    %       BF          : Baseflow object
    %       model       : 'axy', for axisymmetric problem, or '2D', for a
    %                     cartesian problem.
    %       flowauetExp : Optinal argument giving the floquet exponent in
    %                     the z direction for analysis using the Floquet
    %                     anzats.
    % Output : 
    %       Lc,Lw,R     : Matrices (2ndofs x 2ndofs) representing the 
    %                     spatial stability problem in the expanded space 
    %                     [q,alpha q].  
    %       idx         : Structure containing the indexes of the system.
    %                     Indexes for the different varibles fields and 
    %                     the borders for each variable.
    %      L0,R0,R1,R2  : Matrices (ndofs x ndofs) corresponding to the 
    %                     quadratic eigenproblem in the non-expanded space: 
    %                     (w L0 - RO - alpha R1 - alpga^2 R2)q = 0
    
    if ~exist('floquetExp','var'); floquetExp=0 ;end

    tic; 
    disp('Computing spatial stability operator')
    if ~strcmp(model,'2D')
        error('Currently spatial stability is only valid for 2d analysis')
    end
    alpha_0  =  0;
    alpha_p1 =  1;
    alpha_m1 = -1;
    
    
    % RHS : R0  + alpha R1 + alpha^2 R2
    % Using solutions for RHS with alpha 0,-1 and +1 to obtain Ri
    [A_0 ,idx] = GetLinProblem(mesh,BF,model,alpha_0 ,floquetExp);
    [A_p1,~]   = GetLinProblem(mesh,BF,model,alpha_p1,floquetExp);
    [A_m1,~]   = GetLinProblem(mesh,BF,model,alpha_m1,floquetExp);
    n = size(A_0,1);
    Z = sparse(n,n);
    I = speye(n);
    
    
    % Matrices representing the DR 
    %           (w L0 - RO - alpha R1 - alpga^2 R2)q = 0
    L0 = speye(n)*-1i;
    R0 =  A_0;
    R1 = (A_p1 - A_m1)/2;
    R2 = (A_p1 - R0 - R1);
    
    % From -iw q = (R0 + alpha R1 + alpha^2 R2), construct
    % [ 0          I ][        q ]         [  I  0  ]  [       q ] 
    % [ omega I-R0   ][ alpha  q ] = alpha [  R1 R2 ]  [ alpha q ]

    Lc  = [ Z , I ;-R0 , Z ] ;
    Lw = [ Z , Z ; I  , Z ]*-1i ;
    
    R = [I,Z;R1,R2];
    
    %update idx to include extended space
    for f = fields(idx)'
        ids = idx.(f{1});
        idx.(f{1}) = [ids,ids+n];
    end
    
    disp(['    elapsed time - Spatial stability matrices:',datestr(toc/24/3600, 'HH:MM:SS')]);

    