function [A,B,idx] = GetLinProblemIncomp(mesh,BF,model,mkx)
    % [LHS,RHS] = GetLHSRHS(mesh,BF,m,Re,sponge)
    % Creates LHS and RHS operators of the LNS equations, such that
    % LHS (-1i omega) q = RHS q - Sponge q , q=[vr,vtheta,vz,p].
    % Inputs ; 
    %       mesh        : mesh object containing the spatial coordiantes,
    %           diferentiation matrices and integration weights.
    %       BF          : object containing the baseflow to be used, 
    %           including a sponge function
    %       model       : axysimmetric       ('axi') model                                     
    %       mkx         : the azymuthal number (m), if model= 'axy'
    %
    % Outputs :
    %       A,B      : matrices representing the linear system  B omega q = A q 
    %       idx     : structure containing the indexes of the system.
    %                   Indexes for the different varibles fields and the
    %                   borders for each variable.
    %   -----------  NOTE ---------------
    % The current definition of the coordinate system does not mach that
    % used to construct the linear operators (files
    % getCoeffsLaminar_Cartesian and getCoeffsLaminar). See these files for
    % the CS used in them. To be compatible with these files, a change of    
    % coordinate system is thus used here.
    % 
    % Coordinate System definition
    %           +r 
    %           |
    %           |
    %           |______+z 
    %          /
    %         /
    %        /
    %       /
    %      +theta (anti clock-wise) 
    
    if model~='axy' 
        error(' Only axisymmetric model is available ')
    end

    %fprintf( '--- Creating linear operators for model %s', model); tic();
    
    
    % Get mesh and scalar diff. matrices.
    % mesh object is in CS2, and is converted to CS1. 
    z = mesh.X(mesh.usedInd);
    Dz = mesh.Dx;
    D2z = mesh.D2x;
    
    r = mesh.Y(mesh.usedInd);
    Dr = mesh.Dy;
    D2r = mesh.D2y;

    Drz = mesh.Dxy;
    Dzr = mesh.Dyx;

        
    nGridPoints = mesh.ngp;
    
    % Get meanflow in vector form.
    % (Velocity components converted from CS2 to CS1)

    Vz   =  BF.U(mesh.usedInd);
    Vr   =  BF.V(mesh.usedInd);
    Vtheta   =  BF.W(mesh.usedInd);

    %aux matrices
    Z = zeros(nGridPoints, 1);
    I = ones(nGridPoints, 1);
    R = reshape(r, nGridPoints, 1);
    
    %% COMPUTE BASEFLOW DERIVATIVES
    tic
    % Base flow derivatives
    dVzdr    = Dr*Vz;     
    dVrdr    = Dr*Vr;     
    dVthetadr    = Dr*Vtheta;     
    dVzdz    = Dz*Vz;     
    dVrdz    = Dz*Vr;     
    dVthetadz    = Dz*Vtheta;     
   

    %% Coefficient matrices
    %get global diff. matrices
    m=mkx;
    Re = BF.Re;
    [DR,D2R,DZ,D2Z,D2RZ,D2ZR] =   CreateDiffMatrices_Axy(mesh,m,'i');

    %build matrices from coefficients
    ZZ = 0*speye(nGridPoints,nGridPoints); % need a zero diagonal sparse matrix... is there a better way?
    II = speye(nGridPoints,nGridPoints);

    %get coefficients of the linear operator
    getCoeffsLaminarIncomp
   
    LHS = ...
        [ ...
        ZZ ZZ ZZ ZZ
        II ZZ ZZ ZZ
        ZZ II ZZ ZZ
        ZZ ZZ II ZZ
        ];

    A0 =    ...
        [   ...
        diag(sparse(A0_11)) diag(sparse(A0_12)) diag(sparse(A0_13)) diag(sparse(A0_14)) 
        diag(sparse(A0_21)) diag(sparse(A0_22)) diag(sparse(A0_23)) diag(sparse(A0_24)) 
        diag(sparse(A0_31)) diag(sparse(A0_32)) diag(sparse(A0_33)) diag(sparse(A0_34)) 
        diag(sparse(A0_41)) diag(sparse(A0_42)) diag(sparse(A0_43)) diag(sparse(A0_44))      
        ];

    Ar =    ...
        [   ...
        diag(sparse(Ar_11)) diag(sparse(Ar_12)) diag(sparse(Ar_13)) diag(sparse(Ar_14)) 
        diag(sparse(Ar_21)) diag(sparse(Ar_22)) diag(sparse(Ar_23)) diag(sparse(Ar_24)) 
        diag(sparse(Ar_31)) diag(sparse(Ar_32)) diag(sparse(Ar_33)) diag(sparse(Ar_34)) 
        diag(sparse(Ar_41)) diag(sparse(Ar_42)) diag(sparse(Ar_43)) diag(sparse(Ar_44)) 
        ];

    Az =    ...
        [   ...
        diag(sparse(Az_11)) diag(sparse(Az_12)) diag(sparse(Az_13)) diag(sparse(Az_14)) 
        diag(sparse(Az_21)) diag(sparse(Az_22)) diag(sparse(Az_23)) diag(sparse(Az_24)) 
        diag(sparse(Az_31)) diag(sparse(Az_32)) diag(sparse(Az_33)) diag(sparse(Az_34)) 
        diag(sparse(Az_41)) diag(sparse(Az_42)) diag(sparse(Az_43)) diag(sparse(Az_44)) 
        ];

    Arr =   ...
        [   ...
        diag(sparse(Arr_11)) diag(sparse(Arr_12)) diag(sparse(Arr_13)) diag(sparse(Arr_14)) 
        diag(sparse(Arr_21)) diag(sparse(Arr_22)) diag(sparse(Arr_23)) diag(sparse(Arr_24)) 
        diag(sparse(Arr_31)) diag(sparse(Arr_32)) diag(sparse(Arr_33)) diag(sparse(Arr_34)) 
        diag(sparse(Arr_41)) diag(sparse(Arr_42)) diag(sparse(Arr_43)) diag(sparse(Arr_44)) 
        ];

    Azz =   ...
        [   ...
        diag(sparse(Azz_11)) diag(sparse(Azz_12)) diag(sparse(Azz_13)) diag(sparse(Azz_14)) 
        diag(sparse(Azz_21)) diag(sparse(Azz_22)) diag(sparse(Azz_23)) diag(sparse(Azz_24)) 
        diag(sparse(Azz_31)) diag(sparse(Azz_32)) diag(sparse(Azz_33)) diag(sparse(Azz_34)) 
        diag(sparse(Azz_41)) diag(sparse(Azz_42)) diag(sparse(Azz_43)) diag(sparse(Azz_44)) 
        ];

    %% Construct the operator
    RHS = A0 + Ar*DR + Az*DZ + Arr*D2R + Azz*D2Z;

    %% Construct the varibles and boundary indexes.
    % Row indices for primitive variables 
    idx = struct();
    ngp        = mesh.ngp;

    idx = mesh.idx;

    boundaries = fields(mesh.idx);
    variables  = {'vr','vtheta','vz','p'}; 
    
    for ib = 1:length(boundaries)
        b  = boundaries{ib};
        idx.([b '_all']) = [];
        for iv = 1:length(variables)
            v = variables{iv};
            idx.([b '_' v]) =  idx.(b) + ngp*(iv-1);
            idx.([b '_all'])=  [idx.([b '_all']); idx.([b '_' v])] ;
        end
    end
   
    % Column indices for conserved variables
    for iv = 1:length(variables)
        v = variables{iv};
        idx.([v '_j']) =  (1:ngp) + ngp*(iv-1);
    end
    
    A = RHS;
    B = LHS.*(-1i);
    
    %% Add sponge function to the linear operator
    if isfield(mesh,'sponge')
        sponge  = mesh.sponge(mesh.usedInd);
        vec_sponge = [zeros(ngp,1); repmat(sponge(:),3,1)];
        Asponge = spdiags(vec_sponge,0,ngp*4,ngp*4);
        A = A - Asponge;   
    end
    
    time    = toc;
    %fprintf( ' Done in %.0f seconds.\n',toc); 

