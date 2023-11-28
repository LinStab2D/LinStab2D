function [DR,D2R,DZ,D2Z,D2RZ,D2ZR] = CreateDiffMatrices_Axy(mesh,m)
    % [DR,D2R,DZ,D2Z] = CreateDiffMatrices_Axy(mesh,m)
    % Create global differentiatial matricies to be used in an axysimmetric
    % analysis. 
    tic
    pipeBC=false;

    pipe_idx_bottom     = [];
    pipe_idx_top        = [];

    % [Nr,Nz] = size(mesh.X);
    ngp=mesh.ngp;

    z = mesh.X;
    r = mesh.Y;

    Dz = mesh.Dx;
    Dr = mesh.Dy;
    Dr_symm  = mesh.Dy_symm;       
    Dr_asymm = mesh.Dy_asymm;     

    D2z         = mesh.D2x;
    D2r         = mesh.D2y;
    D2r_symm    = mesh.D2y_symm;       
    D2r_asymm   = mesh.D2y_asymm;     

    D2zr         = mesh.Dxy;
    D2rz         = mesh.Dyx;
    D2zr_symm    = mesh.Dxy_symm;       
    D2rz_symm    = mesh.Dyx_symm;       
    D2rz_asymm   = mesh.Dxy_asymm;       
    D2zr_asymm   = mesh.Dyx_asymm;       


    %%
    calcMethod='';
    Z = 0*speye(ngp,ngp); % need a zero diagonal sparse matrix... is there a better way?

    % r dependent derivatives (depend on symmetry)
    if  m==0
        DR   = blkdiag(Dr_symm,Dr_asymm,Dr_asymm,Dr_symm,Dr_symm);
        D2R  = blkdiag(D2r_symm,D2r_asymm,D2r_asymm,D2r_symm,D2r_symm);
        D2RZ = blkdiag(D2rz_symm,D2rz_asymm,D2rz_asymm,D2rz_symm,D2rz_symm);
        D2ZR = blkdiag(D2zr_symm,D2zr_asymm,D2zr_asymm,D2zr_symm,D2zr_symm);
    elseif m==1
        DR   = blkdiag(Dr_asymm,Dr_symm,Dr_symm,Dr_asymm,Dr_asymm);
        D2R  = blkdiag(D2r_asymm,D2r_symm,D2r_symm,D2r_asymm,D2r_asymm);
        D2RZ = blkdiag(D2rz_asymm,D2rz_symm,D2rz_symm,D2rz_asymm,D2rz_asymm);
        D2ZR = blkdiag(D2zr_asymm,D2zr_symm,D2zr_symm,D2zr_asymm,D2zr_asymm);
    elseif m>=2
        DR   = blkdiag(Dr_asymm,Dr_asymm,Dr_asymm,Dr_asymm,Dr_asymm);
        D2R  = blkdiag(D2r_asymm,D2r_asymm,D2r_asymm,D2r_asymm,D2r_asymm);
        D2RZ = blkdiag(D2rz_asymm,D2rz_asymm,D2rz_asymm,D2rz_asymm,D2rz_asymm);
        D2ZR = blkdiag(D2zr_asymm,D2zr_asymm,D2zr_asymm,D2zr_asymm,D2zr_asymm);
    end

    % z dependent derivatives (do not depend on symmetry)
    DZ      = kron(speye(5,5),Dz);
    D2Z     = kron(speye(5,5),D2z);


