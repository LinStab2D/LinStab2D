function mesh = CreateMesh(xrange,yrange,Nx,Ny,FDorder,y_symmetry,x_periodicity,alpha_filter)
    % mesh = CreateMesh(xrange,yrange,Nx,Ny,FDorder,useSymmetry,periodic,alpha_filter)
    % Creates uniform rectangular mesh with finite difference 
    % differentiation schemes with order FDorder are created. If alpha_filter is 
    % provided, a spatial low-pass filter to filter point-to-point
    % oscilations is created. If 'useSymmetry' is true, differentiation 
    % matrices around the y axis are created, and the mesh is defformed as
    % not to have point at the symmetry line (y=0). If 'periodic'=true,
    % periodicity around the x axis is assumed. The points at the right
    % border are removed.
    % Inputs : 
    %   xrange        : domain range in x.
    %   yrange        : domain range in y.
    %   nx            : number of points in x (==1,for '1D' analysis,or >4)
    %   ny            : number of points in y ( > 4 )
    %   FDorder       : finite difference order (2,4 or 8)
    %   y_symmetry    : create diff. matrices assuming periodicity in y
    %   x_periodicity : periodicity in the x axis (true,false)
    %   alpha_filter  : spatial filter coefficient ('none' or >=0 and , 1)
    % output : 
    %   mesh is an object containing the mesh, differentiation
    %   mesh.X      : x coordinates of the gridpoints
    %   mesh.Y      : y coordinates of the gridpoints
    %   mesh.idx    : indexes of mesh borders
    %   mesh.filters: contains the direct (FILTER) and complex transpose
    %   (FILTER_CT) spatial filters
    %   mesh.W      : integration weights
    %   mesh.D*     : first and second order differentiation matrices
    %   mesh.*_sym/asym  : integration weights and differentiation
    %   matrices assuming symmetry/anti-symmetry around the y axis.
    
    fprintf( '--- Creating mesh...  '); tic();

    % Correct domain length to avoid points on the symmetry line
    if (y_symmetry==true && yrange(1)==0)
        % to use symmetry around the bottom y boundary, point at y = 0 is
        % removed
        dy =  yrange(2)/(Ny-.5);
        yrange(1) = yrange(1)+dy/2;       
    end
    
    % Correct domain length to account for periodicity
    if x_periodicity
        xrange(2) = xrange(2)-(xrange(2)-xrange(1)) / (Nx)  ;
    end
    
    
    % Create mesh
    x = linspace(xrange(1),xrange(2),Nx)';
    y = linspace(yrange(1),yrange(2),Ny)';

    [X,Y]  = meshgrid(x,y);
    originalStructure=size(X);    
    
    % Prepare output
    mesh.X                 = X              ;
    mesh.Y                 = Y              ;
    mesh.ngp               = numel(X)       ;
    mesh.usedInd           = (1:numel(X))'  ;
    mesh.y_symmetry        = y_symmetry    ;
    mesh.x_periodicity     = x_periodicity  ;
    mesh.FDorder           = FDorder        ;

    mesh=AddDW2mesh(mesh); % Create dif matrices and int weights
    
    %create low-pass filter
    if ~exist('alpha_filter')
        %no filtering
        mesh.alpha_filter='none'
    else
        mesh.alpha_filter=alpha_filter;        
    end
    mesh.alpga_filter=alpha_filter;
    [filter,filter_ct]  = GetFilter(mesh,alpha_filter);
    mesh.filters.filter    =  filter    ;
    mesh.filters.filter_ct =  filter_ct ;
    
    % get boundary indexes
    idx.bi    = find(Y(mesh.usedInd)==min(Y(:)));
    idx.ti    = find(Y(mesh.usedInd)==max(Y(:)));
    idx.li    = find(X(mesh.usedInd)==min(X(:)));
    idx.ri    = find(X(mesh.usedInd)==max(X(:)));
    mesh.idx  = idx;
    
    fprintf( ' Done in %.0f seconds.\n',toc); 
