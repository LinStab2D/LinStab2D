function [L0,index_set] = BC_Neumann(L0,mesh,idx,borders,variables,penalty)   
    % [L0,index_set] = BC_Neumann(L0,idx,borders,variables)   
    % Imposes Neumann boundary conditions on L0 at the selected borders
    % for the choosen variables
    % Inputs :
    %   L0          : the linear operator to be changed
    %   mesh        : mesh object. Used to get boundary normals and
    %                 derivative matrices
    %   idx         : a structure containing the indexes of each border for
    %                 each variable
    %   borders     : string array listing the borders to which the b.c. is 
    %                 to be applied. For standanrt mesh, 'lrbt' indicate the 
    %                 left, right, bottom and top boundaries
    %   variables   : string array listing the variables to which the b.c. 
    %                 are to te applied. Options, 'ruvwT', for density 
    %                 (rho), the u,v,w velocities and temperature.
    %   penalty     : penalty value used to impose bc.
    % Outputs: 
    %   L0          : linear operator with Dirichlet b.c. on the selected
    %                 dofs: L0(index_set,index_set)= I 
    %   index_set   : ids of the dofs where b.c. were applied.
    %   

    if ~exist('penalty'); penalty=1e4; end

    fprintf( '--- Applying Neumann bondary conditions...'); tic();

    [nx,ny] = size(mesh.X);

    % compute mesh gradient
    gradX_x = zeros(size(mesh.X));
    gradX_y = zeros(size(mesh.X));
    gradY_x = zeros(size(mesh.X));
    gradY_y = zeros(size(mesh.X));

    gradX_x(:) = mesh.Dx*mesh.X(:);
    gradX_y(:) = mesh.Dy*mesh.X(:);
    gradY_x(:) = mesh.Dx*mesh.Y(:);
    gradY_y(:) = mesh.Dy*mesh.Y(:);
    
    %normalize vector
    %n_norm = sqrt(nx.^2 + ny.^2);
    %nx = nx./n_norm;
    %ny = ny./n_norm;

    % Construct index_set, containing the dofs were b.c. will be applied
    index_set  = [] ;
    index_mesh = [] ;
    normal_direction  = [] ;
    for b=borders
        for v=variables 
            if v=='r';  v='rho'; end  % if density, expand the string

            index_str = [b,'i_' v];
            index_set = [index_set ;idx.(index_str)(:)];

            index_str = [b,'i'];
            index_mesh= [index_mesh;idx.(index_str)(:)];

            n_ids = length(idx.(index_str)(:))
            if b == 't' || b == 'b'
                normal_direction = [normal_direction; [zeros(n_ids,1) ones(n_ids,1) ] ] ;
            elseif b == 'l' || b == 'r'
                normal_direction = [normal_direction; [ones(n_ids,1)  zeros(n_ids,1)] ] ;
            else
                error("Neumman B.C. are currently only supported for b,t,l and r.")
            end

        end
    end

    % avoid having the same node more than once (e.g., at corners)
    [index_set,unique_pos] = unique(index_set(:)); % avoid repetitions
    index_mesh = index_mesh(unique_pos  );
    normal_direction  = normal_direction (unique_pos,:);

    % apply b.c.
    
    L0(index_set, :) = 0;
    for i = 1:length(index_set(:))
        imesh = mod(index_set(i)-1,nx*ny)+1  ; 
        i_dof = (index_set(i)-imesh)/(nx*ny) ;

        ii = (1:nx*ny) + (i_dof)*nx*ny;
        L0(index_set(i), : ) = 0; 
        L0(index_set(i), ii) =  (                                       ...
                              ( gradX_x(imesh)*mesh.Dx(imesh,:) +       ...
                                gradX_y(imesh)*mesh.Dy(imesh,:)    )*   ... 
                                        normal_direction(i,1) +     ...
                              ( gradY_x(imesh)*mesh.Dx(imesh,:) +       ...
                                gradY_y(imesh)*mesh.Dy(imesh,:)    )*   ...
                                        normal_direction(i,2)       ...
                            )*penalty; 

    end
    

fprintf( ' Done in %.0f seconds.\n',toc); 

