function [L,R,index_set] = BC_Dirichlet_SpatialStability(L,R,idx,borders,variables)   
    % L0,index_set] = SetBoundaryConditions(L0,idx,borders,variables)   
    % Imposes Dirichllet boundary conditions on L0 at the selected borders
    % for the choosen variables
    % Inputs :
    %   L0          : the linear operator to be changed
    %   idx         : a structure containing the indexes of each border for
    %                 each variable
    %   borders     : string array listing the borders to which the b.c. is 
    %                 to be applied. For standanrt mesh, 'lrbt' indicate the 
    %                 left, right, bottom and top boundaries
    %   variables   : string array listing the variables to which the b.c. 
    %                 are to te applied. Options, 'ruvwT', for density 
    %                 (rho), the u,v,w velocities and temperature.
    % Outputs: 
    %   L0          : linear operator with Dirichlet b.c. on the selected
    %                 dofs: L0(index_set,index_set)= I 
    %   index_set   : ids of the dofs where b.c. were applied.
    %   

    fprintf( '--- Applying dirichlet bondary conditions...'); tic();

    % Construct index_set, containing the dofs were b.c. will be applied
    index_set = [] ;
    for b=borders
        for v=variables 
            if v=='r';  v='rho'; end  % if density, expand the string
            index_str = [b,'i_' v];
            index_set = [index_set;idx.(index_str)(:)];
        end
    end

    % apply b.c.
    index_set= unique(index_set(:));
    
    L(index_set,:) = 0;    L(:,index_set) = 0;
    R(index_set,:) = 0;    R(:,index_set) = 0;
    L(index_set,index_set) = eye(length(index_set))*1e4;
    

fprintf( ' Done in %.0f seconds.\n',toc); 

