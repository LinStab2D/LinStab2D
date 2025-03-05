function [A,B,index_set] = BC_Neumann_incomp(A,B,idx,borders,variables,D)   
    % L0,index_set] = SetBoundaryConditions(L0,idx,borders,variables)   
    % Imposes Neumann boundary conditions on A at the selected borders
    % for the choosen variables
    % Inputs :
    %   A, B         : the linear operators to be changed
    %   idx         : a structure containing the indexes of each border for
    %                 each variable
    %   borders     : string array listing the borders to which the b.c. is 
    %                 to be applied. For standanrt mesh, 'lrbt' indicate the 
    %                 left, right, bottom and top boundaries
    %   variables   : string array listing the variables to which the b.c. 
    %                 are to te applied. Options, 'ruvwT', for density 
    %                 (rho), the u,v,w velocities and temperature.
    % Outputs: 
    %   A, B          : linear operators with Neumann b.c. on the selected
    %                 dofs: L0(index_set,index_set)= I 
    %   index_set   : ids of the dofs where b.c. were applied.
    %   

%     fprintf( '--- Applying dirichlet bondary conditions...'); tic();

    % Construct index_set, containing the dofs were b.c. will be applied
    index_set = [] ;
    for b=borders
        for v=variables 
            if v=='u'; v='vr'; end
            if v=='v'; v='vtheta'; end
            if v=='w'; v='vz'; end
            index_str = [b,'i_' v];
            index_set = [index_set;idx.(index_str)(:)];
        end
    end

    % apply b.c.
    index_set= unique(index_set(:));
    A(index_set, :) = 0;
    B(index_set, :) = 0;

    for i=index_set'
        if i<=length(A)/4
            j=i;
            A(i, 1:length(A)/4) = D(j,:);
        elseif i<=length(A)/2
            j=i-length(A)/4;
            A(i, length(A)/4+1:length(A)/2) = D(j,:);
        elseif i<=3*length(A)/4
            j=i-length(A)/2;
            A(i, length(A)/2+1:3*length(A)/4) = D(j,:);
        else
            j=i-3*length(A)/4;
            A(i, 3*length(A)/4+1:end) = D(j,:);
        end    
    end
 
% fprintf( ' Done in %.0f seconds.\n',toc); 

