function [LHS,RHS] = SetDirichletBC(LHS,RHS,index_set)   
tic

LHS(index_set, :)   = 0;
RHS(index_set, :)   = 0;

% LHS(index_set, index_set)   = LHS(index_set, index_set) + eye(length(index_set));
% RHS(index_set, index_set)   = RHS(index_set, index_set) + eye(length(index_set));


LHS(index_set, index_set)   = eye(length(index_set));
RHS(index_set, index_set)   = eye(length(index_set));
    
time    = toc;
disp(['    elapsed time - Boundary conditions: ' datestr(time/24/3600, 'HH:MM:SS')]);

