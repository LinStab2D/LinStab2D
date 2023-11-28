function LinProb = SetDirBC(LinProb,DirBCind)
%GETDIRICHETEXPCONTMATRICES Summary of this function goes here
%   Detailed explanation goes here
    nDOFS = size(LinProb.A,1);
    
    T_cont = speye(nDOFS);
    T_exp  = speye(nDOFS);
    
    T_cont(DirBCind,:)=[];
    T_exp (:,DirBCind)=[];
    
    LinProb.T_cont = T_cont;
    LinProb.T_exp  = T_exp;
    
    
end

