function [invA_fun,invA_T_fun] = GetInverseFunction(A,opts)%,tol,toliLU,maxIter,verbose)
    if ~exist('opts','var' ); opts=struct(); end
    if ~isfield(opts,'type'); opts.type = 'lu'; end
        
    if strcmp(opts.type,'builtin')
        disp('Computing inverse function using matlab internal function')
        umf         = matlab.internal.decomposition.builtin.UMFPACKWrapper(A , 0.1, 0.001, true, 0);
        umf_T       = matlab.internal.decomposition.builtin.UMFPACKWrapper(A', 0.1, 0.001, true, 0);
        invA_fun    = @(v) solve(umf, v, false);
        invA_T_fun  = @(v) solve(umf_T, v, false);
    elseif strcmp(opts.type,'lu')
        disp('Computing inverse function using a LU decomposition')
        [LL,UU,pp,qq,rr]    = lu(A);
        LLt = LL';  UUt = UU'; ppt = pp'; qqt = qq'; rrt = rr';
        invA_fun   = @(v) qq *(UU \(LL \(pp *(rr \v))));
        invA_T_fun = @(v) rrt\(ppt*(LLt\(UUt\(qqt*v))));
    elseif strcmp(opts.type,'\')
        disp("Computing inverse function using matlab's /")
        invA_fun   = @(v) A \v ;
        invA_T_fun = @(v) A'\v ;
    elseif strcmp(opts.type,'ilu')
        disp('Computing inverse function using an iLU preconditioned iterative method')
        if ~isfield(opts,'maxIter') ; opts.maxIter   =100    ; end
        if ~isfield(opts,'tol'    ) ; opts.tol       =1e-6   ; end
        if ~isfield(opts,'toliLU' ) ; opts.toliLU    =1e-6   ; end
        if ~isfield(opts,'verbose') ; opts.verbose   =false  ; end
        if ~isfield(opts,'solver' ) ; opts.solver    ='gmres'; end
        if ~isfield(opts,'thresh' ) ; opts.thresh    =.8; end


        % Setup preconditioner for Ar
        n=size(A,1);
        %Create pivoting matrices/indexes
        %%
        Pi = speye(n); pi = 1:n;
        Pj = speye(n); pj = 1:n;

        Pi = speye(n)*0;
        Pj = speye(n)*0;

        pj = symrcm(A ); % symrcm ,  amd,symamd,dissect,colperm,colamd
        pi = symrcm(A ); % symrcm ,  amd,symamd
        for i=1:n
            Pi(i,pj(i)) = 1; 
            Pj(pi(i),i) = 1; 
        end

        % Get iLU decomposition to be used as preconditioner
%         [iLL,iUU]      = ilu(A(pi,pj),struct('type','ilutp','droptol',opts.toliLU,'milu','off'));
%         [iLL,iUU]      = ilu(A(pi,pj),struct('type','ilutp','droptol',opts.toliLU,'milu','off'));
        [iLL,iUU]      = ilu(A(pi,pj),struct('type','ilutp','droptol',opts.toliLU,'thresh',opts.thresh));
        
        %%
        
%         prec     = @(x) Pi'*(iUU \(iLL \(Pj'*(x))));
%         prec_H   = @(x) Pj *(iLL'\(iUU'\(Pi *(x)))); 

        % Function to solve Al*x = y for x
        if opts.verbose
            if strcmp('bicg',opts.solver)
                invA_fun    = @(x)   Pi'*bicg(A(pi,pj)  ,Pj'*x,opts.tol,opts.maxIter,iLL ,iUU );
                invA_T_fun  = @(x)   Pj *bicg(A(pi,pj)' ,Pi *x,opts.tol,opts.maxIter,iUU',iLL');            
            elseif  strcmp('cgs',opts.solver)
                invA_fun    = @(x)   Pi'* cgs(A(pi,pj)  ,Pj'*x,opts.tol,opts.maxIter,iLL ,iUU );
                invA_T_fun  = @(x)   Pj * cgs(A(pi,pj)' ,Pi *x,opts.tol,opts.maxIter,iUU',iLL');
            elseif strcmp('gmres',opts.solver)
                invA_fun    = @(x)   Pi *gmres(A(pi,pj) ,Pj'*x,[],opts.tol,opts.maxIter,iLL ,iUU );
                invA_T_fun  = @(x)   Pj'*gmres(A(pi,pj)',Pi *x,[],opts.tol,opts.maxIter,iUU',iLL');
            end
        else
            if strcmp('bicg',opts.solver)
                invA_fun    = @(x) Pi'*quiet_iterative(@bicg,A(pi,pj)  ,Pj'*x,opts.tol,opts.maxIter,iLL ,iUU );
                invA_T_fun  = @(x) Pj *quiet_iterative(@bicg,A(pi,pj)' ,Pi *x,opts.tol,opts.maxIter,iUU',iLL');            
            elseif  strcmp('cgs',opts.solver)
                invA_fun    = @(x) Pi'*quiet_iterative(@cgs,A(pi,pj)   ,Pj'*x,opts.tol,opts.maxIter,iLL ,iUU);
                invA_T_fun  = @(x) Pj *quiet_iterative(@cgs,A(pi,pj)'  ,Pi *x,opts.tol,opts.maxIter,iUU',iLL');
            elseif strcmp('gmres',opts.solver)
                invA_fun    = @(x) Pi'*quiet_iterative(@gmres,A(pi,pj) ,Pj'*x,[],opts.tol,opts.maxIter,iLL ,iUU );
                invA_T_fun  = @(x) Pj *quiet_iterative(@gmres,A(pi,pj)',Pi *x,[],opts.tol,opts.maxIter,iUU',iLL');
            end
        end
    end