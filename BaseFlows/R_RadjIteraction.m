function R_RadjIteraction(f_prev,fList,B,C,W,invW,TM_setup,TM_setup_adj,i_out,n_block,tol,)
            %%
            f = @(t) B*squeeze(f_prev(:,i_m,:))*exp(-2i*pi*fList.'*t);
            q0   = zeros(nq*nq_multi,1);
            
            [q,qhist]    = TM(q0,f ,inf,TM_setup,i_out,n_block,tol);
            %%
            q     = C*q(1:nq,:); 
            q_hat = ifft(q,[],2)*nf*dt_sampling;            
            q_hat = W*q_hat;
            
            %define forcing term for the adjoint problem
            
            % Normalize frequenci components
            norms = zeros(nf,1);
            for i=1:nf
                norms(i)    = norm(q_hat(:,i))   ;
                q_hat(:,i)  = q_hat(:,i)/norms(i);
            end

            
            fadj = @(t) C*q_hat             *exp(2i*pi*fList.'*t);
            [f,fhist]   = TM(q0,fadj ,inf,TM_setup_adj,i_out,n_block,tol);
            f = B*f(1:nq,:); 

            f_hat = fft(f,[],2)*dt_sampling;
            f_hat = invW*f_hat;
            
            % De-normalize frequenci components
            for i=1:nf
                f_hat(:,i)  = f_hat(:,i)*norms(i);
            end
            
            f_out(:,i_m,:,i_iter) = f_hat;
end