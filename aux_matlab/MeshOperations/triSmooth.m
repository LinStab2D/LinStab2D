% TRISMOOTH
%   q = TriSmooth(q,dim,passes) uses a triangular averaging stencile to
%   smooth data in the two-dimentional array q along the dimention specified
%   by dim.  The input 'passes' spefifies the number of times to apply the
%   smoothing function.


function q = TriSmooth(q,dim,passes)

% Define smooting stencile
sten_w = [1/4,1/2,1/4];
sten_i = [-1,0,1];

% If smoothing desired on second dimension, transpose array
if dim == 2
    q = transpose(q);
end

% Length of segments
N = size(q,1);


% Apply #passes of smoothing
for p = 1:passes  
    
    % Set temp array to previous array
    q_tmp = q;
    
    % Apply smoothing to interior points
    for j=2:N-1
        q(j,:) = sten_w*q_tmp(sten_i+j,:);
    end
    
    % Apply smooting to end points (assuming end-point value outside of
    % domain)
    q(1,:)=((3/4).*q_tmp(1,:)+(1/4).*q_tmp(2,:));
    q(N,:)=((3/4).*q_tmp(N,:)+(1/4).*q_tmp(N-1,:));
    
end


% If smoothing on second dimension, transpose array back to original form
if dim == 2
    q = transpose(q);
end


end

