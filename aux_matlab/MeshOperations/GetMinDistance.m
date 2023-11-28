function D = GetMinDistance(mesh,indexes)
%GETMINDISTANCE Summary of this function goes here
%   Detailed explanation goes here
    X = mesh.X;
    Y = mesh.Y;
    
    D = ones(size(X))*inf;
    for i = indexes(:)'
       D = min(D, sqrt( (X-X(i)).^2 + (Y-Y(i)).^2));
    end

end

